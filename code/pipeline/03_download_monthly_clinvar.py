#!/usr/bin/env python3
"""Discover and download one reproducible GRCh38 ClinVar VCF per month."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import os
import shutil
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen

BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/"
VCF_RE = re.compile(r"^clinvar_(\d{8})\.vcf\.gz$")


class Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.hrefs.append(href)


@dataclass(frozen=True)
class Release:
    date: str
    url: str

    @property
    def year(self) -> str:
        return self.date[:4]

    @property
    def month(self) -> str:
        return self.date[:6]

    @property
    def filename(self) -> str:
        return self.url.rsplit("/", 1)[-1]


def get_text(url: str) -> str:
    request = Request(url, headers={"User-Agent": "AG-VUS-research/1.0"})
    with urlopen(request, timeout=90) as response:
        return response.read().decode("utf-8", "replace")


def discover_all(start_year: int, end_year: int) -> list[Release]:
    candidates: dict[str, Release] = {}
    for directory_year in range(start_year, end_year + 1):
        year_url = urljoin(BASE, f"{directory_year}/")
        parser = Links()
        parser.feed(get_text(year_url))
        for href in parser.hrefs:
            match = VCF_RE.match(Path(href).name)
            if not match:
                continue
            date = match.group(1)
            release_year = int(date[:4])
            if not start_year <= release_year <= end_year:
                continue
            release = Release(date=date, url=urljoin(year_url, href))
            candidates[release.date] = release
    return sorted(candidates.values(), key=lambda item: item.date)


def select_monthly(releases: list[Release]) -> list[Release]:
    selected: dict[str, Release] = {}
    for release in releases:
        current = selected.get(release.month)
        if current is None or release.date < current.date:
            selected[release.month] = release
    return sorted(selected.values(), key=lambda item: item.date)


def remote_size(url: str) -> int | None:
    request = Request(url, method="HEAD", headers={"User-Agent": "AG-VUS-research/1.0"})
    with urlopen(request, timeout=90) as response:
        value = response.headers.get("Content-Length")
    return int(value) if value else None


def download(url: str, target: Path) -> None:
    if target.exists() and target.stat().st_size > 0:
        return
    part = target.with_suffix(target.suffix + ".part")
    request = Request(url, headers={"User-Agent": "AG-VUS-research/1.0"})
    with urlopen(request, timeout=300) as response, part.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=8 * 1024 * 1024)
    part.replace(target)


def reuse_monthly_file(root: Path, release: Release, target: Path) -> bool:
    source = root / "data" / "raw" / "clinvar" / "monthly_grch38" / release.year / target.name
    if not source.exists() or target.exists():
        return False
    try:
        os.link(source, target)
    except OSError:
        target.symlink_to(source)
    return True


def md5(path: Path) -> str:
    digest = hashlib.md5()  # nosec: required to verify NCBI-provided checksum
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/path/to/project_root/AG-VUS")
    parser.add_argument("--start-year", type=int, default=2017)
    parser.add_argument("--end-year", type=int, default=datetime.now().year)
    parser.add_argument("--selection", choices=("monthly", "all"), default="monthly")
    parser.add_argument("--manifest-only", action="store_true")
    args = parser.parse_args()

    root = Path(args.root)
    collection = "monthly_grch38" if args.selection == "monthly" else "all_dated_grch38"
    out_dir = root / "data" / "raw" / "clinvar" / collection
    report_dir = root / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    all_releases = discover_all(args.start_year, args.end_year)
    releases = select_monthly(all_releases) if args.selection == "monthly" else all_releases
    manifest = report_dir / f"clinvar_{args.selection}_download_manifest.tsv"
    rows: list[dict[str, object]] = []
    for release in releases:
        size = remote_size(release.url)
        rows.append(
            {
                "release_date": release.date,
                "calendar_month": release.month,
                "url": release.url,
                "remote_bytes": size or "",
                "selection_rule": (
                    "earliest_dated_GRCh38_release_in_calendar_month"
                    if args.selection == "monthly"
                    else "all_available_dated_GRCh38_releases"
                ),
                "discovered_at_utc": datetime.now(timezone.utc).isoformat(),
            }
        )

    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Frozen {len(rows)} monthly releases in {manifest}", flush=True)
    print(f"Estimated bytes: {sum(int(row['remote_bytes'] or 0) for row in rows)}", flush=True)
    if args.manifest_only:
        return

    status_path = report_dir / f"clinvar_{args.selection}_download_status.tsv"
    with status_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["release_date", "vcf", "bytes", "md5", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index, release in enumerate(releases, 1):
            year_dir = out_dir / release.year
            year_dir.mkdir(exist_ok=True)
            vcf = year_dir / release.filename
            checksum_file = year_dir / f"{release.filename}.md5"
            index_file = year_dir / f"{release.filename}.tbi"
            status = ""
            last_error: Exception | None = None
            for attempt in range(1, 4):
                try:
                    if args.selection == "all":
                        reuse_monthly_file(root, release, vcf)
                        reuse_monthly_file(root, release, checksum_file)
                        reuse_monthly_file(root, release, index_file)
                    download(release.url, vcf)
                    download(f"{release.url}.md5", checksum_file)
                    download(f"{release.url}.tbi", index_file)
                    expected = checksum_file.read_text(encoding="utf-8").split()[0]
                    observed = md5(vcf)
                    if observed.lower() != expected.lower():
                        raise RuntimeError(
                            f"MD5 mismatch: expected {expected}, observed {observed}"
                        )
                    status = "verified"
                    break
                except Exception as exc:
                    last_error = exc
                    for path in (vcf, checksum_file, index_file):
                        if path.is_symlink() or path.exists():
                            path.unlink()
                    if attempt < 3:
                        print(
                            f"[{index}/{len(releases)}] {release.date}: "
                            f"retry {attempt}/3 after {type(exc).__name__}",
                            flush=True,
                        )
                        time.sleep(5 * attempt)
            if status != "verified":
                assert last_error is not None
                status = f"failed:{type(last_error).__name__}:{last_error}"
            writer.writerow(
                {
                    "release_date": release.date,
                    "vcf": vcf,
                    "bytes": vcf.stat().st_size if vcf.exists() else 0,
                    "md5": md5(vcf) if status == "verified" else "",
                    "status": status,
                }
            )
            handle.flush()
            print(f"[{index}/{len(releases)}] {release.date}: {status}", flush=True)
            if status != "verified":
                sys.exit(1)


if __name__ == "__main__":
    main()
