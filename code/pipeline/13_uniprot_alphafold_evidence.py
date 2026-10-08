import concurrent.futures,json
from pathlib import Path
import pandas as pd,requests
root=Path(__file__).resolve().parents[1]; out=root/'results/objective2_benchmark'
df=pd.read_parquet(out/'variant_scores_with_external_predictors.parquet'); ids=sorted(set(df.uniprot_id.dropna().astype(str)))
def get(uid):
 r={'uniprot_id':uid}
 try:
  u=requests.get('https://rest.uniprot.org/uniprotkb/'+uid+'.json',timeout=20); r['uniprot_status']=u.status_code
  a=requests.get('https://alphafold.ebi.ac.uk/api/prediction/'+uid,timeout=20); r['alphafold_status']=a.status_code; r['alphafold_available']=a.ok
 except Exception as x:r['error_type']=type(x).__name__
 return r
with concurrent.futures.ThreadPoolExecutor(8) as ex: rows=list(ex.map(get,ids))
pd.DataFrame(rows).to_parquet(out/'uniprot_alphafold_evidence.parquet',index=False)
s={'uniprot_ids':len(ids),'uniprot_resolved':sum(x.get('uniprot_status')==200 for x in rows),'alphafold_available':sum(x.get('alphafold_available',False) for x in rows),'residue_plddt_values':0,'boundary':'Residue pLDDT requires MANE/UniProt coordinate reconciliation.'}
(out/'uniprot_alphafold_evidence_summary.json').write_text(json.dumps(s,indent=2)+'\n'); print(json.dumps(s))
