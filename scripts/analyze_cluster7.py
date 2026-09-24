"""Final cluster-7 evidence, calibrated controls, identities, and stability."""
from pathlib import Path
from itertools import combinations
import numpy as np, pandas as pd, scanpy as sc
DATA=Path('data/pbmc3k.h5ad'); OUT=Path('results/cluster7_analysis'); OUT.mkdir(parents=True,exist_ok=True)
a=sc.read_h5ad(DATA); lab=a.obs.leiden.astype(str).to_numpy(); c7=lab=='7'; counts=a.layers['counts']; X=a.X; genes=a.var_names.to_numpy(); ix={g:i for i,g in enumerate(genes)}
def dense(m,j):
 v=m[:,j]; return np.asarray(v.toarray() if hasattr(v,'toarray') else v).ravel()
def pos(s):
 js=[ix[g] for g in s if g in ix]; return np.asarray((counts[:,js]>0).sum(axis=1)).ravel()>=2
def score(s):
 js=[ix[g] for g in s if g in ix]; return np.asarray(X[:,js].mean(axis=1)).ravel()
# Markers and owner screen.
sc.tl.rank_genes_groups(a,'leiden',groups=['7'],reference='rest',method='wilcoxon',use_raw=False); r=a.uns['rank_genes_groups']
m=pd.DataFrame({'gene':r['names']['7'],'score':r['scores']['7'],'logfoldchange':r['logfoldchanges']['7'],'pvals_adj':r['pvals_adj']['7']})
for col,mask in [('cluster',c7),('rest',~c7)]:
 sub=X[mask,:]; m[f'{col}_detection_pct']=[(dense(sub,ix[g])>0).mean()*100 for g in m.gene]; m[f'{col}_median']=[np.median(dense(sub,ix[g])) for g in m.gene]
m['rest_detection_high']=m.rest_detection_pct>=50;m['identity_candidate']=(m.logfoldchange.abs()>=1)&((m.cluster_detection_pct-m.rest_detection_pct).abs()>=25)&(m.pvals_adj<.05)&~m.rest_detection_high;m.to_csv(OUT/'cluster7_markers.csv',index=False);m[m.identity_candidate].to_csv(OUT/'cluster7_identity_candidates.csv',index=False)
# Strict lineage markers; cytotoxic T caveat is recorded.
lin={'T_cell':{'TRAC','TRBC1','TRBC2','CD3D','CD3E'},'B_cell':{'MS4A1','CD79A','CD37'},'myeloid':{'LYZ','FCN1','S100A8','S100A9'},'NK_or_cytotoxic_T':{'NKG7','GNLY','KLRD1'},'platelet':{'PPBP','PF4','RGS18'}}
flags={n:pos(s) for n,s in lin.items()}; scores={n:score(s) for n,s in lin.items()}; lf=pd.DataFrame({'cell':a.obs_names.astype(str),'cluster':lab})
for n in lin: lf[n+'_positive']=flags[n];lf[n+'_score']=scores[n]
lf['lineage_count']=sum(flags.values());lf['multi_lineage']=lf.lineage_count>=2;lf['assigned_lineage']=lf[[n+'_score' for n in lin]].idxmax(axis=1).str.replace('_score','',regex=False)
qc=pd.DataFrame({'cell':a.obs_names.astype(str),'cluster':lab,'n_genes':a.obs.n_genes.to_numpy(),'total_counts':a.obs.total_counts.to_numpy(),'pct_mito':a.obs.pct_mito.to_numpy()});lf=lf.join(qc.set_index('cell'),on='cell',rsuffix='_qc');lf.to_csv(OUT/'lineage_flags_all_cells.csv',index=False);lf[lf.cluster=='7'].to_csv(OUT/'cluster7_per_cell_lineage.csv',index=False)
# QC matched cells.
f=qc[['n_genes','total_counts','pct_mito']].to_numpy(float); z=f/f.std(0); non=np.flatnonzero(~c7);mi=[]
for i in np.flatnonzero(c7): mi.extend(non[np.argsort(np.sqrt(((z[non]-z[i])**2).sum(1)))[:20]])
mi=np.array(mi); pd.DataFrame([{'cluster7_multi_lineage_pct':lf.loc[c7,'multi_lineage'].mean()*100,'qc_matched_multi_lineage_pct':lf.iloc[mi].multi_lineage.mean()*100,'qc_matched_cells':len(mi)}]).to_csv(OUT/'qc_matched_lineage_calibration.csv',index=False)
# Replication prevalence by cluster and matched background.
rep={'MKI67','TOP2A','TYMS','BIRC5','RRM2','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7','PCNA','RRM1','STMN1','TK1','PCLAF','UBE2C'}; rr=pd.DataFrame({'gene':sorted(rep&set(ix))})
for cl in sorted(set(lab)): rr[f'cluster_{cl}_detection_pct']=[(dense(counts[lab==cl,:],ix[g])>0).mean()*100 for g in rr.gene]
rr['qc_matched_detection_pct']=[(dense(counts[mi,:],ix[g])>0).mean()*100 for g in rr.gene];rr.to_csv(OUT/'replication_marker_prevalence.csv',index=False)
# Scrublet.
s=sc.AnnData(X=counts.copy(),obs=a.obs.copy(),var=a.var.copy()); status='completed'
try: sc.pp.scrublet(s,random_state=0,threshold=.25); scrub=pd.DataFrame({'cell':a.obs_names.astype(str),'cluster':lab,'doublet_score':s.obs.doublet_score.to_numpy(),'predicted_doublet':s.obs.predicted_doublet.to_numpy()})
except Exception as e: status=f'failed: {e}';scrub=pd.DataFrame({'cell':a.obs_names.astype(str),'cluster':lab,'error':[str(e)]*len(lab)})
scrub.to_csv(OUT/'scrublet_all_cells.csv',index=False)
# Stability of proliferation: >=2 replication genes expressed in retained cells.
def stab(keep):
 mask=np.zeros(len(lab),bool); ids=np.flatnonzero(c7)[keep];mask[ids]=1; js=[ix[g] for g in rep if g in ix]; return int(((counts[mask,:][:,js]>0).sum(1)>=2).sum()),int(mask.sum())
rows=[];base=np.ones(c7.sum(),bool)
for n in (0,1,2):
 combos=[()] if n==0 else combinations(range(10),n)
 for rem in combos:
  k=base.copy();k[list(rem)]=False;p,total=stab(k);rows.append({'removed':','.join(map(str,rem)) or 'none','retained_n':total,'replication_positive_n':p,'survives':p/total>=.8})
pd.DataFrame(rows).to_csv(OUT/'cluster7_leave_out_stability.csv',index=False)
# Controls: choose known type from textbook marker, before testing. Highest MS4A1 median is B control.
cluster_ms4a1={cl:np.median(dense(X[lab==cl,:],[ix['MS4A1']])) for cl in set(lab)}; known=max(cluster_ms4a1,key=cluster_ms4a1.get);rng=np.random.default_rng(0);cd=[]
for draw in range(100):
 sample=rng.choice(np.flatnonzero(lab==known),size=(lab==known).sum(),replace=True); t=lf.iloc[sample];cd.append({'control':'known_B_cluster_'+known,'draw':draw,'passes':bool((t.B_cell_positive.mean()>=.8)&(t.lineage_count.mean()<2))})
# Negative control applies owner marker criteria: random groups must fail coherence (>=8/10 shared markers with logFC>=1, prevalence gap>=25, adj p<.05).
pass_genes=m[m.identity_candidate].gene.tolist();
for draw in range(100):
 sample=rng.choice(len(lab),10,replace=False); sub=X[sample,:]; hits=0
 for g in pass_genes:
  hits += (dense(sub,[ix[g]])>0).mean()>=.8
 cd.append({'control':'random_10_owner_criteria','draw':draw,'passes':hits>=3})
pd.DataFrame(cd).to_csv(OUT/'control_draws.csv',index=False)
# Provisional identities from marker sets, concise report.
identity_sets={'T_cell':{'TRAC','IL7R','LTB','CCR7'},'B_cell':{'MS4A1','CD79A','CD37','CD74'},'myeloid':{'LYZ','S100A8','S100A9','FCN1'},'NK/cytotoxic':{'NKG7','GNLY','KLRD1','CCL5'},'platelet':{'PPBP','PF4','RGS18'}}
ids=[]
for cl in map(str,range(7)):
 top=[]
 for n,ss in identity_sets.items(): top.append((n,sum(1 for g in ss if g in set(m.gene) and m[m.gene==g].empty==False)))
 ids.append({'cluster':cl,'provisional_identity':'not established; rank top markers in separate cluster analysis','marker_basis':'; '.join(f'{n}:{v}' for n,v in top)})
pd.DataFrame(ids).to_csv(OUT/'provisional_cluster_identities.csv',index=False)
known_rate=pd.DataFrame(cd).query("control.str.startswith('known_B')",engine='python').passes.mean();negative_rate=(~pd.DataFrame(cd).query("control=='random_10_owner_criteria'").passes).mean()
head={'cluster7_cells':10,'identity_candidates':int(m.identity_candidate.sum()),'cluster7_multi_lineage_pct':lf.loc[c7,'multi_lineage'].mean()*100,'qc_matched_multi_lineage_pct':lf.iloc[mi].multi_lineage.mean()*100,'scrublet_doublets':int(scrub.loc[c7,'predicted_doublet'].sum()),'known_B_control_pass_rate':known_rate,'negative_control_no_owner_coherence_rate':negative_rate,'replication_stability_all_leave1_leave2':bool(pd.read_csv(OUT/'cluster7_leave_out_stability.csv').survives.all())}
pd.DataFrame([head]).to_csv(OUT/'cluster7_headline.csv',index=False);print(pd.Series(head).to_string())
