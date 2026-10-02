"""Fixed-set score contrasts and reference-assignment exclusion diagnostics."""
from __future__ import annotations
from itertools import combinations
import math
import numpy as np
import pandas as pd
import statsmodels.api as sm

GROUPS = ['wt_dry','bmr6_dry','bmr12_dry','wt_wet','bmr6_wet','bmr12_wet']
CONTRASTS = {
    'bmr6_vs_wt_wet': {'bmr6_wet':1, 'wt_wet':-1},
    'bmr6_water_interaction': {'bmr6_wet':1, 'bmr6_dry':-1, 'wt_wet':-1, 'wt_dry':1},
    'bmr12_vs_bmr6_wet': {'bmr12_wet':1, 'bmr6_wet':-1},
}


def standardised_score(logcpm: pd.DataFrame, genes: list[str]) -> pd.Series:
    if not genes or len(set(genes)) != len(genes):
        raise ValueError('A non-empty unique fixed gene set is required')
    if not set(genes).issubset(logcpm.index):
        raise ValueError('Fixed-set gene absent from expression matrix')
    x = logcpm.loc[genes].to_numpy(dtype=float, copy=True)
    if not np.isfinite(x).all():
        raise ValueError('Non-finite fixed-set expression')
    sd = x.std(axis=1, ddof=0)
    if np.any(sd == 0):
        raise ValueError('Zero-variance fixed-set gene; do not silently drop it')
    return pd.Series(((x - x.mean(axis=1, keepdims=True)) / sd[:,None]).mean(axis=0),
                     index=logcpm.columns, name='score')


def aligned_table(score: pd.Series, samples: pd.DataFrame) -> pd.DataFrame:
    if samples.sample_id.duplicated().any() or score.index.duplicated().any():
        raise ValueError('Duplicate sample identifier')
    if not set(samples.sample_id).issubset(score.index):
        raise ValueError('Missing sample scores')
    t = samples.set_index('sample_id').copy()
    t['score'] = score.reindex(t.index).to_numpy(float)
    if not np.isfinite(t.score).all():
        raise ValueError('Non-finite score')
    return t


def difference(t: pd.DataFrame, weights: dict[str,float]) -> float:
    ans=0.0
    for g,w in weights.items():
        v=t.loc[t.group.eq(g),'score'].to_numpy(float)
        if len(v)==0: raise ValueError('Empty contrast cell: '+g)
        ans += w*float(v.mean())
    return ans


def split_null(v: np.ndarray, n_a: int) -> np.ndarray:
    if not 0 < n_a < len(v): raise ValueError('Invalid split sizes')
    total=float(v.sum()); nb=len(v)-n_a
    return np.asarray([float(v[list(c)].sum())/n_a-(total-float(v[list(c)].sum()))/nb
                       for c in combinations(range(len(v)),n_a)])


def label_permutation(t: pd.DataFrame, contrast: str) -> tuple[float,float,int,str]:
    observed=difference(t,CONTRASTS[contrast])
    if contrast=='bmr6_water_interaction':
        nulls={}
        for geno in ['wt','bmr6']:
            sub=t[t.genotype.eq(geno)]
            nulls[geno]=split_null(sub.score.to_numpy(float),int(sub.water.eq('wet').sum()))
        null=(nulls['bmr6'][:,None]-nulls['wt'][None,:]).ravel()
        kind='within_genotype_water_label_permutation_not_interaction_only_null'
    else:
        a=next(g for g,w in CONTRASTS[contrast].items() if w==1)
        b=next(g for g,w in CONTRASTS[contrast].items() if w==-1)
        x=t.loc[t.group.eq(a),'score'].to_numpy(float)
        y=t.loc[t.group.eq(b),'score'].to_numpy(float)
        null=split_null(np.concatenate([x,y]),len(x))
        kind='exact_two_group_label_permutation'
    p=float(np.mean(np.abs(null)>=abs(observed)-1e-12))
    return observed,p,len(null),kind


def hc3(t: pd.DataFrame, weights: dict[str,float], rpm_column: str | None) -> dict:
    parts=[pd.get_dummies(t.group,dtype=float).reindex(columns=GROUPS,fill_value=0).to_numpy()]
    if rpm_column:
        v=pd.to_numeric(t[rpm_column],errors='raise').to_numpy(float)
        if not np.isfinite(v).all() or np.any(v<0): raise ValueError('Invalid burden proxy')
        parts.append(np.log10(v+1)[:,None])
    x=np.column_stack(parts)
    if np.linalg.matrix_rank(x)!=x.shape[1]: raise ValueError('Rank-deficient model')
    fit=sm.OLS(t.score.to_numpy(float),x).fit(cov_type='HC3',use_t=False)
    c=np.zeros(x.shape[1]);
    for g,w in weights.items(): c[GROUPS.index(g)]=w
    test=fit.t_test(c); ci=np.asarray(test.conf_int()).ravel()
    return {'estimate':float(np.asarray(test.effect).item()),
            'se':float(np.asarray(test.sd).item()),'p':float(np.asarray(test.pvalue).item()),
            'ci_low':float(ci[0]),'ci_high':float(ci[1]),'n':len(t)}


def analyse(score: pd.Series, samples: pd.DataFrame, label: str) -> pd.DataFrame:
    t=aligned_table(score,samples); rows=[]
    for c in CONTRASTS:
        effect,p,n,kind=label_permutation(t,c)
        row={'analysis':label,'contrast':c,'score_difference':effect,
             'label_permutation_p':p,'permutation_assignments':n,'permutation_null':kind}
        for tag,col in [('unadjusted',None),('input_rpm','target_unique_rpm_input'),
                        ('post_filter_rpm','target_unique_rpm_post_filter')]:
            for k,v in hc3(t,CONTRASTS[c],col).items(): row[f'{tag}_{k}']=v
        rows.append(row)
    return pd.DataFrame(rows)
