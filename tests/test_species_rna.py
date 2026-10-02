import unittest, math, json, subprocess, sys, tempfile, os
from pathlib import Path
import numpy as np
import pandas as pd
from rna_metrics import count, fastp_read_counts, calculate_metrics
import score_sensitivity as ss

FUNGI=['ft','fv','fpro','ffuj','mp']

def row():
    return dict(pre_fastp_reads=10000,post_fastp_reads=8000,ft_unique=100,
                fv_unique=10,fpro_unique=50,ffuj_unique=20,mp_unique=20)

class MetricTests(unittest.TestCase):
    def test_counts(self): self.assertEqual(count('3.0','x'),3)
    def test_invalid_counts(self):
        for x in [-1,1.2,math.nan,math.inf,None,True]:
            with self.subTest(x=x),self.assertRaises(ValueError):count(x,'x')
    def test_fastp(self):
        p={'summary':{'before_filtering':{'total_reads':100},'after_filtering':{'total_reads':80}}}
        self.assertEqual(fastp_read_counts(p),(100,80))
    def test_fastp_missing_before(self):
        with self.assertRaises(KeyError):fastp_read_counts({'summary':{'after_filtering':{'total_reads':3}}})
    def test_both_denominators(self):
        x=calculate_metrics(row(),'ft',['fv','fpro','ffuj'],FUNGI)
        self.assertEqual(x['target_unique_rpm_input'],10000)
        self.assertEqual(x['target_unique_rpm_post_filter'],12500)
    def test_comparator_excludes_mp_when_named(self):
        x=calculate_metrics(row(),'ft',['fv','fpro','ffuj'],FUNGI)
        self.assertAlmostEqual(x['target_specificity'],100/180)
        self.assertAlmostEqual(x['target_fraction_all_fungal'],100/200)
    def test_mp_has_same_denominator_definitions(self):
        x=calculate_metrics(row(),'mp',['ft','fv','fpro','ffuj'],FUNGI,'input')
        self.assertEqual(x['target_unique_rpm'],2000)
        self.assertEqual(x['target_unique_rpm_post_filter'],2500)
    def test_explicit_post_basis(self):
        x=calculate_metrics(row(),'ft',['fv','fpro','ffuj'],FUNGI,'post_filter')
        self.assertEqual(x['target_unique_rpm'],12500)
    def test_zero_target_comparators_is_missing_fraction(self):
        r=row();r.update({g+'_unique':0 for g in FUNGI})
        self.assertTrue(math.isnan(calculate_metrics(r,'ft',['fv','fpro','ffuj'],FUNGI)['target_specificity']))
    def test_missing_comparator_not_zero_filled(self):
        r=row();del r['fv_unique']
        with self.assertRaises(KeyError):calculate_metrics(r,'ft',['fv','fpro','ffuj'],FUNGI)
    def test_invalid_denominators(self):
        for before,after in [(0,0),(100,101),(-1,80)]:
            r=row();r.update(pre_fastp_reads=before,post_fastp_reads=after)
            with self.subTest(before=before),self.assertRaises(ValueError):calculate_metrics(r,'ft',[],FUNGI)
    def test_target_not_its_own_decoy(self):
        with self.assertRaises(ValueError):calculate_metrics(row(),'ft',['ft'],FUNGI)
    def test_duplicate_groups(self):
        with self.assertRaises(ValueError):calculate_metrics(row(),'ft',['fv','fv'],FUNGI)
    def test_more_unique_than_reads_rejected(self):
        r=row();r['ft_unique']=9000
        with self.assertRaises(ValueError):calculate_metrics(r,'ft',['fv'],FUNGI)
    def test_production_cli(self):
        script=Path(__file__).resolve().parents[1]/'workflow/scripts/summarize_mapping.py'
        if not script.is_file():script=Path(__file__).resolve().parents[2]/'workflow/scripts/summarize_mapping.py'
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);r=row();r.update(sample_id='synthetic',pathogen='pdb',genotype='wt',water='dry',dai=3,
              all_primary=8000,sb_unique=7000,overall_alignment_pct=90,
              target_breadth_1x=0.1,target_breadth_3x=0.01,target_breadth_5x=0.0)
            pd.DataFrame([r]).to_csv(td/'input.tsv',sep='\t',index=False)
            for target,decoys in [('FT',['fv','fpro','ffuj']),('MP',['ft','fv','fpro','ffuj'])]:
                cmd=[sys.executable,str(script),'--summary-table',str(td/'input.tsv'),
                  '--output-dir',str(td/target),'--target-group',target,'--rpm-basis','input','--specificity-groups',*decoys]
                p=subprocess.run(cmd,capture_output=True,text=True)
                self.assertEqual(p.returncode,0,p.stdout+p.stderr)
                x=pd.read_csv(td/target/'sample_mapping.tsv',sep='\t').iloc[0]
                self.assertAlmostEqual(x.target_unique_rpm_input,r[target.lower()+'_unique']*100)
                self.assertAlmostEqual(x.target_unique_rpm_post_filter,r[target.lower()+'_unique']*125)

class ScoreTests(unittest.TestCase):
    def test_population_sd_and_set_mean(self):
        x=pd.DataFrame([[1,2,4],[2,4,8]],index=['g1','g2'],columns=['s1','s2','s3'])
        score=ss.standardised_score(x,['g1','g2'])
        np.testing.assert_allclose(score,(np.array([1,2,4])-7/3)/np.std([1,2,4],ddof=0))
    def test_zero_variance_not_silently_dropped(self):
        with self.assertRaises(ValueError):ss.standardised_score(pd.DataFrame([[1,1]],index=['g']),['g'])
    def test_missing_gene(self):
        with self.assertRaises(ValueError):ss.standardised_score(pd.DataFrame([[1,2]],index=['g']),['x'])
    def test_duplicate_gene(self):
        with self.assertRaises(ValueError):ss.standardised_score(pd.DataFrame([[1,2]],index=['g']),['g','g'])
    def test_split_count_and_complements(self):
        x=ss.split_null(np.array([0,1,3,4.]),2)
        self.assertEqual(len(x),math.comb(4,2));np.testing.assert_allclose(np.sort(x),np.sort(-x))
    def test_pair_permutation(self):
        t=pd.DataFrame({'group':['bmr6_wet']*2+['wt_wet']*2,'score':[3.,4.,0.,1.]})
        effect,p,n,kind=ss.label_permutation(t,'bmr6_vs_wt_wet')
        self.assertEqual(effect,3);self.assertEqual(n,6);self.assertAlmostEqual(p,2/6)
    def test_interaction_count(self):
        t=pd.DataFrame({'genotype':['wt']*4+['bmr6']*4,'water':['dry','dry','wet','wet']*2,
          'group':['wt_dry']*2+['wt_wet']*2+['bmr6_dry']*2+['bmr6_wet']*2,
          'score':[1.,2.,3.,4.,1.,2.,5.,6.]})
        effect,p,n,kind=ss.label_permutation(t,'bmr6_water_interaction')
        self.assertAlmostEqual(effect,2);self.assertEqual(n,36);self.assertIn('not_interaction_only',kind)
    def test_hc3_se_against_matrix_formula(self):
        rng=np.random.default_rng(17);groups=np.repeat(ss.GROUPS,4)
        t=pd.DataFrame({'group':groups,'score':rng.normal(size=24),'rpm':rng.uniform(10,10000,24)})
        ans=ss.hc3(t,ss.CONTRASTS['bmr6_water_interaction'],'rpm')
        x=np.column_stack([pd.get_dummies(t.group,dtype=float).reindex(columns=ss.GROUPS).to_numpy(),np.log10(t.rpm.to_numpy()+1)])
        inv=np.linalg.inv(x.T@x);beta=inv@x.T@t.score.to_numpy();h=np.einsum('ij,jk,ik->i',x,inv,x)
        residual=t.score.to_numpy()-x@beta
        v=inv@x.T@np.diag((residual/(1-h))**2)@x@inv
        c=np.array([1,-1,0,-1,1,0,0])
        self.assertAlmostEqual(ans['estimate'],float(c@beta),places=10)
        self.assertAlmostEqual(ans['se'],float(np.sqrt(c@v@c)),places=10)
    def test_alignment_uses_identifiers(self):
        s=pd.Series([5.,2.],index=['b','a']);t=pd.DataFrame({'sample_id':['a','b']})
        self.assertEqual(ss.aligned_table(s,t).score.tolist(),[2,5])

if __name__=='__main__':unittest.main()
