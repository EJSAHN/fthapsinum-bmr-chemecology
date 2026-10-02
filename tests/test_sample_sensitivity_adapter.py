import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import sample_sensitivity as adapter

ROOT = Path(__file__).resolve().parents[1]


def inputs(tmp_path):
    rng = np.random.default_rng(42)
    groups = [(g, w, n) for w in ['dry', 'wet'] for g, n in
              [('wt', 3 if w == 'wet' else 4), ('bmr6', 4), ('bmr12', 3 if w == 'dry' else 4)]]
    samples = pd.DataFrame([{'sample_id': f'{g}_{w}_{i}', 'genotype': g, 'water': w}
                           for g, w, n in groups for i in range(n)])
    genes = [f'gene{i}' for i in range(30)]
    counts = pd.DataFrame(rng.poisson(50, (30, 22)), columns=samples.sample_id)
    counts.insert(0, 'gene_id', genes)
    sets = pd.DataFrame({'gene_id': genes, 'module': ['aromatic']*10 + ['other']*20})
    mapping = samples.copy()
    mapping['all_primary'] = 10000
    mapping['post_fastp_reads'] = 10000
    mapping['target_unique_rpm_input'] = np.arange(22)*100 + 500
    mapping['target_unique_rpm_post_filter'] = mapping.target_unique_rpm_input * 1.01
    mapping['target_specificity'] = 0.95
    cfg = {'focus_module': 'aromatic', 'exclude_samples': ['bmr6_dry_0'],
           'expected_primary_n': 22, 'expected_reduced_n': 21,
           'expected_genes': 30, 'expected_focus_genes': 10}
    for name, table in [('counts', counts), ('samples', samples), ('gene_sets', sets), ('mapping', mapping)]:
        table.to_csv(tmp_path/f'{name}.tsv', sep='\t', index=False)
    (tmp_path/'settings.json').write_text(json.dumps(cfg))
    return argparse.Namespace(counts=tmp_path/'counts.tsv', samples=tmp_path/'samples.tsv',
        gene_sets=tmp_path/'gene_sets.tsv', mapping=tmp_path/'mapping.tsv',
        settings=tmp_path/'settings.json', output=tmp_path/'prepared')


def test_fixed_prepare(tmp_path):
    args = inputs(tmp_path)
    adapter.prepare(args)
    samples = adapter.read(args.output/'samples.tsv')
    assert len(samples) == 22
    assert samples.all_primary.eq(10000).all()
    cfg = json.loads((args.output/'settings.json').read_text())
    assert cfg['focus_genes'] == [f'gene{i}' for i in range(10)]
    assert adapter.read(args.output/'species_assignment_qc.tsv').excluded_in_sensitivity.sum() == 1


def test_no_unlisted_exclusion(tmp_path):
    args = inputs(tmp_path)
    cfg = json.loads(args.settings.read_text()); cfg['exclude_samples'] = ['absent']
    args.settings.write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match='exclusions'):
        adapter.prepare(args)


def test_tissue_offset_not_burden_denominator(tmp_path):
    args = inputs(tmp_path)
    mapping = adapter.read(args.mapping); mapping.loc[0, 'all_primary'] = 10100
    adapter.save(mapping, args.mapping)
    with pytest.raises(ValueError, match='Tissue offsets'):
        adapter.prepare(args)


def test_duplicate_count_identifier(tmp_path):
    args = inputs(tmp_path)
    counts = adapter.read(args.counts); counts.loc[1, 'gene_id'] = counts.loc[0, 'gene_id']
    adapter.save(counts, args.counts)
    with pytest.raises(ValueError, match='duplicate'):
        adapter.prepare(args)


def test_requires_refit_completion(tmp_path):
    args = inputs(tmp_path); adapter.prepare(args)
    with pytest.raises(FileNotFoundError, match='not completed'):
        adapter.summarise(argparse.Namespace(prepared=args.output, refit=tmp_path/'models', output=tmp_path/'out'))


def test_complete_target_is_declared():
    text = (ROOT/'workflow/rules/sample_sensitivity.smk').read_text()
    assert 'rule complete_study:' in text and 'rule sample_sensitivity:' in text
    assert 'workflow/rules/sample_sensitivity.smk' in (ROOT/'Snakefile').read_text()


def test_cli_help():
    subprocess.run([sys.executable, str(ROOT/'workflow/scripts/sample_sensitivity.py'), '--help'], check=True, capture_output=True)
