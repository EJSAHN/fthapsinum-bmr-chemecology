"""Explicit read denominators and reference-assignment fractions.

No species identity or biological co-infection is inferred from these metrics.
"""
from __future__ import annotations
import math
from collections.abc import Mapping, Sequence
from typing import Any


def count(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f'{label} is not a read count')
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'Missing/non-numeric {label}') from exc
    if not math.isfinite(x) or x < 0 or x != math.floor(x):
        raise ValueError(f'Invalid integer count {label}: {value}')
    return int(x)


def fastp_read_counts(payload: Mapping[str, Any]) -> tuple[int, int]:
    summary = payload['summary']
    before = count(summary['before_filtering']['total_reads'], 'pre_fastp_reads')
    after = count(summary['after_filtering']['total_reads'], 'post_fastp_reads')
    if before <= 0 or after <= 0 or after > before:
        raise ValueError('Require positive post-filter <= pre-filter reads')
    return before, after


def calculate_metrics(row: Mapping[str, Any], target: str, decoys: Sequence[str],
                      fungal_groups: Sequence[str], rpm_basis: str = 'input') -> dict[str, Any]:
    target = target.lower()
    decoys = tuple(x.lower() for x in decoys)
    fungi = tuple(x.lower() for x in fungal_groups)
    if len(set(fungi)) != len(fungi) or len(set(decoys)) != len(decoys):
        raise ValueError('Duplicate reference groups')
    if target not in fungi or target in decoys or not set(decoys).issubset(fungi):
        raise ValueError('Target and comparator group definitions are inconsistent')
    before = count(row['pre_fastp_reads'], 'pre_fastp_reads')
    after = count(row['post_fastp_reads'], 'post_fastp_reads')
    if not 0 < after <= before:
        raise ValueError('Require positive post-filter <= pre-filter reads')
    amounts = {g: count(row[f'{g}_unique'], f'{g}_unique') for g in fungi}
    if sum(amounts.values()) > after:
        raise ValueError('Fungal unique alignments exceed post-filter reads')
    target_n = amounts[target]
    denom = target_n + sum(amounts[g] for g in decoys)
    total_fungi = sum(amounts.values())
    pre_rpm = 1e6 * target_n / before
    post_rpm = 1e6 * target_n / after
    if rpm_basis not in {'input', 'post_filter'}:
        raise ValueError('rpm_basis must be input or post_filter')
    return {
        'pre_fastp_reads': before,
        'post_fastp_reads': after,
        'target_unique': target_n,
        'target_unique_rpm_input': pre_rpm,
        'target_unique_rpm_post_filter': post_rpm,
        'target_unique_rpm': pre_rpm if rpm_basis == 'input' else post_rpm,
        'rpm_basis': rpm_basis,
        'target_specificity': target_n / denom if denom else math.nan,
        'specificity_comparators': ';'.join(decoys),
        'target_fraction_all_fungal': target_n / total_fungi if total_fungi else math.nan,
    }
