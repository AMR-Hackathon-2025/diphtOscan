#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summary statistics generation for diphtOscan.

This module provides functions to calculate and format summary statistics
from batch analysis results.
"""

from typing import Dict, Any, List
import pandas as pd


def calculate_summary_statistics(results: pd.DataFrame) -> Dict[str, Any]:
    """
    Calculate batch summary statistics from analysis results.

    Parameters
    ----------
    results : pd.DataFrame
        DataFrame containing analysis results with strains as rows.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - total_samples: Number of samples analyzed
        - species_distribution: Dict of species counts (if species column exists)
        - biovar_distribution: Dict of biovar counts (if biovar column exists)
        - tox_positive: Number of tox-positive samples (if tox_allele column exists)
        - tox_percent: Percentage of tox-positive samples
        - amr_prevalence: Dict of AMR family prevalences (if AMR columns exist)
    """
    summary: Dict[str, Any] = {'total_samples': len(results)}

    if len(results) == 0:
        return summary

    # Species distribution
    if 'species' in results.columns:
        species_counts = results['species'].value_counts().to_dict()
        summary['species_distribution'] = species_counts

    # Biovar distribution
    if 'biovar' in results.columns:
        biovar_counts = results['biovar'].value_counts().to_dict()
        # Remove entries with '-' or empty values
        biovar_counts = {k: v for k, v in biovar_counts.items() if k and k != '-'}
        if biovar_counts:
            summary['biovar_distribution'] = biovar_counts

    # ST distribution (top 5)
    if 'ST' in results.columns:
        st_counts = results['ST'].value_counts().head(5).to_dict()
        if st_counts:
            summary['st_distribution'] = st_counts

    # Tox prevalence
    if 'tox_allele' in results.columns:
        tox_positive = len(results[results['tox_allele'] != '-'])
        summary['tox_positive'] = tox_positive
        summary['tox_negative'] = len(results) - tox_positive
        summary['tox_percent'] = round((tox_positive / len(results)) * 100, 1) if len(results) > 0 else 0.0

    # AMR prevalence
    amr_families = ['AMINOGLYCOSIDE', 'MACROLIDE', 'TETRACYCLINE', 'PHENICOL',
                    'SULFONAMIDE', 'TRIMETHOPRIM', 'BETA-LACTAM', 'FLUOROQUINOLONE']
    amr_summary: Dict[str, Dict[str, Any]] = {}

    for family in amr_families:
        if family in results.columns:
            positive = len(results[results[family] != '-'])
            if positive > 0:
                amr_summary[family] = {
                    'count': positive,
                    'percent': round((positive / len(results)) * 100, 1)
                }

    if amr_summary:
        summary['amr_prevalence'] = amr_summary

    # QC summary (if available)
    qc_columns = ['qc_total_length', 'qc_num_contigs', 'qc_n50', 'qc_gc_percent']
    qc_present = [col for col in qc_columns if col in results.columns]
    if qc_present:
        qc_summary = {}
        if 'qc_total_length' in results.columns:
            try:
                lengths = pd.to_numeric(results['qc_total_length'], errors='coerce')
                qc_summary['avg_length'] = int(lengths.mean())
                qc_summary['min_length'] = int(lengths.min())
                qc_summary['max_length'] = int(lengths.max())
            except (ValueError, TypeError):
                pass
        if 'qc_num_contigs' in results.columns:
            try:
                contigs = pd.to_numeric(results['qc_num_contigs'], errors='coerce')
                qc_summary['avg_contigs'] = int(contigs.mean())
            except (ValueError, TypeError):
                pass
        if qc_summary:
            summary['qc_summary'] = qc_summary

    return summary


def format_summary_for_console(summary: Dict[str, Any]) -> str:
    """
    Format summary statistics for console output.

    Parameters
    ----------
    summary : Dict[str, Any]
        Summary statistics dictionary from calculate_summary_statistics().

    Returns
    -------
    str
        Formatted multi-line string for console display.
    """
    lines: List[str] = []
    lines.append("")
    lines.append("=" * 50)
    lines.append("ANALYSIS SUMMARY")
    lines.append("=" * 50)
    lines.append(f"Total samples: {summary.get('total_samples', 0)}")

    # Species distribution
    if 'species_distribution' in summary:
        lines.append("")
        lines.append("Species:")
        for species, count in summary['species_distribution'].items():
            lines.append(f"  {species}: {count}")

    # Biovar distribution
    if 'biovar_distribution' in summary:
        lines.append("")
        lines.append("Biovar:")
        for biovar, count in summary['biovar_distribution'].items():
            lines.append(f"  {biovar}: {count}")

    # ST distribution
    if 'st_distribution' in summary:
        lines.append("")
        lines.append("Top Sequence Types:")
        for st, count in summary['st_distribution'].items():
            lines.append(f"  {st}: {count}")

    # Tox prevalence
    if 'tox_positive' in summary:
        lines.append("")
        lines.append(f"Tox positive: {summary['tox_positive']} ({summary.get('tox_percent', 0)}%)")

    # AMR prevalence
    if 'amr_prevalence' in summary:
        lines.append("")
        lines.append("AMR Prevalence:")
        for family, data in summary['amr_prevalence'].items():
            lines.append(f"  {family}: {data['count']} ({data['percent']}%)")

    # QC summary
    if 'qc_summary' in summary:
        lines.append("")
        lines.append("Assembly Quality:")
        qc = summary['qc_summary']
        if 'avg_length' in qc:
            lines.append(f"  Avg length: {qc['avg_length']:,} bp")
        if 'avg_contigs' in qc:
            lines.append(f"  Avg contigs: {qc['avg_contigs']}")

    lines.append("=" * 50)
    return "\n".join(lines)


def add_summary_to_json(json_output: Dict[str, Any], summary: Dict[str, Any]) -> Dict[str, Any]:
    """
    Add summary statistics to JSON output structure.

    Parameters
    ----------
    json_output : Dict[str, Any]
        Existing JSON output dictionary.
    summary : Dict[str, Any]
        Summary statistics from calculate_summary_statistics().

    Returns
    -------
    Dict[str, Any]
        Updated JSON output with summary field.
    """
    json_output['summary'] = summary
    return json_output
