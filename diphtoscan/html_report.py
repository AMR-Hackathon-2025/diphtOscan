#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML report generation for diphtOscan.

This module generates self-contained HTML reports from analysis results
using Jinja2 templates.
"""

import os
from typing import Dict, Any, List
from jinja2 import Environment, PackageLoader, select_autoescape
import pandas as pd


def generate_html_report(
    results: pd.DataFrame,
    summary: Dict[str, Any],
    outdir: str,
    version: str
) -> str:
    """
    Generate a self-contained HTML report from analysis results.

    Parameters
    ----------
    results : pd.DataFrame
        Analysis results DataFrame with strains as rows.
    summary : Dict[str, Any]
        Summary statistics dictionary.
    outdir : str
        Output directory path.
    version : str
        diphtOscan version string.

    Returns
    -------
    str
        Path to the generated HTML report file.
    """
    env = Environment(
        loader=PackageLoader('diphtoscan', 'templates'),
        autoescape=select_autoescape(['html', 'xml'])
    )
    template = env.get_template('report.html')

    # Prepare data for template
    samples_list: List[Dict[str, Any]] = []
    for idx, row in results.iterrows():
        sample_data = {'name': str(idx)}
        sample_data.update(row.to_dict())
        samples_list.append(sample_data)

    # Get column names (excluding the index name)
    columns = list(results.columns)

    # Identify column groups for styling
    qc_columns = [c for c in columns if c.startswith('qc_')]
    amr_columns = ['AMINOGLYCOSIDE', 'MACROLIDE', 'TETRACYCLINE', 'PHENICOL',
                   'SULFONAMIDE', 'TRIMETHOPRIM', 'BETA-LACTAM', 'FLUOROQUINOLONE']
    amr_columns = [c for c in amr_columns if c in columns]

    context = {
        'version': version,
        'summary': summary,
        'samples': samples_list,
        'columns': columns,
        'qc_columns': qc_columns,
        'amr_columns': amr_columns,
        'total_samples': len(results),
    }

    html_content = template.render(**context)

    # Write HTML file
    output_basename = os.path.basename(outdir)
    html_path = os.path.join(outdir, f"{output_basename}.html")

    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    return html_path


def get_species_color(species: str) -> str:
    """
    Get a color for a given species for visualization.

    Parameters
    ----------
    species : str
        Species name.

    Returns
    -------
    str
        Hex color code.
    """
    colors = {
        'C. diphtheriae': '#4CAF50',
        'C. ulcerans': '#2196F3',
        'C. pseudotuberculosis': '#FF9800',
        'C. belfantii': '#9C27B0',
        'C. rouxii': '#E91E63',
        'C. ramonii': '#00BCD4',
        'unknown': '#9E9E9E',
    }
    return colors.get(species, '#9E9E9E')


def get_tox_color(tox_allele: str) -> str:
    """
    Get a color for tox status visualization.

    Parameters
    ----------
    tox_allele : str
        Tox allele value ('-' for negative).

    Returns
    -------
    str
        Hex color code.
    """
    if tox_allele == '-':
        return '#E0E0E0'  # Grey for negative
    return '#F44336'  # Red for positive


def get_amr_color(value: str) -> str:
    """
    Get a color for AMR gene presence visualization.

    Parameters
    ----------
    value : str
        AMR gene value ('-' for absent).

    Returns
    -------
    str
        Hex color code.
    """
    if value == '-':
        return '#E8F5E9'  # Light green for absent
    return '#FFCDD2'  # Light red for present
