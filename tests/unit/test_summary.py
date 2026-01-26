#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Unit tests for diphtoscan.summary module."""

import pytest
import pandas as pd
from diphtoscan.summary import (
    calculate_summary_statistics,
    format_summary_for_console,
    add_summary_to_json
)


class TestCalculateSummaryStatistics:
    """Tests for calculate_summary_statistics function."""

    def test_empty_dataframe(self):
        """Empty DataFrame should return minimal summary."""
        df = pd.DataFrame()
        summary = calculate_summary_statistics(df)
        assert summary['total_samples'] == 0

    def test_total_samples_count(self):
        """Should correctly count total samples."""
        df = pd.DataFrame({'species': ['A', 'B', 'C']})
        summary = calculate_summary_statistics(df)
        assert summary['total_samples'] == 3

    def test_species_distribution(self):
        """Should calculate species distribution correctly."""
        df = pd.DataFrame({
            'species': ['C. diphtheriae', 'C. diphtheriae', 'C. ulcerans']
        })
        summary = calculate_summary_statistics(df)
        assert 'species_distribution' in summary
        assert summary['species_distribution']['C. diphtheriae'] == 2
        assert summary['species_distribution']['C. ulcerans'] == 1

    def test_biovar_distribution(self):
        """Should calculate biovar distribution correctly."""
        df = pd.DataFrame({
            'biovar': ['gravis', 'mitis', 'gravis', '-']
        })
        summary = calculate_summary_statistics(df)
        assert 'biovar_distribution' in summary
        assert summary['biovar_distribution']['gravis'] == 2
        assert summary['biovar_distribution']['mitis'] == 1
        assert '-' not in summary['biovar_distribution']

    def test_tox_prevalence(self):
        """Should calculate tox prevalence correctly."""
        df = pd.DataFrame({
            'tox_allele': ['1', '-', '2', '-', '-']
        })
        summary = calculate_summary_statistics(df)
        assert summary['tox_positive'] == 2
        assert summary['tox_negative'] == 3
        assert summary['tox_percent'] == 40.0

    def test_tox_prevalence_all_negative(self):
        """Should handle all tox-negative samples."""
        df = pd.DataFrame({
            'tox_allele': ['-', '-', '-']
        })
        summary = calculate_summary_statistics(df)
        assert summary['tox_positive'] == 0
        assert summary['tox_percent'] == 0.0

    def test_amr_prevalence(self):
        """Should calculate AMR prevalence correctly."""
        df = pd.DataFrame({
            'AMINOGLYCOSIDE': ['aph(3\')', '-', 'aph(3\')'],
            'MACROLIDE': ['-', '-', 'ermX'],
            'TETRACYCLINE': ['-', '-', '-']
        })
        summary = calculate_summary_statistics(df)
        assert 'amr_prevalence' in summary
        assert summary['amr_prevalence']['AMINOGLYCOSIDE']['count'] == 2
        assert summary['amr_prevalence']['MACROLIDE']['count'] == 1
        assert 'TETRACYCLINE' not in summary['amr_prevalence']  # All negative

    def test_amr_prevalence_percent(self):
        """AMR prevalence should include percentage."""
        df = pd.DataFrame({
            'MACROLIDE': ['ermX', '-', 'ermX', 'ermB']
        })
        summary = calculate_summary_statistics(df)
        assert summary['amr_prevalence']['MACROLIDE']['percent'] == 75.0

    def test_qc_summary(self):
        """Should calculate QC summary statistics."""
        df = pd.DataFrame({
            'qc_total_length': ['2500000', '2400000', '2600000'],
            'qc_num_contigs': ['50', '40', '60']
        })
        summary = calculate_summary_statistics(df)
        assert 'qc_summary' in summary
        assert summary['qc_summary']['avg_length'] == 2500000
        assert summary['qc_summary']['avg_contigs'] == 50

    def test_missing_columns_handled(self):
        """Should handle missing columns gracefully."""
        df = pd.DataFrame({'other_column': [1, 2, 3]})
        summary = calculate_summary_statistics(df)
        assert 'species_distribution' not in summary
        assert 'tox_positive' not in summary
        assert 'amr_prevalence' not in summary


class TestFormatSummaryForConsole:
    """Tests for format_summary_for_console function."""

    def test_returns_string(self):
        """Should return a string."""
        summary = {'total_samples': 5}
        result = format_summary_for_console(summary)
        assert isinstance(result, str)

    def test_includes_total_samples(self):
        """Should include total samples count."""
        summary = {'total_samples': 10}
        result = format_summary_for_console(summary)
        assert 'Total samples: 10' in result

    def test_includes_header(self):
        """Should include ANALYSIS SUMMARY header."""
        summary = {'total_samples': 5}
        result = format_summary_for_console(summary)
        assert 'ANALYSIS SUMMARY' in result

    def test_includes_species_distribution(self):
        """Should include species distribution when present."""
        summary = {
            'total_samples': 2,
            'species_distribution': {'C. diphtheriae': 2}
        }
        result = format_summary_for_console(summary)
        assert 'Species:' in result
        assert 'C. diphtheriae: 2' in result

    def test_includes_tox_info(self):
        """Should include tox information when present."""
        summary = {
            'total_samples': 5,
            'tox_positive': 2,
            'tox_percent': 40.0
        }
        result = format_summary_for_console(summary)
        assert 'Tox positive: 2 (40.0%)' in result

    def test_includes_amr_prevalence(self):
        """Should include AMR prevalence when present."""
        summary = {
            'total_samples': 10,
            'amr_prevalence': {
                'MACROLIDE': {'count': 3, 'percent': 30.0}
            }
        }
        result = format_summary_for_console(summary)
        assert 'AMR Prevalence:' in result
        assert 'MACROLIDE: 3 (30.0%)' in result


class TestAddSummaryToJson:
    """Tests for add_summary_to_json function."""

    def test_adds_summary_field(self):
        """Should add summary to JSON output."""
        json_output = {'version': '1.0'}
        summary = {'total_samples': 5}
        result = add_summary_to_json(json_output, summary)
        assert 'summary' in result
        assert result['summary'] == summary

    def test_preserves_existing_fields(self):
        """Should preserve existing JSON fields."""
        json_output = {'version': '1.0', 'samples': {}}
        summary = {'total_samples': 5}
        result = add_summary_to_json(json_output, summary)
        assert result['version'] == '1.0'
        assert result['samples'] == {}
