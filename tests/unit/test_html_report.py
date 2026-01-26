#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Unit tests for diphtoscan.html_report module."""

import os
import tempfile
import pytest
import pandas as pd
from diphtoscan.html_report import (
    generate_html_report,
    get_species_color,
    get_tox_color,
    get_amr_color
)


class TestGenerateHtmlReport:
    """Tests for generate_html_report function."""

    @pytest.fixture
    def sample_results(self):
        """Create sample results DataFrame."""
        return pd.DataFrame({
            'species': ['C. diphtheriae', 'C. ulcerans'],
            'ST': ['ST1', 'ST2'],
            'tox_allele': ['1', '-'],
            'biovar': ['gravis', '-']
        }, index=['sample1', 'sample2'])

    @pytest.fixture
    def sample_summary(self):
        """Create sample summary dictionary."""
        return {
            'total_samples': 2,
            'species_distribution': {
                'C. diphtheriae': 1,
                'C. ulcerans': 1
            },
            'tox_positive': 1,
            'tox_negative': 1,
            'tox_percent': 50.0
        }

    def test_creates_html_file(self, sample_results, sample_summary):
        """Should create an HTML file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                sample_results, sample_summary, tmpdir, '1.8.0'
            )
            assert os.path.exists(html_path)
            assert html_path.endswith('.html')

    def test_html_contains_version(self, sample_results, sample_summary):
        """HTML should contain version string."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                sample_results, sample_summary, tmpdir, '1.8.0'
            )
            with open(html_path, 'r') as f:
                content = f.read()
            assert '1.8.0' in content

    def test_html_contains_sample_names(self, sample_results, sample_summary):
        """HTML should contain sample names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                sample_results, sample_summary, tmpdir, '1.8.0'
            )
            with open(html_path, 'r') as f:
                content = f.read()
            assert 'sample1' in content
            assert 'sample2' in content

    def test_html_contains_summary(self, sample_results, sample_summary):
        """HTML should contain summary section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                sample_results, sample_summary, tmpdir, '1.8.0'
            )
            with open(html_path, 'r') as f:
                content = f.read()
            assert 'Summary' in content
            assert 'Total Samples' in content

    def test_html_is_self_contained(self, sample_results, sample_summary):
        """HTML should be self-contained with embedded CSS."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                sample_results, sample_summary, tmpdir, '1.8.0'
            )
            with open(html_path, 'r') as f:
                content = f.read()
            assert '<style>' in content
            assert '</style>' in content

    def test_html_contains_table(self, sample_results, sample_summary):
        """HTML should contain results table."""
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                sample_results, sample_summary, tmpdir, '1.8.0'
            )
            with open(html_path, 'r') as f:
                content = f.read()
            assert '<table' in content
            assert '</table>' in content

    def test_empty_results(self, sample_summary):
        """Should handle empty results DataFrame."""
        empty_df = pd.DataFrame()
        with tempfile.TemporaryDirectory() as tmpdir:
            html_path = generate_html_report(
                empty_df, sample_summary, tmpdir, '1.8.0'
            )
            assert os.path.exists(html_path)


class TestGetSpeciesColor:
    """Tests for get_species_color function."""

    def test_diphtheriae_color(self):
        """C. diphtheriae should return green."""
        color = get_species_color('C. diphtheriae')
        assert color == '#4CAF50'

    def test_ulcerans_color(self):
        """C. ulcerans should return blue."""
        color = get_species_color('C. ulcerans')
        assert color == '#2196F3'

    def test_pseudotuberculosis_color(self):
        """C. pseudotuberculosis should return orange."""
        color = get_species_color('C. pseudotuberculosis')
        assert color == '#FF9800'

    def test_unknown_species_color(self):
        """Unknown species should return grey."""
        color = get_species_color('unknown')
        assert color == '#9E9E9E'

    def test_unrecognized_species_color(self):
        """Unrecognized species should return grey."""
        color = get_species_color('Some Other Species')
        assert color == '#9E9E9E'


class TestGetToxColor:
    """Tests for get_tox_color function."""

    def test_negative_tox_color(self):
        """Tox-negative should return grey."""
        color = get_tox_color('-')
        assert color == '#E0E0E0'

    def test_positive_tox_color(self):
        """Tox-positive should return red."""
        color = get_tox_color('1')
        assert color == '#F44336'

    def test_any_positive_value(self):
        """Any non-dash value should return red."""
        assert get_tox_color('2') == '#F44336'
        assert get_tox_color('allele_X') == '#F44336'


class TestGetAmrColor:
    """Tests for get_amr_color function."""

    def test_absent_amr_color(self):
        """Absent AMR gene should return light green."""
        color = get_amr_color('-')
        assert color == '#E8F5E9'

    def test_present_amr_color(self):
        """Present AMR gene should return light red."""
        color = get_amr_color('ermX')
        assert color == '#FFCDD2'

    def test_any_present_value(self):
        """Any non-dash value should return light red."""
        assert get_amr_color('aph(3\')') == '#FFCDD2'
        assert get_amr_color('gene1;gene2') == '#FFCDD2'
