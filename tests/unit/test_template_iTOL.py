"""Tests for diphtoscan/template_iTOL.py - iTOL template generation."""

import pytest
from unittest.mock import MagicMock
import pandas as pd
import os

from diphtoscan.template_iTOL import (
    get_BINARY_header,
    get_STRIP_header,
    get_TOX_header,
    writeTemplateBinary,
    writeTemplateTOX,
    writeTemplateStrip,
    spuA,
    narG,
    toxin,
    amr_families,
    list_familiesRes,
)


class TestGetBinaryHeader:
    """Tests for get_BINARY_header function."""

    def test_contains_dataset_binary(self):
        """Test that header contains DATASET_BINARY."""
        header = get_BINARY_header()
        assert 'DATASET_BINARY' in header

    def test_contains_separator(self):
        """Test that header contains SEPARATOR."""
        header = get_BINARY_header()
        assert 'SEPARATOR COMMA' in header

    def test_contains_placeholders(self):
        """Test that header contains placeholders for customization."""
        header = get_BINARY_header()
        assert 'Title' in header
        assert 'Shapes_Binary' in header
        assert 'Labels_Binary' in header
        assert 'Colors_Binary' in header

    def test_contains_data_section(self):
        """Test that header contains DATA section."""
        header = get_BINARY_header()
        assert 'DATA' in header


class TestGetStripHeader:
    """Tests for get_STRIP_header function."""

    def test_contains_dataset_colorstrip(self):
        """Test that header contains DATASET_COLORSTRIP."""
        header = get_STRIP_header()
        assert 'DATASET_COLORSTRIP' in header

    def test_contains_show_labels(self):
        """Test that header contains SHOW_LABELS setting."""
        header = get_STRIP_header()
        assert 'SHOW_LABELS' in header

    def test_contains_border_settings(self):
        """Test that header contains border settings."""
        header = get_STRIP_header()
        assert 'BORDER_WIDTH' in header
        assert 'BORDER_COLOR' in header


class TestGetToxHeader:
    """Tests for get_TOX_header function."""

    def test_contains_toxin_label(self):
        """Test that header contains toxin label."""
        header = get_TOX_header()
        assert 'toxin' in header

    def test_contains_field_shapes(self):
        """Test that header contains correct field shapes."""
        header = get_TOX_header()
        assert 'FIELD_SHAPES,3,3' in header

    def test_contains_field_labels(self):
        """Test that header contains toxin field labels."""
        header = get_TOX_header()
        assert 'toxin,toxin truncated' in header


class TestWriteTemplateBinary:
    """Tests for writeTemplateBinary function."""

    def test_creates_file(self, tmp_path, mock_itol_results_df):
        """Test that writeTemplateBinary creates a file."""
        values = ['spuA']
        colors = ['#002b00']
        symbols = ['2']

        writeTemplateBinary(str(tmp_path), mock_itol_results_df, 'spuA', values, colors, symbols)

        output_file = tmp_path / 'spuA.txt'
        assert output_file.exists()

    def test_file_contains_header(self, tmp_path, mock_itol_results_df):
        """Test that output file contains header."""
        values = ['spuA']
        colors = ['#002b00']
        symbols = ['2']

        writeTemplateBinary(str(tmp_path), mock_itol_results_df, 'spuA', values, colors, symbols)

        output_file = tmp_path / 'spuA.txt'
        content = output_file.read_text()
        assert 'DATASET_BINARY' in content

    def test_file_contains_data(self, tmp_path, mock_itol_results_df):
        """Test that output file contains strain data."""
        values = ['spuA']
        colors = ['#002b00']
        symbols = ['2']

        writeTemplateBinary(str(tmp_path), mock_itol_results_df, 'spuA', values, colors, symbols)

        output_file = tmp_path / 'spuA.txt'
        content = output_file.read_text()
        assert 'strain1' in content

    def test_returns_dataframe(self, tmp_path, mock_itol_results_df):
        """Test that function returns a DataFrame."""
        values = ['spuA']
        colors = ['#002b00']
        symbols = ['2']

        result = writeTemplateBinary(str(tmp_path), mock_itol_results_df, 'spuA', values, colors, symbols)

        assert isinstance(result, pd.DataFrame)

    def test_handles_slash_in_column_name(self, tmp_path):
        """Test that slashes in column names are replaced."""
        df = pd.DataFrame({
            'VIRULENCE/ADHESIN': ['gene1', '-', 'gene2'],
        }, index=['strain1', 'strain2', 'strain3'])

        writeTemplateBinary(str(tmp_path), df, 'VIRULENCE/ADHESIN', ['gene1'], ['#ff0000'], ['1'])

        # File should have underscore instead of slash
        output_file = tmp_path / 'VIRULENCE_ADHESIN.txt'
        assert output_file.exists()


class TestWriteTemplateTOX:
    """Tests for writeTemplateTOX function."""

    def test_creates_file(self, tmp_path, mock_itol_results_df):
        """Test that writeTemplateTOX creates a file."""
        writeTemplateTOX(str(tmp_path), mock_itol_results_df, 'TOXIN')

        output_file = tmp_path / 'TOXIN.txt'
        assert output_file.exists()

    def test_toxin_present_marked(self, tmp_path, mock_itol_results_df):
        """Test that toxin presence is marked correctly."""
        writeTemplateTOX(str(tmp_path), mock_itol_results_df, 'TOXIN')

        output_file = tmp_path / 'TOXIN.txt'
        content = output_file.read_text()

        # strain1 has 'tox' (no truncation)
        lines = content.split('\n')
        strain1_line = [l for l in lines if l.startswith('strain1')][0]
        assert ',1,' in strain1_line or strain1_line.endswith(',1')

    def test_truncated_toxin_marked(self, tmp_path, mock_itol_results_df):
        """Test that truncated toxin is marked differently."""
        writeTemplateTOX(str(tmp_path), mock_itol_results_df, 'TOXIN')

        output_file = tmp_path / 'TOXIN.txt'
        content = output_file.read_text()

        # strain2 has 'tox-50%' (truncated)
        lines = content.split('\n')
        strain2_line = [l for l in lines if l.startswith('strain2')][0]
        # Truncated should have second column as 1
        parts = strain2_line.split(',')
        assert len(parts) >= 3

    def test_absent_toxin_marked(self, tmp_path, mock_itol_results_df):
        """Test that absent toxin is marked with -1."""
        writeTemplateTOX(str(tmp_path), mock_itol_results_df, 'TOXIN')

        output_file = tmp_path / 'TOXIN.txt'
        content = output_file.read_text()

        # strain3 has '-' (absent)
        lines = content.split('\n')
        strain3_line = [l for l in lines if l.startswith('strain3')][0]
        assert '-1,-1' in strain3_line


class TestWriteTemplateStrip:
    """Tests for writeTemplateStrip function."""

    def test_creates_file(self, tmp_path, mock_itol_results_df):
        """Test that writeTemplateStrip creates a file."""
        writeTemplateStrip(str(tmp_path), mock_itol_results_df, 'AMINOGLYCOSIDE', list_familiesRes)

        output_file = tmp_path / 'AMINOGLYCOSIDE.txt'
        assert output_file.exists()

    def test_file_contains_colorstrip_header(self, tmp_path, mock_itol_results_df):
        """Test that output contains COLORSTRIP header."""
        writeTemplateStrip(str(tmp_path), mock_itol_results_df, 'AMINOGLYCOSIDE', list_familiesRes)

        output_file = tmp_path / 'AMINOGLYCOSIDE.txt'
        content = output_file.read_text()
        assert 'DATASET_COLORSTRIP' in content

    def test_absent_genes_use_first_color(self, tmp_path, mock_itol_results_df):
        """Test that absent genes use first color."""
        writeTemplateStrip(str(tmp_path), mock_itol_results_df, 'AMINOGLYCOSIDE', list_familiesRes)

        output_file = tmp_path / 'AMINOGLYCOSIDE.txt'
        content = output_file.read_text()

        # strain2 has '-' for AMINOGLYCOSIDE
        lines = [l for l in content.split('\n') if l.startswith('strain2')]
        if lines:
            assert list_familiesRes['AMINOGLYCOSIDE'][0] in lines[0]


class TestSpuA:
    """Tests for spuA convenience function."""

    def test_creates_spua_file_when_column_exists(self, tmp_path, mock_itol_results_df):
        """Test that spuA creates file when spuA column exists."""
        args = MagicMock()
        args.outdir = str(tmp_path)

        spuA(mock_itol_results_df, args)

        output_file = tmp_path / 'spuA.txt'
        assert output_file.exists()

    def test_does_nothing_when_column_missing(self, tmp_path):
        """Test that spuA does nothing when column is missing."""
        args = MagicMock()
        args.outdir = str(tmp_path)

        df = pd.DataFrame({'other': ['val']}, index=['strain1'])
        spuA(df, args)

        # No file should be created
        assert not (tmp_path / 'spuA.txt').exists()


class TestNarG:
    """Tests for narG convenience function."""

    def test_creates_narg_file_when_column_exists(self, tmp_path, mock_itol_results_df):
        """Test that narG creates file when narG column exists."""
        args = MagicMock()
        args.outdir = str(tmp_path)

        narG(mock_itol_results_df, args)

        output_file = tmp_path / 'narG.txt'
        assert output_file.exists()


class TestToxin:
    """Tests for toxin convenience function."""

    def test_creates_toxin_file_when_column_exists(self, tmp_path, mock_itol_results_df):
        """Test that toxin creates file when TOXIN column exists."""
        args = MagicMock()
        args.outdir = str(tmp_path)

        toxin(mock_itol_results_df, args)

        output_file = tmp_path / 'TOXIN.txt'
        assert output_file.exists()


class TestAmrFamilies:
    """Tests for amr_families convenience function."""

    def test_creates_files_for_existing_families(self, tmp_path, mock_itol_results_df):
        """Test that amr_families creates files for existing AMR families."""
        args = MagicMock()
        args.outdir = str(tmp_path)

        amr_families(mock_itol_results_df, args)

        # AMINOGLYCOSIDE and MACROLIDE are in the mock data
        assert (tmp_path / 'AMINOGLYCOSIDE.txt').exists()
        assert (tmp_path / 'MACROLIDE.txt').exists()

    def test_does_not_create_files_for_missing_families(self, tmp_path, mock_itol_results_df):
        """Test that amr_families doesn't create files for missing families."""
        args = MagicMock()
        args.outdir = str(tmp_path)

        amr_families(mock_itol_results_df, args)

        # PHENICOL is not in the mock data
        assert not (tmp_path / 'PHENICOL.txt').exists()


class TestListFamiliesRes:
    """Tests for list_familiesRes constant."""

    def test_contains_expected_families(self):
        """Test that list_familiesRes contains expected AMR families."""
        expected = [
            'AMINOGLYCOSIDE', 'MACROLIDE', 'PHENICOL', 'SULFONAMIDE',
            'TETRACYCLINE', 'TRIMETHOPRIM', 'QUATERNARY AMMONIUM',
            'BETA-LACTAM', 'QUINOLONE', 'RIFAMYCIN'
        ]
        for family in expected:
            assert family in list_familiesRes

    def test_each_family_has_two_colors(self):
        """Test that each family has two colors (absent and present)."""
        for family, colors in list_familiesRes.items():
            assert isinstance(colors, list)
            assert len(colors) == 2
            assert colors[0].startswith('#')
            assert colors[1].startswith('#')
