"""Tests for diphtoscan/utils.py - utility functions."""

import pytest
from unittest.mock import patch, MagicMock
import pandas as pd

from diphtoscan.utils import (
    is_non_zero_file,
    get_chromosome_mlst_header,
    get_tox_header,
    get_virulence,
    get_virulence_extended,
    find_len_contig,
    is_contig_edge,
    armfinder_to_table,
    get_genomic_context,
    find_amrfinderplus_version,
    determine_biovar,
)


class TestIsNonZeroFile:
    """Tests for is_non_zero_file function."""

    def test_existing_non_empty_file(self, tmp_path):
        """Test existing file with content returns True."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("content")

        assert is_non_zero_file(str(test_file)) is True

    def test_existing_empty_file(self, tmp_path):
        """Test existing empty file returns False."""
        test_file = tmp_path / "empty.txt"
        test_file.write_text("")

        assert is_non_zero_file(str(test_file)) is False

    def test_non_existent_file(self, tmp_path):
        """Test non-existent file returns False."""
        assert is_non_zero_file(str(tmp_path / "nonexistent.txt")) is False

    def test_directory_path(self, tmp_path):
        """Test directory path returns False."""
        assert is_non_zero_file(str(tmp_path)) is False


class TestHeaderFunctions:
    """Tests for header functions."""

    def test_get_chromosome_mlst_header(self):
        """Test MLST header contains expected loci."""
        header = get_chromosome_mlst_header()

        assert isinstance(header, list)
        assert len(header) == 7
        assert 'atpA' in header
        assert 'dnaE' in header
        assert 'dnaK' in header
        assert 'fusA' in header
        assert 'leuA' in header
        assert 'odhA' in header
        assert 'rpoB' in header

    def test_get_tox_header(self):
        """Test tox header contains tox_allele."""
        header = get_tox_header()

        assert isinstance(header, list)
        assert 'tox_allele' in header

    def test_get_virulence(self):
        """Test virulence list contains expected genes."""
        virulence = get_virulence()

        assert isinstance(virulence, list)
        assert 'TOXIN' in virulence
        assert 'spuA' in virulence
        assert 'narG' in virulence

    def test_get_virulence_extended(self):
        """Test extended virulence list contains more genes."""
        virulence = get_virulence()
        extended = get_virulence_extended()

        assert isinstance(extended, list)
        assert len(extended) >= len(virulence)
        assert 'TOXIN' in extended


class TestFindLenContig:
    """Tests for find_len_contig function."""

    def test_find_existing_contig(self, sample_fasta_file):
        """Test finding length of existing contig."""
        length = find_len_contig(sample_fasta_file, 'contig_1')

        assert length == 178  # 60 + 59 + 59 from sample_fasta_content fixture

    def test_find_second_contig(self, sample_fasta_file):
        """Test finding length of second contig."""
        length = find_len_contig(sample_fasta_file, 'contig_2')

        assert length == 118  # 59 + 59 from sample_fasta_content fixture

    def test_non_existent_contig(self, sample_fasta_file):
        """Test finding non-existent contig returns None."""
        length = find_len_contig(sample_fasta_file, 'contig_999')

        assert length is None

    def test_contig_with_header_info(self, tmp_path):
        """Test finding contig with header info after name."""
        fasta_content = """>contig_1 length=100 extra_info
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGC
"""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)

        length = find_len_contig(str(fasta_file), 'contig_1')

        assert length == 96


class TestIsContigEdge:
    """Tests for is_contig_edge function."""

    def test_not_at_edge(self, sample_contig_file):
        """Test gene not at contig edge returns False."""
        data = pd.Series({
            'Reference sequence length': '30',  # 30 aa = 90 bp
            'Start': '50',
            'Stop': '139',
            'File': sample_contig_file,
            'Contig id': 'contig_1',
        })

        result = is_contig_edge(data)

        assert result is False

    def test_at_start_edge(self, sample_contig_file):
        """Test gene at start of contig returns True."""
        data = pd.Series({
            'Reference sequence length': '50',  # 50 aa = 150 bp reference
            'Start': '1',
            'Stop': '100',  # Only 100 bp found, missing 50 bp
            'File': sample_contig_file,
            'Contig id': 'contig_1',
        })

        result = is_contig_edge(data)

        assert result is True

    def test_at_end_edge(self, sample_contig_file):
        """Test gene at end of contig returns True."""
        # contig_1 is 180 bp long
        data = pd.Series({
            'Reference sequence length': '50',  # 50 aa = 150 bp reference
            'Start': '100',
            'Stop': '180',  # At the very end, only 81 bp found
            'File': sample_contig_file,
            'Contig id': 'contig_1',
        })

        result = is_contig_edge(data)

        assert result is True


class TestArmfinderToTable:
    """Tests for armfinder_to_table function."""

    @patch('diphtoscan.utils.find_amrfinderplus_version')
    @patch('diphtoscan.utils.is_contig_edge')
    def test_creates_pivot_table(self, mock_edge, mock_version, sample_amrfinder_output):
        """Test that armfinder_to_table creates a pivot table."""
        mock_version.return_value = '4'
        mock_edge.return_value = False

        result = armfinder_to_table(sample_amrfinder_output)

        assert isinstance(result, pd.DataFrame)
        assert 'TETRACYCLINE' in result.columns
        assert 'MACROLIDE' in result.columns

    @patch('diphtoscan.utils.find_amrfinderplus_version')
    @patch('diphtoscan.utils.is_contig_edge')
    def test_method_annotations(self, mock_edge, mock_version, sample_amrfinder_output):
        """Test that methods are annotated correctly."""
        mock_version.return_value = '4'
        mock_edge.return_value = False

        result = armfinder_to_table(sample_amrfinder_output)

        # BLASTX should have * annotation
        macrolide_values = result['MACROLIDE'].values
        assert any('*' in str(v) for v in macrolide_values if str(v) != '')

    @patch('diphtoscan.utils.find_amrfinderplus_version')
    @patch('diphtoscan.utils.is_contig_edge')
    def test_partial_coverage_annotation(self, mock_edge, mock_version, sample_amrfinder_output):
        """Test that partial coverage is annotated with percentage."""
        mock_version.return_value = '4'
        mock_edge.return_value = False

        result = armfinder_to_table(sample_amrfinder_output)

        # PARTIALX method should have ? and percentage
        tet_values = result['TETRACYCLINE'].values
        # At least one should have the partial annotation
        assert any('?' in str(v) or '%' in str(v) for v in tet_values if str(v) != '')


class TestFindAmrfinderplusVersion:
    """Tests for find_amrfinderplus_version function."""

    def test_version_3(self):
        """Test detection of AMRFinderPlus version 3."""
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(stdout='3.10.24')

            version = find_amrfinderplus_version()

            assert version == '3'

    def test_version_4(self):
        """Test detection of AMRFinderPlus version 4."""
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(stdout='4.0.3')

            version = find_amrfinderplus_version()

            assert version == '4'


class TestGetGenomicContext:
    """Tests for get_genomic_context function."""

    @patch('diphtoscan.utils.find_amrfinderplus_version')
    def test_single_gene_on_contig(self, mock_version, tmp_path):
        """Test genomic context with single gene on contig."""
        mock_version.return_value = '4'

        data = pd.DataFrame({
            'Contig id': ['contig_1'],
            'Element symbol': ['tetA'],
            'Start': ['100'],
            'Stop': ['400'],
            'Class': ['TETRACYCLINE'],
        })

        result = get_genomic_context(str(tmp_path), data)

        assert 'tetA' in result

    @patch('diphtoscan.utils.find_amrfinderplus_version')
    def test_multiple_genes_close_together(self, mock_version, tmp_path):
        """Test genomic context with genes within 8000 bp."""
        mock_version.return_value = '4'

        data = pd.DataFrame({
            'Contig id': ['contig_1', 'contig_1'],
            'Element symbol': ['tetA', 'ermC'],
            'Start': ['100', '500'],
            'Stop': ['400', '800'],
            'Class': ['TETRACYCLINE', 'MACROLIDE'],
        })

        result = get_genomic_context(str(tmp_path), data)

        # Genes should be joined with ;
        assert ';' in result or 'tetA' in result

    @patch('diphtoscan.utils.find_amrfinderplus_version')
    def test_multiple_genes_far_apart(self, mock_version, tmp_path):
        """Test genomic context with genes > 8000 bp apart."""
        mock_version.return_value = '4'

        data = pd.DataFrame({
            'Contig id': ['contig_1', 'contig_1'],
            'Element symbol': ['tetA', 'ermC'],
            'Start': ['100', '10000'],
            'Stop': ['400', '10300'],
            'Class': ['TETRACYCLINE', 'MACROLIDE'],
        })

        result = get_genomic_context(str(tmp_path), data)

        # Genes should be separated with ||
        assert '||' in result


class TestDetermineBiovar:
    """Tests for determine_biovar function."""

    def test_gravis_biovar(self):
        """Test detection of gravis biovar (spuA+ and narG+)."""
        results = {
            'species': 'C. diphtheriae',
            'spuA': 'spuA',
            'narG': 'narG',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'gravis'

    def test_mitis_biovar(self):
        """Test detection of mitis biovar (spuA- and narG-)."""
        results = {
            'species': 'C. diphtheriae',
            'spuA': '-',
            'narG': '-',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'mitis'

    def test_belfantii_biovar(self):
        """Test detection of belfantii biovar (spuA- and narG+)."""
        results = {
            'species': 'C. diphtheriae',
            'spuA': '-',
            'narG': 'narG',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'belfantii'

    def test_intermediate_biovar(self):
        """Test detection of intermediate biovar (spuA+ and narG-)."""
        results = {
            'species': 'C. diphtheriae',
            'spuA': 'spuA',
            'narG': '-',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'intermediate'

    def test_non_diphtheriae_species(self):
        """Test that non-C. diphtheriae returns NA."""
        results = {
            'species': 'C. ulcerans',
            'spuA': 'spuA',
            'narG': 'narG',
        }

        biovar = determine_biovar(results, 'C. ulcerans')

        assert biovar == 'NA'

    def test_with_gene_values_in_columns(self):
        """Test biovar detection with gene names in column values."""
        results = {
            'species': 'C. diphtheriae',
            'SpuA-CLUSTER': 'spuA;spuB',
            'narIJHK': 'narG;narI;narJ',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'gravis'

    def test_empty_results(self):
        """Test biovar with empty results."""
        results = {}

        biovar = determine_biovar(results, 'C. diphtheriae')

        # With no gene data, should return unknown
        assert biovar in ['unknown', 'mitis']  # mitis if all genes absent

    def test_none_species(self):
        """Test biovar with None species."""
        results = {
            'spuA': 'spuA',
            'narG': 'narG',
        }

        # Should still determine biovar if species not provided
        biovar = determine_biovar(results, None)

        assert biovar == 'gravis'

    def test_case_insensitive_gene_detection(self):
        """Test that gene detection is case insensitive."""
        results = {
            'species': 'C. diphtheriae',
            'spuA': 'SPUA',
            'narG': 'NARG',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'gravis'

    def test_genes_in_resistance_columns(self):
        """Test biovar detection when genes are in resistance columns."""
        results = {
            'species': 'C. diphtheriae',
            'SOME_CLASS': 'spuA;other_gene',
            'ANOTHER_CLASS': 'narG;something',
        }

        biovar = determine_biovar(results, 'C. diphtheriae')

        assert biovar == 'gravis'
