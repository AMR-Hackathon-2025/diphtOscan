"""Tests for diphtoscan/assembly_qc.py - assembly quality control metrics."""

import pytest
from pathlib import Path

from diphtoscan.assembly_qc import (
    parse_fasta_sequences,
    calculate_n50_l50,
    calculate_gc_content,
    calculate_n_content,
    calculate_assembly_stats,
    get_qc_header,
    format_qc_results,
)


class TestParseFastaSequences:
    """Tests for parse_fasta_sequences function."""

    def test_single_contig(self, tmp_path):
        """Test parsing FASTA with single contig."""
        fasta_content = """>contig_1
ATGCATGCATGC
ATGCATGCATGC
"""
        fasta_file = tmp_path / "single.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta_sequences(str(fasta_file))

        assert len(sequences) == 1
        assert sequences[0][0] == 'contig_1'
        assert sequences[0][1] == 'ATGCATGCATGCATGCATGCATGC'

    def test_multiple_contigs(self, tmp_path):
        """Test parsing FASTA with multiple contigs."""
        fasta_content = """>contig_1
ATGCATGC
>contig_2
GCTAGCTA
>contig_3
AAAAAAA
"""
        fasta_file = tmp_path / "multi.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta_sequences(str(fasta_file))

        assert len(sequences) == 3
        assert sequences[0][0] == 'contig_1'
        assert sequences[1][0] == 'contig_2'
        assert sequences[2][0] == 'contig_3'

    def test_empty_file(self, tmp_path):
        """Test parsing empty FASTA file."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        sequences = parse_fasta_sequences(str(fasta_file))

        assert len(sequences) == 0

    def test_converts_to_uppercase(self, tmp_path):
        """Test that sequences are converted to uppercase."""
        fasta_content = """>contig_1
atgcATGC
"""
        fasta_file = tmp_path / "lowercase.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta_sequences(str(fasta_file))

        assert sequences[0][1] == 'ATGCATGC'

    def test_header_with_description(self, tmp_path):
        """Test parsing header with description."""
        fasta_content = """>contig_1 length=100 description
ATGCATGC
"""
        fasta_file = tmp_path / "desc.fasta"
        fasta_file.write_text(fasta_content)

        sequences = parse_fasta_sequences(str(fasta_file))

        # Should only keep the first word of header
        assert sequences[0][0] == 'contig_1'


class TestCalculateN50L50:
    """Tests for calculate_n50_l50 function."""

    def test_simple_case(self):
        """Test N50/L50 with simple input."""
        # Total: 100, half: 50
        # Sorted desc: [40, 30, 20, 10]
        # Cumulative: 40, 70 (>=50)
        # N50 = 30, L50 = 2
        lengths = [10, 20, 30, 40]
        n50, l50 = calculate_n50_l50(lengths)

        assert n50 == 30
        assert l50 == 2

    def test_single_contig(self):
        """Test N50/L50 with single contig."""
        lengths = [1000]
        n50, l50 = calculate_n50_l50(lengths)

        assert n50 == 1000
        assert l50 == 1

    def test_empty_list(self):
        """Test N50/L50 with empty list."""
        n50, l50 = calculate_n50_l50([])

        assert n50 == 0
        assert l50 == 0

    def test_all_equal_lengths(self):
        """Test N50/L50 with equal length contigs."""
        # 5 contigs of 100 each = 500 total, half = 250
        # Need 3 contigs to reach 300 >= 250
        lengths = [100, 100, 100, 100, 100]
        n50, l50 = calculate_n50_l50(lengths)

        assert n50 == 100
        assert l50 == 3

    def test_two_contigs(self):
        """Test N50/L50 with two contigs."""
        # Total: 150, half: 75
        # Sorted: [100, 50], cumulative: 100 >= 75
        # N50 = 100, L50 = 1
        lengths = [50, 100]
        n50, l50 = calculate_n50_l50(lengths)

        assert n50 == 100
        assert l50 == 1


class TestCalculateGcContent:
    """Tests for calculate_gc_content function."""

    def test_all_gc(self):
        """Test 100% GC content."""
        result = calculate_gc_content("GGGGCCCC")
        assert result == 100.0

    def test_no_gc(self):
        """Test 0% GC content."""
        result = calculate_gc_content("AAAAATTTT")
        assert result == 0.0

    def test_half_gc(self):
        """Test 50% GC content."""
        result = calculate_gc_content("ATGC")
        assert result == 50.0

    def test_empty_sequence(self):
        """Test empty sequence."""
        result = calculate_gc_content("")
        assert result == 0.0

    def test_realistic_gc(self):
        """Test realistic GC content calculation."""
        # C. diphtheriae typically has ~53% GC
        sequence = "GCGCGCGCATATATGCGCGC"  # 14 GC out of 20 = 70%
        result = calculate_gc_content(sequence)
        assert result == 70.0


class TestCalculateNContent:
    """Tests for calculate_n_content function."""

    def test_no_ns(self):
        """Test sequence with no Ns."""
        result = calculate_n_content("ATGCATGC")
        assert result == 0.0

    def test_all_ns(self):
        """Test sequence with all Ns."""
        result = calculate_n_content("NNNNNNNN")
        assert result == 100.0

    def test_some_ns(self):
        """Test sequence with some Ns."""
        result = calculate_n_content("ATGNNNNN")  # 5 Ns out of 8 = 62.5%
        assert result == 62.5

    def test_empty_sequence(self):
        """Test empty sequence."""
        result = calculate_n_content("")
        assert result == 0.0


class TestCalculateAssemblyStats:
    """Tests for calculate_assembly_stats function."""

    def test_basic_assembly(self, tmp_path):
        """Test calculating stats for basic assembly."""
        fasta_content = """>contig_1
ATGCATGCATGCATGCATGC
>contig_2
GCTAGCTAGCTA
>contig_3
AAAAAAAAAA
"""
        fasta_file = tmp_path / "basic.fasta"
        fasta_file.write_text(fasta_content)

        stats = calculate_assembly_stats(str(fasta_file))

        assert stats['total_length'] == 42  # 20 + 12 + 10
        assert stats['num_contigs'] == 3
        assert stats['largest_contig'] == 20
        assert 'n50' in stats
        assert 'l50' in stats
        assert 'gc_percent' in stats
        assert 'n_percent' in stats

    def test_single_contig_assembly(self, tmp_path):
        """Test stats for single contig assembly."""
        fasta_content = """>chromosome
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
"""
        fasta_file = tmp_path / "single.fasta"
        fasta_file.write_text(fasta_content)

        stats = calculate_assembly_stats(str(fasta_file))

        assert stats['num_contigs'] == 1
        assert stats['total_length'] == 60
        assert stats['largest_contig'] == 60
        assert stats['n50'] == 60
        assert stats['l50'] == 1

    def test_assembly_with_ns(self, tmp_path):
        """Test stats for assembly containing Ns."""
        fasta_content = """>contig_1
ATGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAT
"""
        fasta_file = tmp_path / "with_ns.fasta"
        fasta_file.write_text(fasta_content)

        stats = calculate_assembly_stats(str(fasta_file))

        assert stats['n_percent'] > 0

    def test_empty_file_raises_error(self, tmp_path):
        """Test that empty file raises ValueError."""
        fasta_file = tmp_path / "empty.fasta"
        fasta_file.write_text("")

        with pytest.raises(ValueError):
            calculate_assembly_stats(str(fasta_file))

    def test_realistic_assembly(self, tmp_path):
        """Test stats for realistic assembly size."""
        # Create a small assembly with multiple contigs of varying sizes
        contigs = []
        for i, size in enumerate([50000, 30000, 20000, 15000, 10000, 5000]):
            seq = "ATGC" * (size // 4)
            contigs.append(f">contig_{i+1}\n{seq}")

        fasta_file = tmp_path / "realistic.fasta"
        fasta_file.write_text("\n".join(contigs))

        stats = calculate_assembly_stats(str(fasta_file))

        assert stats['total_length'] == 130000
        assert stats['num_contigs'] == 6
        assert stats['largest_contig'] == 50000
        # GC should be 50% since we used equal ATGC
        assert abs(stats['gc_percent'] - 50.0) < 0.1


class TestGetQcHeader:
    """Tests for get_qc_header function."""

    def test_returns_list(self):
        """Test that get_qc_header returns a list."""
        header = get_qc_header()
        assert isinstance(header, list)

    def test_contains_expected_fields(self):
        """Test that header contains expected metric names."""
        header = get_qc_header()

        expected_fields = [
            'total_length',
            'num_contigs',
            'n50',
            'l50',
            'largest_contig',
            'gc_percent',
            'n_percent'
        ]

        for field in expected_fields:
            assert field in header


class TestFormatQcResults:
    """Tests for format_qc_results function."""

    def test_adds_prefix(self):
        """Test that results are prefixed with 'qc_'."""
        stats = {
            'total_length': 1000000,
            'num_contigs': 50,
            'n50': 50000,
        }

        formatted = format_qc_results(stats)

        assert 'qc_total_length' in formatted
        assert 'qc_num_contigs' in formatted
        assert 'qc_n50' in formatted

    def test_converts_to_strings(self):
        """Test that values are converted to strings."""
        stats = {
            'total_length': 1000000,
            'gc_percent': 53.2,
        }

        formatted = format_qc_results(stats)

        assert formatted['qc_total_length'] == '1000000'
        assert formatted['qc_gc_percent'] == '53.2'


class TestIntegration:
    """Integration tests for assembly QC module."""

    def test_full_workflow(self, sample_fasta_file):
        """Test complete QC workflow with sample file."""
        stats = calculate_assembly_stats(sample_fasta_file)
        formatted = format_qc_results(stats)

        # Check all expected fields are present
        expected_fields = [
            'qc_total_length',
            'qc_num_contigs',
            'qc_n50',
            'qc_l50',
            'qc_largest_contig',
            'qc_gc_percent',
            'qc_n_percent'
        ]

        for field in expected_fields:
            assert field in formatted
            assert isinstance(formatted[field], str)
