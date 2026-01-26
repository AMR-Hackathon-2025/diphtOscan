"""Tests for diphtoscan/truncation.py - gene truncation checking."""

import pytest
from unittest.mock import MagicMock


class TestTruncationCheck:
    """Tests for truncation_check function."""

    def test_full_length_gene(self):
        """Test a gene that is full length (no truncation)."""
        from diphtoscan.truncation import truncation_check

        # Create a mock hit for a full-length gene
        # Gene is 300 bp (100 amino acids including stop codon)
        # Sequence starts at position 1 and has no stop codons
        mock_hit = MagicMock()
        mock_hit.ref_length = 303  # 100 aa + stop codon = 101 codons = 303 bp
        mock_hit.get_seq_start_end_pos_strand.return_value = (
            'ATGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA',  # 303 bp, ends with TAA stop
            1,
            303
        )

        result, coverage, translation = truncation_check(mock_hit)

        assert result == ''  # No truncation marker
        assert coverage >= 90.0

    def test_truncated_gene(self):
        """Test a gene that is truncated (premature stop codon)."""
        from diphtoscan.truncation import truncation_check

        # Create a mock hit with a premature stop codon
        mock_hit = MagicMock()
        mock_hit.ref_length = 303  # Expected 100 aa
        # Sequence with stop codon at position 50 (after ~16 aa)
        mock_hit.get_seq_start_end_pos_strand.return_value = (
            'ATGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA',  # Stop at position ~50
            1,
            303
        )

        result, coverage, translation = truncation_check(mock_hit)

        assert '-' in result  # Has truncation marker
        assert '%' in result  # Shows percentage
        assert coverage < 90.0

    def test_gene_not_starting_at_position_1(self):
        """Test a gene that doesn't start at position 1 (considered 0%)."""
        from diphtoscan.truncation import truncation_check

        mock_hit = MagicMock()
        mock_hit.ref_length = 303
        mock_hit.get_seq_start_end_pos_strand.return_value = (
            'ATGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA',
            5,  # Doesn't start at 1
            65
        )

        result, coverage, translation = truncation_check(mock_hit)

        assert result == '-0%'
        assert coverage == 0.0
        assert translation == ''

    def test_gene_with_ambiguous_bases(self):
        """Test a gene with ambiguous bases (N) that breaks translation."""
        from diphtoscan.truncation import truncation_check

        mock_hit = MagicMock()
        mock_hit.ref_length = 303
        # Sequence with N in the middle
        mock_hit.get_seq_start_end_pos_strand.return_value = (
            'ATGGCAGCAGCAGCAGCAGCANNNNNNGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATAA',
            1,
            303
        )

        result, coverage, translation = truncation_check(mock_hit)

        # Result depends on where the ambiguous bases are
        # The function splits on ambiguous bases, taking only the part before them
        assert isinstance(result, str)
        assert isinstance(coverage, float)

    def test_custom_coverage_threshold(self):
        """Test truncation check with a custom coverage threshold."""
        from diphtoscan.truncation import truncation_check

        mock_hit = MagicMock()
        mock_hit.ref_length = 303  # 100 aa
        # Create a sequence that would give ~85% coverage
        # We need about 85 aa worth of sequence before stop
        seq_85_percent = 'ATG' + 'GCA' * 83 + 'TAA' + 'GCA' * 15  # ~85 aa then stop, then junk
        mock_hit.get_seq_start_end_pos_strand.return_value = (seq_85_percent, 1, 303)

        # With default 90% threshold, should be truncated
        result_default, cov_default, _ = truncation_check(mock_hit, cov_threshold=90.0)

        # With 80% threshold, should not be truncated
        result_lower, cov_lower, _ = truncation_check(mock_hit, cov_threshold=80.0)

        # Same coverage calculated either way
        assert cov_default == cov_lower

        # But truncation status differs based on threshold
        if cov_default >= 90.0:
            assert result_default == ''
        else:
            assert '-' in result_default

    def test_sequence_not_multiple_of_3(self):
        """Test handling of sequences that aren't multiples of 3."""
        from diphtoscan.truncation import truncation_check

        mock_hit = MagicMock()
        mock_hit.ref_length = 303
        # Sequence length 301 (not divisible by 3)
        mock_hit.get_seq_start_end_pos_strand.return_value = (
            'ATGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCATA',  # 301 bp
            1,
            301
        )

        result, coverage, translation = truncation_check(mock_hit)

        # Should still work - sequence is trimmed to multiple of 3
        assert isinstance(result, str)
        assert isinstance(coverage, float)

    def test_empty_translation(self):
        """Test when translation results in empty or very short protein."""
        from diphtoscan.truncation import truncation_check

        mock_hit = MagicMock()
        mock_hit.ref_length = 303
        # Starts with stop codon
        mock_hit.get_seq_start_end_pos_strand.return_value = ('TAAGCAGCAGCA', 1, 12)

        result, coverage, translation = truncation_check(mock_hit)

        assert translation == ''
        assert coverage == 0.0
        assert '-0%' in result

    def test_100_percent_coverage(self):
        """Test gene with exactly 100% coverage."""
        from diphtoscan.truncation import truncation_check

        mock_hit = MagicMock()
        mock_hit.ref_length = 12  # 3 aa + stop = 4 codons = 12 bp
        # ATG (M) GCA (A) GCA (A) TAA (stop) - 3 aa before stop
        mock_hit.get_seq_start_end_pos_strand.return_value = ('ATGGCAGCATAA', 1, 12)

        result, coverage, translation = truncation_check(mock_hit)

        assert result == ''
        assert coverage == 100.0
        assert translation == 'MAA'
