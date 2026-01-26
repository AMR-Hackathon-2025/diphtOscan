"""Tests for diphtoscan/misc.py - complement and reverse complement functions."""

import pytest
from diphtoscan.misc import complement_base, reverse_complement, REV_COMP_DICT


class TestComplementBase:
    """Tests for complement_base function."""

    def test_standard_bases_uppercase(self):
        """Test complement of standard uppercase bases."""
        assert complement_base('A') == 'T'
        assert complement_base('T') == 'A'
        assert complement_base('G') == 'C'
        assert complement_base('C') == 'G'

    def test_standard_bases_lowercase(self):
        """Test complement of standard lowercase bases."""
        assert complement_base('a') == 't'
        assert complement_base('t') == 'a'
        assert complement_base('g') == 'c'
        assert complement_base('c') == 'g'

    def test_ambiguous_bases_uppercase(self):
        """Test complement of IUPAC ambiguous bases (uppercase)."""
        assert complement_base('R') == 'Y'  # A or G -> T or C
        assert complement_base('Y') == 'R'  # T or C -> A or G
        assert complement_base('S') == 'S'  # G or C -> C or G (self-complementary)
        assert complement_base('W') == 'W'  # A or T -> T or A (self-complementary)
        assert complement_base('K') == 'M'  # G or T -> C or A
        assert complement_base('M') == 'K'  # A or C -> T or G
        assert complement_base('B') == 'V'  # C, G, or T -> G, C, or A
        assert complement_base('V') == 'B'  # A, C, or G -> T, G, or C
        assert complement_base('D') == 'H'  # A, G, or T -> T, C, or A
        assert complement_base('H') == 'D'  # A, C, or T -> T, G, or A
        assert complement_base('N') == 'N'  # Any base (self-complementary)

    def test_ambiguous_bases_lowercase(self):
        """Test complement of IUPAC ambiguous bases (lowercase)."""
        assert complement_base('r') == 'y'
        assert complement_base('y') == 'r'
        assert complement_base('s') == 's'
        assert complement_base('w') == 'w'
        assert complement_base('k') == 'm'
        assert complement_base('m') == 'k'
        assert complement_base('n') == 'n'

    def test_gap_characters(self):
        """Test complement of gap characters."""
        assert complement_base('.') == '.'
        assert complement_base('-') == '-'
        assert complement_base('?') == '?'

    def test_unknown_character(self):
        """Test complement of unknown characters returns 'N'."""
        assert complement_base('X') == 'N'
        assert complement_base('Z') == 'N'
        assert complement_base('1') == 'N'
        assert complement_base('!') == 'N'
        assert complement_base(' ') == 'N'


class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_simple_sequence(self):
        """Test reverse complement of a simple sequence."""
        assert reverse_complement('ATGC') == 'GCAT'

    def test_single_base(self):
        """Test reverse complement of a single base."""
        assert reverse_complement('A') == 'T'
        assert reverse_complement('T') == 'A'
        assert reverse_complement('G') == 'C'
        assert reverse_complement('C') == 'G'

    def test_empty_sequence(self):
        """Test reverse complement of empty sequence."""
        assert reverse_complement('') == ''

    def test_palindrome(self):
        """Test reverse complement of a palindromic sequence."""
        # GAATTC is an EcoRI site - a biological palindrome
        assert reverse_complement('GAATTC') == 'GAATTC'

    def test_mixed_case(self):
        """Test reverse complement with mixed case."""
        assert reverse_complement('AtGc') == 'gCaT'

    def test_longer_sequence(self):
        """Test reverse complement of a longer sequence."""
        seq = 'ATGAAACGCATTAGCACC'
        expected = 'GGTGCTAATGCGTTTCAT'
        assert reverse_complement(seq) == expected

    def test_with_ambiguous_bases(self):
        """Test reverse complement with ambiguous bases."""
        assert reverse_complement('ATGN') == 'NCAT'
        assert reverse_complement('RYSWKM') == 'KMWSRY'

    def test_with_gaps(self):
        """Test reverse complement with gap characters."""
        assert reverse_complement('ATG-C') == 'G-CAT'
        assert reverse_complement('A.TG') == 'CA.T'

    def test_preserves_length(self):
        """Test that reverse complement preserves sequence length."""
        seq = 'ATGCATGCATGC'
        result = reverse_complement(seq)
        assert len(result) == len(seq)

    def test_double_reverse_complement_returns_original(self):
        """Test that applying reverse complement twice returns the original."""
        seq = 'ATGCNNRYATGC'
        assert reverse_complement(reverse_complement(seq)) == seq


class TestRevCompDict:
    """Tests for the REV_COMP_DICT constant."""

    def test_dict_contains_standard_bases(self):
        """Test that REV_COMP_DICT contains all standard bases."""
        for base in 'ATGCatgc':
            assert base in REV_COMP_DICT

    def test_dict_contains_ambiguous_bases(self):
        """Test that REV_COMP_DICT contains IUPAC ambiguous bases."""
        for base in 'RYSWKMBVDHNryswkmbvdhn':
            assert base in REV_COMP_DICT

    def test_dict_contains_gap_characters(self):
        """Test that REV_COMP_DICT contains gap characters."""
        for char in '.-?':
            assert char in REV_COMP_DICT
