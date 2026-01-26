"""Tests for diphtoscan/species.py - species identification functions."""

import pytest
from unittest.mock import patch, MagicMock

from diphtoscan.species import (
    get_corynebacterium_species,
    get_species_results,
    is_cd_complex,
)


class TestGetCorynebacteriumSpecies:
    """Tests for get_corynebacterium_species function."""

    def test_strong_match(self, mock_mash_output):
        """Test strong match (distance <= 0.05)."""
        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mock_mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            assert species == 'C. diphtheriae'
            assert strength == 'strong'

    def test_weak_match(self):
        """Test weak match (0.05 < distance <= 0.1)."""
        mash_output = [
            "Corynebacterium/C.ulcerans_ref.fasta\tsample.fasta\t0.08\t1.0e-30\t800/1000\n",
        ]

        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            assert species == 'C. ulcerans'
            assert strength == 'weak'

    def test_no_match(self):
        """Test no match (distance > 0.1)."""
        mash_output = [
            "Corynebacterium/C.other_ref.fasta\tsample.fasta\t0.25\t1.0e-5\t200/1000\n",
        ]

        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            assert species == 'unknown'
            assert strength == ''

    def test_species_name_expansion(self):
        """Test that abbreviated species names are expanded."""
        mash_output = [
            "Corynebacterium/C.pseudotub_ref.fasta\tsample.fasta\t0.02\t1.0e-50\t1000/1000\n",
        ]

        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            assert 'pseudotuberculosis' in species
            assert strength == 'strong'

    def test_species_name_formatting(self):
        """Test that C. is expanded to C. with space."""
        mash_output = [
            "Corynebacterium/C.diphtheriae_ref.fasta\tsample.fasta\t0.01\t1.0e-50\t1000/1000\n",
        ]

        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            assert species.startswith('C. ')

    def test_empty_mash_output(self):
        """Test handling of empty Mash output."""
        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter([])
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            assert species == 'unknown'
            assert strength == ''

    def test_malformed_mash_output_line(self):
        """Test handling of malformed Mash output lines."""
        mash_output = [
            "incomplete\tline\n",  # Less than 4 fields
            "Corynebacterium/C.diphtheriae_ref.fasta\tsample.fasta\t0.01\t1.0e-50\t1000/1000\n",
        ]

        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            species, strength = get_corynebacterium_species(
                'sample.fasta', '/path/to/folder', '4'
            )

            # Should skip malformed line and use valid one
            assert species == 'C. diphtheriae'


class TestGetSpeciesResults:
    """Tests for get_species_results function."""

    def test_returns_dict_format(self, mock_mash_output):
        """Test that results are returned as a dict with correct keys."""
        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mock_mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            result = get_species_results('sample.fasta', '/path/to/folder', '4')

            assert isinstance(result, dict)
            assert 'species' in result
            assert 'species_match' in result

    def test_contains_species_and_strength(self, mock_mash_output):
        """Test that result contains species name and match strength."""
        with patch('os.popen') as mock_popen:
            mock_file = MagicMock()
            mock_file.__iter__ = lambda self: iter(mock_mash_output)
            mock_file.close = MagicMock()
            mock_popen.return_value = mock_file

            result = get_species_results('sample.fasta', '/path/to/folder', '4')

            assert result['species'] == 'C. diphtheriae'
            assert result['species_match'] == 'strong'


class TestIsCdComplex:
    """Tests for is_cd_complex function."""

    def test_c_diphtheriae(self):
        """Test that C. diphtheriae is in CD-complex."""
        result = {'species': 'C. diphtheriae'}
        assert is_cd_complex(result) is True

    def test_c_belfantii(self):
        """Test that C. belfantii is in CD-complex."""
        result = {'species': 'C. belfantii'}
        assert is_cd_complex(result) is True

    def test_c_rouxii(self):
        """Test that C. rouxii is in CD-complex."""
        result = {'species': 'C. rouxii'}
        assert is_cd_complex(result) is True

    def test_c_ulcerans(self):
        """Test that C. ulcerans is in CD-complex."""
        result = {'species': 'C. ulcerans'}
        assert is_cd_complex(result) is True

    def test_c_pseudotuberculosis(self):
        """Test that C. pseudotuberculosis is in CD-complex."""
        result = {'species': 'C. pseudotuberculosis'}
        assert is_cd_complex(result) is True

    def test_c_ramonii(self):
        """Test that C. ramonii is in CD-complex."""
        result = {'species': 'C. ramonii'}
        assert is_cd_complex(result) is True

    def test_non_cd_complex_species(self):
        """Test that non-CD-complex species return False."""
        result = {'species': 'C. glutamicum'}
        assert is_cd_complex(result) is False

    def test_unknown_species(self):
        """Test that unknown species return False."""
        result = {'species': 'unknown'}
        assert is_cd_complex(result) is False

    def test_partial_match(self):
        """Test that species name with prefix still matches."""
        result = {'species': 'C. diphtheriae bv. gravis'}
        assert is_cd_complex(result) is True

    def test_missing_species_key_raises(self):
        """Test that missing species key raises AssertionError."""
        result = {}
        with pytest.raises(AssertionError):
            is_cd_complex(result)
