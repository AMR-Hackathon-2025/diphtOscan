"""Tests for diphtoscan/download_alleles_st.py - BIGSdb API interactions."""

import pytest
from unittest.mock import patch, MagicMock
import os

from diphtoscan.download_alleles_st import (
    download_alleles,
    create_db,
    download_profiles_st,
    download_profiles_tox,
    BASE_URI,
)


class TestDownloadAlleles:
    """Tests for download_alleles function."""

    def test_creates_directory_if_not_exists(self, tmp_path):
        """Test that the function creates the output directory."""
        output_dir = tmp_path / "new_dir"

        with patch('requests.get') as mock_get:
            mock_get.return_value = MagicMock(status_code=200, json=lambda: {'loci': []})
            download_alleles('pubmlst_diphtheria_seqdef', '3', str(output_dir))
            assert output_dir.exists()

    def test_handles_404_database(self, tmp_path, capsys):
        """Test that 404 for database prints message and exits."""
        with patch('requests.get') as mock_get:
            mock_get.return_value = MagicMock(status_code=404)

            with patch('diphtoscan.download_alleles_st.sys.exit', side_effect=SystemExit(1)) as mock_exit:
                with pytest.raises(SystemExit):
                    download_alleles('nonexistent_db', '3', str(tmp_path))
                mock_exit.assert_called_with(1)
                captured = capsys.readouterr()
                assert 'does not exist' in captured.out

    def test_handles_404_scheme(self, tmp_path, capsys):
        """Test that 404 for scheme prints message and exits."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=404),
            ]

            with patch('diphtoscan.download_alleles_st.sys.exit', side_effect=SystemExit(1)) as mock_exit:
                with pytest.raises(SystemExit):
                    download_alleles('pubmlst_diphtheria_seqdef', '999', str(tmp_path))
                mock_exit.assert_called_with(1)
                captured = capsys.readouterr()
                assert 'does not exist' in captured.out

    def test_downloads_loci_from_scheme(self, tmp_path):
        """Test downloading loci from a scheme."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, json=lambda: {
                    'loci': ['https://bigsdb.pasteur.fr/api/db/test/loci/atpA']
                }),
                MagicMock(status_code=200, json=lambda: {
                    'id': 'atpA',
                    'alleles_fasta': 'https://bigsdb.pasteur.fr/api/db/test/loci/atpA/alleles_fasta'
                }),
                MagicMock(status_code=200, text='>atpA_1\nATGC\n'),
            ]

            result = download_alleles('test_db', '3', str(tmp_path))

            assert 'atpA' in result
            assert (tmp_path / 'atpA.fas').exists()

    def test_returns_list_of_loci_names(self, tmp_path):
        """Test that function returns list of downloaded loci names."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, json=lambda: {
                    'loci': ['https://bigsdb.pasteur.fr/api/db/test/loci/atpA']
                }),
                MagicMock(status_code=200, json=lambda: {
                    'id': 'atpA',
                    'alleles_fasta': 'https://bigsdb.pasteur.fr/api/db/test/loci/atpA/alleles_fasta'
                }),
                MagicMock(status_code=200, text='>atpA_1\nATGC\n'),
            ]

            result = download_alleles('test_db', '3', str(tmp_path))

            assert isinstance(result, list)
            assert 'atpA' in result


class TestCreateDb:
    """Tests for create_db function."""

    def test_creates_sequences_directory(self, tmp_path):
        """Test that create_db creates a sequences subdirectory."""
        with patch('diphtoscan.download_alleles_st.download_alleles') as mock_download:
            mock_download.return_value = ['atpA']

            create_db('test_db', '3', str(tmp_path))

            mock_download.assert_called_once()
            call_args = mock_download.call_args[0]
            assert 'sequences' in call_args[2]

    def test_returns_database_path_and_loci(self, tmp_path):
        """Test that create_db returns database path and loci list."""
        with patch('diphtoscan.download_alleles_st.download_alleles') as mock_download:
            mock_download.return_value = ['atpA', 'dnaE']

            db_path, loci = create_db('test_db', '3', str(tmp_path))

            assert 'test_db' in db_path
            assert 'scheme_3' in db_path
            assert loci == ['atpA', 'dnaE']

    def test_creates_combined_database_file(self, tmp_path):
        """Test that a combined database file is created."""
        # Create mock sequence files
        seq_dir = tmp_path / "sequences"
        seq_dir.mkdir()
        (seq_dir / "atpA.fas").write_text(">atpA_1\nATGC\n")
        (seq_dir / "dnaE.fas").write_text(">dnaE_1\nGCTA\n")

        with patch('diphtoscan.download_alleles_st.download_alleles') as mock_download:
            mock_download.return_value = ['atpA', 'dnaE']

            db_path, loci = create_db('test_db', '3', str(tmp_path))

            assert os.path.exists(db_path)


class TestDownloadProfilesSt:
    """Tests for download_profiles_st function."""

    def test_creates_directory_if_not_exists(self, tmp_path):
        """Test that directory is created if it doesn't exist."""
        output_dir = tmp_path / "new_profiles_dir"

        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, text="ST\tatpA\tdnaE\n1\t1\t1\n"),
            ]

            download_profiles_st('test_db', '3', str(output_dir), ['atpA', 'dnaE'])

            assert output_dir.exists()

    def test_handles_404_database(self, tmp_path, capsys):
        """Test that 404 for database prints message and exits."""
        with patch('requests.get') as mock_get:
            mock_get.return_value = MagicMock(status_code=404)

            with patch('diphtoscan.download_alleles_st.sys.exit', side_effect=SystemExit(1)) as mock_exit:
                with pytest.raises(SystemExit):
                    download_profiles_st('nonexistent_db', '3', str(tmp_path), ['atpA'])
                mock_exit.assert_called_with(1)
                captured = capsys.readouterr()
                assert 'does not exist' in captured.out

    def test_creates_profiles_file(self, tmp_path):
        """Test that ST profiles file is created."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, text="ST\tatpA\tdnaE\n1\t1\t1\n2\t1\t2\n"),
            ]

            result = download_profiles_st('test_db', '3', str(tmp_path), ['atpA', 'dnaE'])

            assert 'st_profiles.txt' in result
            assert os.path.exists(result)

    def test_returns_profiles_file_path(self, tmp_path):
        """Test that function returns the profiles file path."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, text="ST\tatpA\n1\t1\n"),
            ]

            result = download_profiles_st('test_db', '3', str(tmp_path), ['atpA'])

            assert result == str(tmp_path) + '/st_profiles.txt'


class TestDownloadProfilesTox:
    """Tests for download_profiles_tox function."""

    def test_creates_directory_if_not_exists(self, tmp_path):
        """Test that directory is created if it doesn't exist."""
        output_dir = tmp_path / "new_tox_dir"

        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, text="ST\ttox\n1\t1\n"),
            ]

            download_profiles_tox('test_db', '4', str(output_dir))

            assert output_dir.exists()

    def test_creates_tox_profiles_file(self, tmp_path):
        """Test that tox profiles file is created."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, text="ST\ttox\n1\t1\n2\t2\n"),
            ]

            result = download_profiles_tox('test_db', '4', str(tmp_path))

            assert 'tox_profiles.txt' in result

    def test_returns_tox_profiles_file_path(self, tmp_path):
        """Test that function returns the tox profiles file path."""
        with patch('requests.get') as mock_get:
            mock_get.side_effect = [
                MagicMock(status_code=200),
                MagicMock(status_code=200, text="ST\ttox\n1\t1\n"),
            ]

            result = download_profiles_tox('test_db', '4', str(tmp_path))

            assert result == str(tmp_path) + '/tox_profiles.txt'


class TestBaseUri:
    """Tests for BASE_URI constant."""

    def test_base_uri_is_bigsdb(self):
        """Test that BASE_URI points to BIGSdb API."""
        assert 'bigsdb.pasteur.fr' in BASE_URI
        assert '/api' in BASE_URI
