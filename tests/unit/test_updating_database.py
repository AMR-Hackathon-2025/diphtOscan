"""Tests for diphtoscan/updating_database.py - database update functions."""

import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

from diphtoscan.updating_database import (
    node_class,
    complete_missing_classification,
    download_amrfinder_database,
    update_amrfinderplus_db_file,
    remove_mlst_database,
    update_database,
)


class TestNodeClass:
    """Tests for node_class dictionary."""

    def test_contains_toxin_mapping(self):
        """Test that node_class contains tox -> TOXIN mapping."""
        assert 'tox' in node_class
        assert node_class['tox'] == 'TOXIN'

    def test_contains_pili_mappings(self):
        """Test that node_class contains pili gene mappings."""
        assert 'spaA' in node_class
        assert node_class['spaA'] == 'SpaA-type_pili_diphtheriae'
        assert 'spaD' in node_class
        assert node_class['spaD'] == 'SpaD-type_pili_diphtheriae'
        assert 'spaH' in node_class
        assert node_class['spaH'] == 'SpaH-type_pili_diphtheriae'

    def test_contains_virulence_adhesin_mappings(self):
        """Test that node_class contains virulence/adhesin mappings."""
        assert 'cbpA' in node_class
        assert node_class['cbpA'] == 'VIRULENCE/ADHESIN'
        assert 'nanH' in node_class
        assert node_class['nanH'] == 'VIRULENCE/ADHESIN'

    def test_contains_other_toxins_mapping(self):
        """Test that node_class contains pld -> OTHER_TOXINS mapping."""
        assert 'pld' in node_class
        assert node_class['pld'] == 'OTHER_TOXINS'


class TestCompleteMissingClassification:
    """Tests for complete_missing_classification function."""

    def test_fills_missing_class_values(self, tmp_path):
        """Test that missing class values are filled."""
        # Create a test TSV file
        tsv_content = """#node_id\tparent_node_id\tclass\tsubclass
tox\tVIRULENCE_Cdiphth\t\t
other_gene\tOTHER_PARENT\t\t
"""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text(tsv_content)

        complete_missing_classification(str(tsv_file))

        # Read the updated file
        import pandas as pd
        df = pd.read_csv(str(tsv_file), sep="\t", escapechar="\\")

        # Check that tox got classified
        tox_row = df[df['#node_id'] == 'tox']
        assert tox_row['class'].values[0] == 'TOXIN'

    def test_only_fills_virulence_cdiphth_parent(self, tmp_path):
        """Test that only VIRULENCE_Cdiphth parent nodes are filled."""
        tsv_content = """#node_id\tparent_node_id\tclass\tsubclass
tox\tVIRULENCE_Cdiphth\t\t
other_gene\tOTHER_PARENT\t\t
"""
        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text(tsv_content)

        complete_missing_classification(str(tsv_file))

        import pandas as pd
        df = pd.read_csv(str(tsv_file), sep="\t", escapechar="\\")

        # other_gene should still have empty class
        other_row = df[df['#node_id'] == 'other_gene']
        assert pd.isna(other_row['class'].values[0]) or other_row['class'].values[0] == ''


class TestDownloadAmrfinderDatabase:
    """Tests for download_amrfinder_database function."""

    def test_creates_output_directory(self, tmp_path):
        """Test that output directory is created."""
        output_dir = tmp_path / "new_amr_db"

        with patch('requests.get') as mock_get:
            mock_get.return_value = MagicMock(
                status_code=200,
                iter_lines=lambda: iter([b'<pre>Name', b'<hr></pre>'])
            )

            download_amrfinder_database('https://example.com/', str(output_dir))

            assert output_dir.exists()

    def test_raises_on_failed_download(self, tmp_path):
        """Test that RuntimeError is raised on failed download."""
        with patch('requests.get') as mock_get:
            mock_get.return_value = MagicMock(status_code=500)

            with pytest.raises(RuntimeError) as exc_info:
                download_amrfinder_database('https://example.com/', str(tmp_path))

            assert 'Failed to download' in str(exc_info.value)

    def test_downloads_files_from_listing(self, tmp_path):
        """Test that files are downloaded from the directory listing."""
        html_content = b"""<pre>Name
<a href="file1.txt">file1.txt</a>
<a href="file2.txt">file2.txt</a>
<a href="subdir/">subdir/</a>
<hr></pre>"""

        with patch('requests.get') as mock_get:
            # First call returns directory listing
            # Subsequent calls return file contents
            mock_get.side_effect = [
                MagicMock(status_code=200, iter_lines=lambda: iter(html_content.split(b'\n'))),
                MagicMock(status_code=200, iter_content=lambda chunk_size: [b'file1 content']),
                MagicMock(status_code=200, iter_content=lambda chunk_size: [b'file2 content']),
            ]

            download_amrfinder_database('https://example.com/', str(tmp_path))

            assert (tmp_path / 'file1.txt').exists()
            assert (tmp_path / 'file2.txt').exists()
            # Directories should be skipped
            assert not (tmp_path / 'subdir').exists()


class TestUpdateAmrfinderplusDbFile:
    """Tests for update_amrfinderplus_db_file function."""

    def test_appends_content_to_output(self, tmp_path):
        """Test that content is appended to output file."""
        input_file = tmp_path / "input.txt"
        input_file.write_text("line1\nline2\nline3\n")

        output_file = tmp_path / "output.txt"
        output_file.write_text("existing content\n")

        update_amrfinderplus_db_file(str(input_file), str(output_file))

        content = output_file.read_text()
        assert "existing content" in content
        assert "line1" in content
        assert "line2" in content
        assert "line3" in content

    def test_skip_first_line_option(self, tmp_path):
        """Test that first line can be skipped."""
        input_file = tmp_path / "input.txt"
        input_file.write_text("header\nline1\nline2\n")

        output_file = tmp_path / "output.txt"
        output_file.write_text("")

        update_amrfinderplus_db_file(str(input_file), str(output_file), skip_first_line=True)

        content = output_file.read_text()
        assert "header" not in content
        assert "line1" in content
        assert "line2" in content

    def test_creates_new_file_if_not_exists(self, tmp_path):
        """Test that output file is created if it doesn't exist."""
        input_file = tmp_path / "input.txt"
        input_file.write_text("content\n")

        output_file = tmp_path / "new_output.txt"

        update_amrfinderplus_db_file(str(input_file), str(output_file))

        assert output_file.exists()
        assert "content" in output_file.read_text()


class TestRemoveMlstDatabase:
    """Tests for remove_mlst_database function."""

    def test_removes_files_in_parent_directory(self, tmp_path):
        """Test that files in parent directory are removed."""
        # Create directory structure
        mlst_dir = tmp_path / "mlst"
        mlst_dir.mkdir()
        sequences_dir = mlst_dir / "sequences"
        sequences_dir.mkdir()

        # Create some files
        (mlst_dir / "database.fas").write_text("content")
        (mlst_dir / "profiles.txt").write_text("content")
        (sequences_dir / "atpA.fas").write_text("content")

        db_path = str(mlst_dir / "database.fas")

        remove_mlst_database(db_path)

        # Files should be removed
        assert not (mlst_dir / "database.fas").exists()
        assert not (mlst_dir / "profiles.txt").exists()
        assert not sequences_dir.exists()

    def test_removes_sequences_directory(self, tmp_path):
        """Test that sequences subdirectory is removed."""
        mlst_dir = tmp_path / "mlst"
        mlst_dir.mkdir()
        sequences_dir = mlst_dir / "sequences"
        sequences_dir.mkdir()

        (sequences_dir / "atpA.fas").write_text("content")

        db_path = str(mlst_dir / "database.fas")

        remove_mlst_database(db_path)

        assert not sequences_dir.exists()


class TestUpdateDatabase:
    """Tests for update_database function."""

    def test_does_nothing_when_update_false(self, tmp_path):
        """Test that nothing happens when update flag is False."""
        args = MagicMock()
        args.update = False
        args.path = str(tmp_path)

        mlst_db = ('mlst_header', str(tmp_path / 'mlst.fas'), str(tmp_path / 'st_profiles.txt'))
        tox_db = ('tox_header', str(tmp_path / 'tox.fas'), str(tmp_path / 'tox_profiles.txt'))

        with patch('diphtoscan.updating_database.remove_mlst_database') as mock_remove:
            update_database(args, mlst_db, tox_db)

            mock_remove.assert_not_called()

    def test_calls_update_functions_when_update_true(self, tmp_path):
        """Test that update functions are called when update flag is True."""
        from datetime import datetime
        import os

        args = MagicMock()
        args.update = True
        args.path = str(tmp_path)

        # Create required directories
        (tmp_path / "data" / "mlst").mkdir(parents=True)
        (tmp_path / "data" / "tox").mkdir(parents=True)
        (tmp_path / "data" / "resistance" / "Corynebacterium_diphtheriae").mkdir(parents=True)

        # Create date-stamped directory that update_database expects
        date = datetime.now().strftime('%Y-%m-%d')
        amr_db_path = tmp_path / "data" / "resistance" / date
        amr_db_path.mkdir(parents=True)

        mlst_db = ('mlst_header', str(tmp_path / 'data' / 'mlst' / 'mlst.fas'), str(tmp_path / 'st_profiles.txt'))
        tox_db = ('tox_header', str(tmp_path / 'data' / 'tox' / 'tox.fas'), str(tmp_path / 'tox_profiles.txt'))

        with patch('diphtoscan.updating_database.remove_mlst_database') as mock_remove:
            with patch('diphtoscan.updating_database.create_db') as mock_create:
                mock_create.return_value = ('path', ['atpA'])
                with patch('diphtoscan.updating_database.download_profiles_st'):
                    with patch('diphtoscan.updating_database.download_profiles_tox'):
                        with patch('diphtoscan.updating_database.download_amrfinder_database'):
                            with patch('diphtoscan.updating_database.update_amrfinderplus_db_file'):
                                with patch('diphtoscan.updating_database.complete_missing_classification'):
                                    with patch('subprocess.run') as mock_run:
                                        mock_run.return_value = MagicMock(stdout='4.0.0')

                                        update_database(args, mlst_db, tox_db)

                                        # Should call remove_mlst_database twice (for MLST and TOX)
                                        assert mock_remove.call_count == 2
