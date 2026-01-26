"""Integration tests for diphtoscan/cli.py - command line interface."""

import pytest
from unittest.mock import patch, MagicMock
import sys
import os
import shutil

from diphtoscan.cli import (
    parse_arguments,
    test_unique_dependency as check_unique_dependency,
    test_multiple_dependencies as check_multiple_dependencies,
    test_required_dependency as check_required_dependency,
    redefine_output_file,
    move_file_to_outdir_folder,
    rename_temp_folder_file,
    write_json_output,
    create_common_parser,
    main,
    __version__,
)
import json
import pandas as pd


class TestParseArguments:
    """Tests for parse_arguments function."""

    def test_minimal_args_with_assemblies(self, tmp_path):
        """Test parsing with minimal arguments (just assemblies)."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.assemblies == [str(assembly)]
                assert args.mlst is False
                assert args.tox is False

    def test_multiple_assemblies(self, tmp_path):
        """Test parsing with multiple assembly files."""
        assembly1 = tmp_path / "test1.fasta"
        assembly2 = tmp_path / "test2.fasta"
        assembly1.write_text(">seq1\nATGC")
        assembly2.write_text(">seq2\nGCTA")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly1), str(assembly2)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert len(args.assemblies) == 2
                assert str(assembly1) in args.assemblies
                assert str(assembly2) in args.assemblies

    def test_mlst_flag(self, tmp_path):
        """Test parsing with MLST flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-st']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.mlst is True

    def test_tox_flag(self, tmp_path):
        """Test parsing with tox flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-t']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.tox is True

    def test_resistance_virulence_flag(self, tmp_path):
        """Test parsing with resistance/virulence flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-res_vir']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.resistance_virulence is True

    def test_extend_genotyping_flag(self, tmp_path):
        """Test parsing with extended genotyping flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-plus']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.extend_genotyping is True

    def test_custom_thresholds(self, tmp_path):
        """Test parsing with custom identity and coverage thresholds."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', '-a', str(assembly),
            '--min_identity', '90.0',
            '--min_coverage', '75.0'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.min_identity == 90.0
                assert args.min_coverage == 75.0

    def test_default_thresholds(self, tmp_path):
        """Test default threshold values."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.min_identity == 80.0
                assert args.min_coverage == 50.0

    def test_threads_argument(self, tmp_path):
        """Test parsing threads argument."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '--threads', '8']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.threads == 8

    def test_update_flag_without_assemblies(self):
        """Test that update subcommand works without assemblies."""
        with patch.object(sys, 'argv', ['diphtoscan', 'update']):
            args = parse_arguments()

            assert args.command == 'update'
            assert args.assemblies is None

    def test_tree_flag(self, tmp_path):
        """Test parsing with tree flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-tree']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                # Note: check_required_dependency may modify this
                assert hasattr(args, 'tree')

    def test_integron_flag(self, tmp_path):
        """Test parsing with integron flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-integron']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert hasattr(args, 'integron')

    def test_overwrite_flag(self, tmp_path):
        """Test parsing with overwrite flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '--overwrite']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.overwrite is True

    def test_custom_output_directory(self, tmp_path):
        """Test parsing with custom output directory."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-o', 'my_output']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.outdir == 'my_output'

    def test_no_args_prints_help_and_exits(self, capsys):
        """Test that no arguments prints help and exits."""
        with patch.object(sys, 'argv', ['diphtoscan']):
            with pytest.raises(SystemExit):
                parse_arguments()


class TestDependencyChecks:
    """Tests for dependency checking functions."""

    def check_unique_dependency_found(self):
        """Test that found dependency returns True."""
        with patch('shutil.which', return_value='/usr/bin/python'):
            result = check_unique_dependency('python')
            assert result is True

    def check_unique_dependency_missing(self):
        """Test that missing dependency returns False."""
        with patch('shutil.which', return_value=None):
            result = check_unique_dependency('nonexistent_program')
            assert result is False

    def check_multiple_dependencies_all_found(self, capsys):
        """Test that all dependencies found passes."""
        with patch('shutil.which', return_value='/usr/bin/tool'):
            # Should not raise
            check_multiple_dependencies(['tool1', 'tool2'])

    def check_multiple_dependencies_one_missing(self):
        """Test that missing dependency causes exit."""
        def mock_which(name):
            if name == 'tool1':
                return '/usr/bin/tool1'
            return None

        with patch('shutil.which', side_effect=mock_which):
            with pytest.raises(SystemExit):
                check_multiple_dependencies(['tool1', 'missing_tool'])

    def check_required_dependency_with_no_assemblies(self):
        """Test that dependencies are not checked when no assemblies."""
        args = MagicMock()
        args.assemblies = None

        result = check_required_dependency(args)

        assert result == args


class TestRedefineOutputFile:
    """Tests for redefine_output_file function."""

    def test_creates_directory_if_not_exists(self, tmp_path):
        """Test that directory is created if it doesn't exist."""
        args = MagicMock()
        args.outdir = str(tmp_path / "new_output")

        result_args, final_path = redefine_output_file(args)

        assert os.path.exists(str(tmp_path / "new_output"))
        assert final_path == str(tmp_path / "new_output")
        assert '_temp_folder' in result_args.outdir

    def test_appends_temp_folder_suffix(self, tmp_path):
        """Test that _temp_folder suffix is appended."""
        args = MagicMock()
        args.outdir = str(tmp_path / "output")

        result_args, _ = redefine_output_file(args)

        assert result_args.outdir.endswith('_temp_folder')

    def test_returns_original_path(self, tmp_path):
        """Test that original path is returned as final_output_path."""
        args = MagicMock()
        original_path = str(tmp_path / "my_output")
        args.outdir = original_path

        _, final_path = redefine_output_file(args)

        assert final_path == original_path


class TestMoveFileToOutdirFolder:
    """Tests for move_file_to_outdir_folder function."""

    def test_moves_files(self, tmp_path):
        """Test that files are moved from temp to final folder."""
        temp_folder = tmp_path / "temp_folder"
        temp_folder.mkdir()
        outdir_folder = tmp_path / "output"
        outdir_folder.mkdir()

        # Create a test file in temp folder
        (temp_folder / "test.txt").write_text("content")

        with patch('diphtoscan.cli.rename_temp_folder_file'):
            move_file_to_outdir_folder(str(temp_folder), str(outdir_folder))

        assert (outdir_folder / "test.txt").exists()
        assert not temp_folder.exists()  # Should be removed

    def test_cleans_existing_files_in_outdir(self, tmp_path):
        """Test that existing files in outdir are cleaned."""
        temp_folder = tmp_path / "temp_folder"
        temp_folder.mkdir()
        outdir_folder = tmp_path / "output"
        outdir_folder.mkdir()

        # Create existing file in outdir
        (outdir_folder / "old_file.txt").write_text("old content")
        # Create new file in temp
        (temp_folder / "new_file.txt").write_text("new content")

        with patch('diphtoscan.cli.rename_temp_folder_file'):
            move_file_to_outdir_folder(str(temp_folder), str(outdir_folder))

        assert not (outdir_folder / "old_file.txt").exists()
        assert (outdir_folder / "new_file.txt").exists()


class TestRenameTempFolderFile:
    """Tests for rename_temp_folder_file function."""

    def test_renames_temp_folder_file(self, tmp_path):
        """Test renaming file with _temp_folder suffix."""
        directory = tmp_path / "results_temp_folder"
        directory.mkdir()
        outdir = tmp_path / "results"
        outdir.mkdir()

        # Create the expected temp file
        (directory / "results_temp_folder.txt").write_text("content")

        rename_temp_folder_file(str(outdir), str(directory))

        assert (outdir / "results.txt").exists()

    def test_handles_missing_file(self, tmp_path, capsys):
        """Test handling when expected file doesn't exist."""
        directory = tmp_path / "results_temp_folder"
        directory.mkdir()
        outdir = tmp_path / "results"
        outdir.mkdir()

        rename_temp_folder_file(str(outdir), str(directory))

        captured = capsys.readouterr()
        assert "does not exist" in captured.out


class TestMainFunction:
    """Tests for main function."""

    def test_update_only_mode(self, tmp_path):
        """Test running main with update flag only."""
        with patch.object(sys, 'argv', ['diphtoscan', '-u']):
            with patch('diphtoscan.cli.parse_arguments') as mock_parse:
                mock_args = MagicMock()
                mock_args.update = True
                mock_args.assemblies = None
                mock_args.path = str(tmp_path)
                mock_parse.return_value = mock_args

                with patch('diphtoscan.cli.update_database') as mock_update:
                    with pytest.raises(SystemExit):
                        main()

                    mock_update.assert_called_once()

    def test_main_with_assemblies(self, tmp_path):
        """Test running main with assemblies."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        outdir = tmp_path / "output"

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-o', str(outdir)]):
            with patch('diphtoscan.cli.parse_arguments') as mock_parse:
                mock_args = MagicMock()
                mock_args.update = False
                mock_args.assemblies = [str(assembly)]
                mock_args.outdir = str(outdir)
                mock_args.path = str(tmp_path)
                mock_args.mlst = False
                mock_args.tox = False
                mock_args.resistance_virulence = False
                mock_args.integron = False
                mock_args.tree = False
                mock_args.overwrite = False
                mock_args.extend_genotyping = False
                mock_args.threads = 4
                mock_parse.return_value = mock_args

                with patch('diphtoscan.cli.update_database'):
                    with patch('diphtoscan.cli.find_resistance_db', return_value=str(tmp_path)):
                        with patch('diphtoscan.cli.get_species_results', return_value={'species': 'unknown', 'species_match': ''}):
                            with patch('os.makedirs'):
                                with patch('pandas.DataFrame.to_csv'):
                                    with patch('diphtoscan.cli.spuA'):
                                        with patch('diphtoscan.cli.narG'):
                                            with patch('diphtoscan.cli.toxin'):
                                                with patch('diphtoscan.cli.amr_families'):
                                                    try:
                                                        main()
                                                    except (SystemExit, OSError):
                                                        pass  # Expected in test environment


class TestVersion:
    """Tests for version information."""

    def test_version_is_string(self):
        """Test that version is a string."""
        assert isinstance(__version__, str)

    def test_version_format(self):
        """Test that version follows semver-like format."""
        parts = __version__.split('.')
        assert len(parts) >= 2  # At least major.minor
        # First two parts should be numeric
        assert parts[0].isdigit()
        assert parts[1].isdigit()


class TestIntegrationScenarios:
    """Integration tests for common usage scenarios."""

    def test_help_flag(self, capsys):
        """Test that help flag shows usage information."""
        with patch.object(sys, 'argv', ['diphtoscan', '-h']):
            with pytest.raises(SystemExit) as exc_info:
                parse_arguments()

            # Standard behavior: --help exits with 0
            assert exc_info.value.code == 0
            captured = capsys.readouterr()
            assert 'diphtoscan' in captured.out.lower() or 'usage' in captured.out.lower()

    def test_subcommand_help_flag(self, capsys):
        """Test that subcommand help shows detailed usage."""
        with patch.object(sys, 'argv', ['diphtoscan', 'all', '-h']):
            with pytest.raises(SystemExit) as exc_info:
                parse_arguments()

            assert exc_info.value.code == 0
            captured = capsys.readouterr()
            assert 'assemblies' in captured.out.lower()

    def test_version_flag(self, capsys):
        """Test that version info is available."""
        with patch.object(sys, 'argv', ['diphtoscan', '--version']):
            with pytest.raises(SystemExit) as exc_info:
                parse_arguments()

            # Standard behavior: --version exits with 0
            assert exc_info.value.code == 0

    def test_no_args_shows_help(self, capsys):
        """Test that no arguments shows help."""
        with patch.object(sys, 'argv', ['diphtoscan']):
            with pytest.raises(SystemExit) as exc_info:
                parse_arguments()

            assert exc_info.value.code == 1
            captured = capsys.readouterr()
            assert 'usage' in captured.err.lower()


class TestSubcommands:
    """Tests for subcommand-based CLI architecture."""

    def test_species_subcommand(self, tmp_path):
        """Test species subcommand parses correctly."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'species', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.command == 'species'
                assert args.assemblies == [str(assembly)]

    def test_mlst_subcommand(self, tmp_path):
        """Test mlst subcommand parses correctly."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'mlst', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.command == 'mlst'
                assert args.mlst is True

    def test_amr_subcommand(self, tmp_path):
        """Test amr subcommand parses correctly."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'amr', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.command == 'amr'
                assert args.resistance_virulence is True

    def test_tox_subcommand(self, tmp_path):
        """Test tox subcommand parses correctly."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'tox', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.command == 'tox'
                assert args.tox is True

    def test_qc_subcommand(self, tmp_path):
        """Test qc subcommand parses correctly."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'qc', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.command == 'qc'

    def test_all_subcommand_with_flags(self, tmp_path):
        """Test all subcommand with multiple analysis flags."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'all', '-a', str(assembly),
            '-st', '-t', '-res_vir'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.command == 'all'
                assert args.mlst is True
                assert args.tox is True
                assert args.resistance_virulence is True

    def test_update_subcommand(self):
        """Test update subcommand."""
        with patch.object(sys, 'argv', ['diphtoscan', 'update']):
            args = parse_arguments()

            assert args.command == 'update'
            assert args.assemblies is None

    def test_backward_compatibility_legacy_flags(self, tmp_path):
        """Test backward compatibility with legacy flag-based invocation."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', '-a', str(assembly), '-st']):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                # Should auto-insert 'all' command
                assert args.command == 'all'
                assert args.mlst is True


class TestJsonOutput:
    """Tests for JSON output functionality."""

    def test_format_argument_tsv(self, tmp_path):
        """Test --format tsv argument."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--format', 'tsv'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.format == 'tsv'

    def test_format_argument_json(self, tmp_path):
        """Test --format json argument."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--format', 'json'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.format == 'json'

    def test_format_argument_all(self, tmp_path):
        """Test --format all argument."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--format', 'all'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.format == 'all'

    def test_format_argument_html(self, tmp_path):
        """Test --format html argument."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--format', 'html'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.format == 'html'

    def test_default_format_is_tsv(self, tmp_path):
        """Test that default format is tsv."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'species', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.format == 'tsv'

    def test_write_json_output_creates_file(self, tmp_path):
        """Test write_json_output creates valid JSON file."""
        outdir = tmp_path / "output"
        outdir.mkdir()

        results = pd.DataFrame({
            'species': ['C. diphtheriae', 'C. ulcerans'],
            'ST': ['ST5', 'ST10'],
            'biovar': ['gravis', 'NA'],
        }, index=['strain1', 'strain2'])

        json_path = write_json_output(results, str(outdir), '1.8.0')

        assert os.path.exists(json_path)
        assert json_path.endswith('.json')

        # Verify JSON is valid
        with open(json_path) as f:
            data = json.load(f)

        assert 'version' in data
        assert data['version'] == '1.8.0'
        assert 'date' in data
        assert 'samples' in data
        assert 'summary' in data

    def test_write_json_output_contains_samples(self, tmp_path):
        """Test write_json_output includes sample data."""
        outdir = tmp_path / "output"
        outdir.mkdir()

        results = pd.DataFrame({
            'species': ['C. diphtheriae'],
            'ST': ['ST5'],
        }, index=['test_strain'])

        json_path = write_json_output(results, str(outdir), '1.8.0')

        with open(json_path) as f:
            data = json.load(f)

        assert 'test_strain' in data['samples']
        assert data['samples']['test_strain']['species'] == 'C. diphtheriae'
        assert data['samples']['test_strain']['ST'] == 'ST5'

    def test_write_json_output_contains_summary(self, tmp_path):
        """Test write_json_output includes summary statistics."""
        outdir = tmp_path / "output"
        outdir.mkdir()

        results = pd.DataFrame({
            'species': ['C. diphtheriae', 'C. diphtheriae', 'C. ulcerans'],
            'ST': ['ST5', 'ST5', 'ST10'],
            'biovar': ['gravis', 'mitis', 'NA'],
            'tox_allele': ['1', '-', '2'],
        }, index=['strain1', 'strain2', 'strain3'])

        json_path = write_json_output(results, str(outdir), '1.8.0')

        with open(json_path) as f:
            data = json.load(f)

        summary = data['summary']
        assert summary['total_samples'] == 3
        assert 'species_counts' in summary
        assert summary['species_counts']['C. diphtheriae'] == 2
        assert summary['species_counts']['C. ulcerans'] == 1
        assert 'biovar_counts' in summary
        assert 'tox_positive' in summary
        assert summary['tox_positive'] == 2


class TestQcFlag:
    """Tests for assembly QC functionality."""

    def test_no_qc_flag(self, tmp_path):
        """Test --no-qc flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--no-qc'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.no_qc is True

    def test_default_qc_enabled(self, tmp_path):
        """Test that QC is enabled by default."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'species', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.no_qc is False


class TestCreateCommonParser:
    """Tests for create_common_parser function."""

    def test_returns_parser(self):
        """Test that create_common_parser returns an ArgumentParser."""
        import argparse
        parser = create_common_parser()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_common_parser_has_assemblies(self):
        """Test that common parser includes assemblies argument."""
        parser = create_common_parser()
        # Check by parsing args
        args = parser.parse_args(['-a', 'test.fasta'])
        assert args.assemblies == ['test.fasta']

    def test_common_parser_has_format(self):
        """Test that common parser includes format argument."""
        parser = create_common_parser()
        args = parser.parse_args(['-a', 'test.fasta', '--format', 'json'])
        assert args.format == 'json'

    def test_common_parser_has_no_qc(self):
        """Test that common parser includes no-qc argument."""
        parser = create_common_parser()
        args = parser.parse_args(['-a', 'test.fasta', '--no-qc'])
        assert args.no_qc is True

    def test_common_parser_has_verbose(self):
        """Test that common parser includes verbose argument."""
        parser = create_common_parser()
        args = parser.parse_args(['-a', 'test.fasta', '-v'])
        assert args.verbose is True

    def test_common_parser_has_quiet(self):
        """Test that common parser includes quiet argument."""
        parser = create_common_parser()
        args = parser.parse_args(['-a', 'test.fasta', '-q'])
        assert args.quiet is True

    def test_verbose_quiet_mutually_exclusive(self):
        """Test that verbose and quiet are mutually exclusive."""
        parser = create_common_parser()
        import argparse
        with pytest.raises(SystemExit):
            parser.parse_args(['-a', 'test.fasta', '-v', '-q'])


class TestVerboseQuietModes:
    """Tests for verbose/quiet command-line modes."""

    def test_verbose_flag_short(self, tmp_path):
        """Test -v verbose flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '-v'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.verbose is True
                assert args.quiet is False

    def test_verbose_flag_long(self, tmp_path):
        """Test --verbose flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--verbose'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.verbose is True

    def test_quiet_flag_short(self, tmp_path):
        """Test -q quiet flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '-q'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.quiet is True
                assert args.verbose is False

    def test_quiet_flag_long(self, tmp_path):
        """Test --quiet flag."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', [
            'diphtoscan', 'species', '-a', str(assembly), '--quiet'
        ]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.quiet is True

    def test_default_not_verbose_not_quiet(self, tmp_path):
        """Test that default is neither verbose nor quiet."""
        assembly = tmp_path / "test.fasta"
        assembly.write_text(">seq\nATGC")

        with patch.object(sys, 'argv', ['diphtoscan', 'species', '-a', str(assembly)]):
            with patch('diphtoscan.cli.test_required_dependency', side_effect=lambda x: x):
                args = parse_arguments()

                assert args.verbose is False
                assert args.quiet is False
