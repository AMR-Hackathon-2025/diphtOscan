"""Tests for diphtoscan/jolytree_generation.py - JolyTree phylogenetic tree generation."""

import pytest
from unittest.mock import patch, MagicMock
import os

from diphtoscan.jolytree_generation import generate_jolytree


class TestGenerateJolytree:
    """Tests for generate_jolytree function."""

    def test_creates_folder_structure(self, tmp_path):
        """Test that JolyTree folder is created."""
        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = []
        args.threads = 4

        # Create output dir
        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                joly_folder = tmp_path / "output" / "FolderJolyTree"
                assert joly_folder.exists()

    def test_copies_assemblies(self, tmp_path):
        """Test that assembly files are copied to JolyTree folder."""
        # Create test assembly files
        assembly1 = tmp_path / "assembly1.fasta"
        assembly2 = tmp_path / "assembly2.fasta"
        assembly1.write_text(">seq1\nATGC")
        assembly2.write_text(">seq2\nGCTA")

        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(assembly1), str(assembly2)]
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                joly_folder = tmp_path / "output" / "FolderJolyTree"
                assert (joly_folder / "assembly1.fasta").exists()
                assert (joly_folder / "assembly2.fasta").exists()

    def test_calls_jolytree_command(self, tmp_path):
        """Test that JolyTree.sh command is invoked."""
        assembly = tmp_path / "assembly.fasta"
        assembly.write_text(">seq\nATGC")

        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(assembly)]
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                mock_run.assert_called_once()
                call_args = mock_run.call_args[0][0]
                assert 'JolyTree.sh' in call_args
                assert '-i' in call_args
                assert '-b' in call_args
                assert '-t' in call_args

    def test_passes_threads_argument(self, tmp_path):
        """Test that threads argument is passed to JolyTree."""
        assembly = tmp_path / "assembly.fasta"
        assembly.write_text(">seq\nATGC")

        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(assembly)]
        args.threads = 8

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                call_args = mock_run.call_args[0][0]
                # Check that threads value is in the command
                assert '8' in call_args

    def test_cleanup_after_completion(self, tmp_path):
        """Test that JolyTree folder is cleaned up after completion."""
        assembly = tmp_path / "assembly.fasta"
        assembly.write_text(">seq\nATGC")

        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(assembly)]
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run'):
            with patch('shutil.rmtree') as mock_rmtree:
                generate_jolytree(args)

                # rmtree should be called with the JolyTree folder
                mock_rmtree.assert_called_once()
                call_arg = mock_rmtree.call_args[0][0]
                assert 'FolderJolyTree' in call_arg

    def test_raises_on_missing_assembly(self, tmp_path):
        """Test that FileNotFoundError is raised for missing assembly."""
        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(tmp_path / "nonexistent.fasta")]
        args.threads = 4

        os.makedirs(args.outdir)

        with pytest.raises(FileNotFoundError):
            generate_jolytree(args)

    def test_output_path_includes_jolytree(self, tmp_path):
        """Test that output path includes 'jolytree' prefix."""
        assembly = tmp_path / "assembly.fasta"
        assembly.write_text(">seq\nATGC")

        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(assembly)]
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                call_args = mock_run.call_args[0][0]
                # Find the -b argument value
                b_index = call_args.index('-b')
                output_path = call_args[b_index + 1]
                assert 'jolytree' in output_path


class TestGenerateJolytreeEmptyAssemblies:
    """Tests for generate_jolytree with empty assemblies."""

    def test_handles_empty_assembly_list(self, tmp_path):
        """Test handling of empty assemblies list."""
        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = []
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                # Should still call JolyTree even with empty folder
                mock_run.assert_called_once()


class TestGenerateJolytreePaths:
    """Tests for path handling in generate_jolytree."""

    def test_handles_paths_with_spaces(self, tmp_path):
        """Test handling of paths with spaces."""
        space_dir = tmp_path / "path with spaces"
        space_dir.mkdir()
        assembly = space_dir / "assembly.fasta"
        assembly.write_text(">seq\nATGC")

        args = MagicMock()
        args.outdir = str(tmp_path / "output dir")
        args.assemblies = [str(assembly)]
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                # Should handle paths with spaces correctly
                mock_run.assert_called_once()

    def test_handles_absolute_paths(self, tmp_path):
        """Test handling of absolute assembly paths."""
        assembly = tmp_path / "assembly.fasta"
        assembly.write_text(">seq\nATGC")

        args = MagicMock()
        args.outdir = str(tmp_path / "output")
        args.assemblies = [str(assembly.absolute())]
        args.threads = 4

        os.makedirs(args.outdir)

        with patch('subprocess.run') as mock_run:
            with patch('shutil.rmtree'):
                generate_jolytree(args)

                joly_folder = tmp_path / "output" / "FolderJolyTree"
                assert (joly_folder / "assembly.fasta").exists()
