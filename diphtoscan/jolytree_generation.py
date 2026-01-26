import os
import shutil
import subprocess


def generate_jolytree(arguments):
    """
    Generate a phylogenetic tree using JolyTree.

    Copies input assemblies to a temporary folder, runs JolyTree,
    and cleans up the temporary files.

    Parameters
    ----------
    arguments : argparse.Namespace
        Parsed arguments containing:
        - assemblies: List of assembly file paths
        - outdir: Output directory path
        - threads: Number of processing threads

    Raises
    ------
    FileNotFoundError
        If any assembly file does not exist.

    Notes
    -----
    Requires at least 4 assemblies for meaningful tree construction.
    JolyTree and its dependencies (gawk, fastme, REQ) must be in PATH.
    Output tree is written to {outdir}jolytree.* files.
    """
    print("\nGenerating a phylogenetic tree from JolyTree \n")
    os.makedirs(arguments.outdir + "/FolderJolyTree")
    for assembly in arguments.assemblies:
        if not os.path.exists(assembly):
            raise FileNotFoundError(f"Assembly file {assembly} does not exist.")
        shutil.copy(assembly, arguments.outdir + "/FolderJolyTree/")
    subprocess.run(['JolyTree.sh', '-i', arguments.outdir + "/FolderJolyTree",
                    '-b', arguments.outdir + 'jolytree', '-t', str(arguments.threads)])
    shutil.rmtree(arguments.outdir + "/FolderJolyTree/")
