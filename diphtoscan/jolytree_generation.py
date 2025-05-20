import os
import shutil
import subprocess

def generate_jolytree(arguments): 
        print ("\nGenerating a phylogenetic tree from JolyTree \n")
        os.makedirs(arguments.outdir+"/FolderJolyTree" )
        for assembly in arguments.assemblies:
                if not os.path.exists(assembly):
                        raise FileNotFoundError(f"Assembly file {assembly} does not exist.")
                shutil.copy(assembly, arguments.outdir+"/FolderJolyTree/")
        subprocess.run(['JolyTree.sh', '-i', arguments.outdir+"/FolderJolyTree", 
                        '-b', arguments.outdir + 'jolytree', '-t', str(arguments.threads)])
        shutil.rmtree(arguments.outdir+"/FolderJolyTree/")
