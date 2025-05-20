import os

def generate_jolytree(arguments): 
        print ("\nGenerating a phylogenetic tree from JolyTree \n")
        os.makedirs(arguments.outdir+"/FolderJolyTree" )
        os.system("cp "+ " ".join(arguments.assemblies) + " " + arguments.outdir+"/FolderJolyTree/")   
        os.system('bash ' + arguments.path + '/script/JolyTree/JolyTree.sh -i '+ arguments.outdir+ "/FolderJolyTree -b "+ arguments.outdir+"/"+arguments.outdir + ".jolytree -t " + str(arguments.threads) + " > "
                  + arguments.outdir+"/"+arguments.outdir + ".jolytree.log")
        os.system("rm -R "+ arguments.outdir+"/FolderJolyTree/")