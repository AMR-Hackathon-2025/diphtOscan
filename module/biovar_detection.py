import pandas as pd

def spuA(results:pd.DataFrame, arguments):       
    if "spuA" in results.columns:
        SpuA_CLUSTER = ["spuA"]
        SpuA_CLUSTER_color = ['#002b00']
        SpuA_CLUSTER_symbol = ["2"]
        writeTemplateBinary(arguments.outdir , results, "spuA", SpuA_CLUSTER, SpuA_CLUSTER_color, SpuA_CLUSTER_symbol)

        
def narG(results:pd.DataFrame, arguments):  
    if "narG" in results.columns:
        narIJHGK = ["narG"]
        narIJHGK_color = ['#f1c40f']
        narIJHGK_symbol = ["2"]
        writeTemplateBinary(arguments.outdir, results, "narG", narIJHGK, narIJHGK_color, narIJHGK_symbol)


def toxin(results:pd.DataFrame, arguments):  
    if "TOXIN" in results.columns:
        writeTemplateTOX(arguments.outdir, results, 'TOXIN')