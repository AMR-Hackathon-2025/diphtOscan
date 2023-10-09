import os
import glob
import pandas as pd

from module.mlstBLAST import mlst_blast

def get_chromosome_mlst_results(infoMLST:tuple, contigs:str, cd_complex:bool, args) -> dict:
    chromosome_mlst_header = infoMLST[0]
    if cd_complex:
        seqs = infoMLST[1]
        database = infoMLST[2]
        chr_st, chr_st_detail, _, _ = \
             mlst_blast(seqs, database, 'no', [contigs], min_cov=args.min_coverage,
                       min_ident=args.min_identity, max_missing=3, allow_multiple=False)
        if chr_st != '0':
            chr_st = 'ST' + chr_st
        
        assert len(chromosome_mlst_header) == len(chr_st_detail)
        results = {'ST': chr_st}

    else:
        results = {'ST': "NA"}
        chr_st_detail = ['-'] * len(chromosome_mlst_header)

    results.update(dict(zip(infoMLST[0], chr_st_detail)))
    return results

def get_tox_results(infoTOX:tuple, contigs:str, args) -> dict:
    tox_header = infoTOX[0]
    seqs = infoTOX[1]
    database = infoTOX[2]
    chr_st, chr_st_detail, _, _ = \
         mlst_blast(seqs, database, 'no', [contigs], min_cov=args.min_coverage,
                   min_ident=args.min_identity, max_missing=3, allow_multiple=False)
    if chr_st != '0':
        chr_st = 'TOX' + chr_st
    
    assert len(tox_header) == len(chr_st_detail)
    results = dict(zip(infoTOX[0], chr_st_detail))
    #results = {'TOX': chr_st}
    
    #results.update(dict(zip(infoTOX[0], chr_st_detail)))
    return results

def get_chromosome_mlst_header() -> list:
    return ['atpA', 'dnaE', 'dnaK', 'fusA', 'leuA', 'odhA', 'rpoB']

def get_tox_header() -> list:
    return ['tox']

def get_virulence() -> list:
    return ['REPRESSOR','TOXIN','OTHERS_TOXIN', 
            'spuA', 'narG',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN', 'iusABCDE',
            'chtAB','htaA-hmuTUV-htaBC', 'cdtQP-sidBA-ddpABCD']
    
def get_virulence_extended() -> list: 
    return ['REPRESSOR','TOXIN','OTHERS_TOXIN', 
            'SpuA-CLUSTER', 'narIJHK',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD', 
            'ciuABCD',  'ciuEFG', 'chtAB','chtC','cdtQP-sidBA-ddpABCD','HbpA']
def delete_virulence_extended() -> list:
    return [ 'SpuA-CLUSTER','narIJHK','SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN', 
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD', 
            'ciuABCD',  'ciuEFG', 'chtAB','chtC','cdtQP-sidBA-ddpABCD','HbpA']
       
       
def is_non_zero_file(fpath:str):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def armfinder_to_table(data_resistance:pd.DataFrame) ->  pd.DataFrame:
    dico_Method = {'ALLELEX' : "",
                   'EXACTX' :  "",
                   'POINTX' : "!",
                   'BLASTX' : "*",
                   'PARTIALX' : "?",
                   'PARTIAL_CONTIG_ENDX' : "?$",
                   'INTERNAL_STOP' :  "#"}
    data_resistance['Class'] = data_resistance['Class'].fillna ('NoClass')
    Class = data_resistance['Class'].value_counts().keys()
    Strains = data_resistance['Name'].value_counts().keys()
    table = pd.DataFrame('',index=Strains, columns=Class)
    for res in data_resistance.index : 
        gene = data_resistance['Gene symbol'][res] + dico_Method[data_resistance['Method'][res]]      
        if 'tox' in data_resistance['Gene symbol'][res] :    
            if float(data_resistance['% Coverage of reference sequence'][res]) != 100.00 : 
                if (data_resistance['Method'][res] == 'BLASTX') :
                    gene = data_resistance['Gene symbol'][res] + "-NTTB?-"+str(round(100-float(data_resistance['% Coverage of reference sequence'][res])))+"%"              
                else :
                    gene = data_resistance['Gene symbol'][res] + "-NTTB" + dico_Method[data_resistance['Method'][res]]             
        if (data_resistance['Method'][res] == 'PARTIALX') or \
           (data_resistance['Method'][res] == 'PARTIAL_CONTIG_ENDX') or \
           (data_resistance['Method'][res] == 'INTERNAL_STOP') : 
               gene += "-"+str(round(100-float(data_resistance['% Coverage of reference sequence'][res])))+"%"
        
        strain = data_resistance['Name'][res]
        family = data_resistance['Class'][res]
        
        if table[family][strain] != '' :
               table[family][strain] += ";"
        table[family][strain] += gene
    return table

def get_genomic_context(outdir:str, data:pd.DataFrame) :
    d = []
    data_AMR = data[~data['Class'].isin( list(set(get_virulence_extended())| set(get_virulence())))]
    fi = open(outdir+'/distance_context.txt', 'a', encoding='utf-8')
    for contigs in data_AMR['Contig id'].value_counts().keys() :            
        table_contigs  = data_AMR[data_AMR['Contig id'] == contigs]
        
        if len(table_contigs) == 1 : 
            d.append(table_contigs['Gene symbol'].value_counts().keys()[0])
        else:  
            t = table_contigs['Gene symbol'].iloc[0]
            for i in range(0,len(table_contigs)-1):
                dis = int(table_contigs['Start'].iloc[i+1]) - int(table_contigs['Stop'].iloc[i])
                fi.write(table_contigs['Gene symbol'].iloc[i]+'\t'+table_contigs['Gene symbol'].iloc[i+1]+'\t'+str(abs(dis))+'\n')
                if abs(dis) <=  8000 :  
                    t +=  ";" + table_contigs['Gene symbol'].iloc[i+1]
                else :
                    t +=  " || " + table_contigs['Gene symbol'].iloc[i+1]                
            d.append(t)
    fi.close()        
    return " || ".join(d)

def find_resistance_db(args):
            files = [ name for name in glob.glob(args.path+'/data/resistance/*') if os.path.isdir(name) ]
            max_file = max(files, key = os.path.getctime)    
            return max_file 