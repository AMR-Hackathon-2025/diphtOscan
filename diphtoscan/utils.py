import os
import glob
import subprocess

import pandas as pd

from .mlstBLAST import mlst_blast


def find_amrfinderplus_version() -> str:
    """
    Detect the major version of AMRFinderPlus installed on the system.

    Returns
    -------
    str
        Major version number as a string (e.g., '3' or '4').

    Notes
    -----
    This is used to handle differences in database file formats and
    column names between AMRFinderPlus v3 and v4.
    """
    amrfinderplus_version = subprocess.run(['amrfinder', '--version'], capture_output=True, text=True).stdout.split('.')[0]
    return amrfinderplus_version


def find_resistance_db(args) -> str:
    """
    Find the most recent AMRFinderPlus database directory.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments containing the path attribute for diphtOscan installation.

    Returns
    -------
    str
        Path to the most recently created resistance database directory.

    Notes
    -----
    Database directories are dated (YYYY-MM-DD format) and this function
    returns the one with the most recent creation time.
    """
    files = [ name for name in glob.glob(args.path+'/data/resistance/*') if os.path.isdir(name) ]
    max_file = max(files, key = os.path.getctime)
    return max_file


def is_non_zero_file(fpath: str) -> bool:
    """
    Check if a file exists and has non-zero size.

    Parameters
    ----------
    fpath : str
        Path to the file to check.

    Returns
    -------
    bool
        True if file exists and has size > 0 bytes, False otherwise.
    """
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_chromosome_mlst_header() -> list:
    """
    Return the list of MLST locus names for C. diphtheriae.

    Returns
    -------
    list
        List of 7 MLST housekeeping gene names.
    """
    return ['atpA', 'dnaE', 'dnaK', 'fusA', 'leuA', 'odhA', 'rpoB']


def get_tox_header() -> list:
    """
    Return the header for tox allele results.

    Returns
    -------
    list
        Single-element list containing 'tox_allele'.
    """
    return ['tox_allele']


def get_virulence() -> list:
    """
    Return the list of standard virulence gene categories.

    Returns
    -------
    list
        List of virulence factor categories for basic screening.
    """
    return ['REPRESSOR','TOXIN','OTHER_TOXINS',
            'spuA', 'narG',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN',
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN', 'iusABCDE',
            'chtAB','htaA-hmuTUV-htaBC', 'cdtQP-sidBA-ddpABCD']


def get_virulence_extended() -> list:
    """
    Return the extended list of virulence gene categories.

    Used when --extend_genotyping (-plus) flag is enabled.

    Returns
    -------
    list
        Comprehensive list of all virulence factor categories including
        additional iron uptake systems and pili gene clusters.
    """
    return ['REPRESSOR','TOXIN','OTHER_TOXINS',
            'SpuA-CLUSTER', 'narIJHK',
            'SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN',
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD',
            'ciuABCD',  'ciuEFG', 'chtAB','chtC','cdtQP-sidBA-ddpABCD','HbpA']


def delete_virulence_extended() -> list:
    """
    Return categories to exclude from results when extended genotyping is disabled.

    Returns
    -------
    list
        Categories that are only included with --extend_genotyping flag.
    """
    return [ 'SpuA-CLUSTER','narIJHK','SpaA-type_pili_diphtheriae', 'SpaD-type_pili_diphtheriae',
            'SpaH-type_pili_diphtheriae', 'SapADE_diphtheriae',
            'VIRULENCE/ADHESIN',
            'irp1ABCD','irp2ABCDEFGHI', 'irp2JKLMN','irp6ABC', 'iusABCDE','iutABCDE',
            'htaA-hmuTUV-htaBC','hmuO','frgCBAD',
            'ciuABCD',  'ciuEFG', 'chtAB','chtC','cdtQP-sidBA-ddpABCD','HbpA']


def get_chromosome_mlst_results(infoMLST: tuple, contigs: str, cd_complex: bool, args) -> dict:
    """
    Perform MLST typing on an assembly and return results.

    Parameters
    ----------
    infoMLST : tuple
        Tuple of (header_list, sequence_database_path, profile_database_path).
    contigs : str
        Path to the assembly FASTA file.
    cd_complex : bool
        Whether the species is in the C. diphtheriae species complex.
        MLST is only performed if True.
    args : argparse.Namespace
        Parsed arguments containing min_coverage and min_identity thresholds.

    Returns
    -------
    dict
        Dictionary with 'ST' key and individual allele keys.
        ST format: 'ST{number}', 'ST{number}-{n}LV', '0', or 'NA'.
    """
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


def get_tox_results(infoTOX: tuple, contigs: str, args) -> dict:
    """
    Detect and type the tox gene in an assembly.

    Parameters
    ----------
    infoTOX : tuple
        Tuple of (header_list, sequence_database_path, profile_database_path).
    contigs : str
        Path to the assembly FASTA file.
    args : argparse.Namespace
        Parsed arguments containing min_coverage and min_identity thresholds.

    Returns
    -------
    dict
        Dictionary with tox allele detection results.
    """
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

def is_contig_edge(data_resistance: pd.DataFrame) -> bool:
    """
    Determine if a partial gene match is truncated at a contig boundary.

    Checks whether the missing portion of a gene could extend beyond the
    start or end of the contig it's located on.

    Parameters
    ----------
    data_resistance : pd.DataFrame
        Single row from AMRFinderPlus results containing:
        - 'Reference sequence length': Expected protein length
        - 'Start': Start position on contig
        - 'Stop': Stop position on contig
        - 'File': Path to assembly file
        - 'Contig id': Contig identifier

    Returns
    -------
    bool
        True if the gene appears truncated at a contig boundary,
        False otherwise.

    Notes
    -----
    Reference sequence length is in amino acids; multiplied by 3 for
    nucleotide comparison. Returns True if the missing portion of
    the gene would extend past position 0 or the contig length.
    """
    len_seq_ref = int(data_resistance['Reference sequence length'])*3
    pos_start = int(data_resistance['Start'])
    pos_stop = int(data_resistance['Stop'])
    len_seq_found = pos_stop - (pos_start-1)

    if len_seq_found < len_seq_ref :
        missing_nucleotides = len_seq_ref - len_seq_found
        over_start = (pos_start-missing_nucleotides) < 0
        over_stop = (find_len_contig(data_resistance['File'], data_resistance['Contig id']) - (pos_stop + missing_nucleotides)) < 0

        if over_start or over_stop :
            return True

    return False
                

def find_len_contig(file:str, contig :str):
    """Finds and returns the length of a specific contig in a FASTA file.

    :param file: Path to the FASTA file.
    :param contig: Contig number.
    :return: Length of the specified contig or None if not found.
    """
    with open(file, 'r') as fichier:
        line = fichier.readline()
        while line:
            if line.startswith('>' + contig):
                length = 0
                line = fichier.readline()
                while line and not line.startswith('>'):
                    length += len(line.strip())
                    line = fichier.readline()
                return length
            else:
                line = fichier.readline()
    return None #TODO to change  


def armfinder_to_table(data_resistance: pd.DataFrame) -> pd.DataFrame:
    """
    Convert AMRFinderPlus output to a summary table by gene class.

    Transforms raw AMRFinderPlus results into a strain-by-class matrix
    with annotated gene names indicating match quality and coverage.

    Parameters
    ----------
    data_resistance : pd.DataFrame
        Combined AMRFinderPlus output for all strains with columns including
        'Name', 'Class', 'Method', gene symbol, coverage, and positions.

    Returns
    -------
    pd.DataFrame
        Matrix with strains as rows and gene classes as columns.
        Cell values contain semicolon-separated gene names with annotations:
        - '*': BLAST match (not exact)
        - '!': Point mutation
        - '?': Partial match
        - '#': Internal stop codon
        - '_end_of_contig': Truncated at contig boundary
        - '-NTTB': Non-Toxigenic Tox-Bearing (for tox gene)
        - '-X.X%': Percentage of missing coverage

    Notes
    -----
    Handles differences between AMRFinderPlus v3 and v4 column names.
    NTTB prediction is only applied when tox gene has < 100% coverage
    and is NOT truncated at a contig boundary.
    """
    amrfinderplus_version = find_amrfinderplus_version()
    if amrfinderplus_version == '3':
        coverage_key = '% Coverage of reference sequence'
        gene_symbol_key = 'Gene symbol'
    else:
        coverage_key = '% Coverage of reference'
        gene_symbol_key = 'Element symbol'

    dico_Method = {'ALLELEX' : "",
                   'EXACTX' :  "",
                   'POINTX' : "!",
                   'BLASTX' : "*",
                   'PARTIALX' : "?",
                   'PARTIAL_CONTIG_ENDX' : "_end_of_contig", #The PARTIAL_CONTIG_ENDX method is only attributedd when the start or end position of the sequence being searched coincides exactly with the start or end of the contig.
                   'CTRL_CONTIG_END' : "_end_of_contig",
                   'INTERNAL_STOP' :  "#"}
    
    avoid_NTTB_prediction = ['PARTIAL_CONTIG_ENDX',
                             'CTRL_CONTIG_END']

    data_resistance['Class'] = data_resistance['Class'].fillna ('NoClass')
    Class = data_resistance['Class'].value_counts().keys()
    Strains = data_resistance['Name'].value_counts().keys()
    table = pd.DataFrame('',index=Strains, columns=Class)

    for res in data_resistance.index :
        gene = data_resistance[gene_symbol_key][res] + dico_Method[data_resistance['Method'][res]]
        # Search for certain cases of interruption due to a contig end that AMRfinder is unable to find.    
        if is_contig_edge(data_resistance.iloc[res]) : 
            data_resistance.loc[res, 'Method'] = "CTRL_CONTIG_END" 

        if ('tox' in data_resistance[gene_symbol_key][res]) and \
           (float(data_resistance[coverage_key][res]) != 100.00) and \
           (data_resistance['Method'][res] not in avoid_NTTB_prediction) :
                gene = data_resistance[gene_symbol_key][res] + "-NTTB"             

        # For all methods where coverage can be < 100%, display the %age of missing coverage
        if (data_resistance['Method'][res] == 'PARTIALX') or \
           (data_resistance['Method'][res] == 'BLASTX') or \
           (data_resistance['Method'][res] == 'PARTIAL_CONTIG_ENDX') or \
           (data_resistance['Method'][res] == 'CTRL_CONTIG_END') or \
           (data_resistance['Method'][res] == 'INTERNAL_STOP') :
                missing_coverage = round(100-float(data_resistance[coverage_key][res]),1)
                if (100 - missing_coverage) < 100 :
                    gene = f"{gene}-{missing_coverage}%" 
        strain = data_resistance['Name'][res]
        family = data_resistance['Class'][res]
        
        if table[family][strain] != '' :
               table.loc[strain, family] += ";"
        table.loc[strain, family] += gene
    return table


def get_genomic_context(outdir: str, data: pd.DataFrame) -> str:
    """
    Analyze and report genomic context of resistance genes.

    Groups AMR genes by their proximity on contigs. Genes within 8000 bp
    are considered co-located (potentially on the same mobile element).

    Parameters
    ----------
    outdir : str
        Output directory path for writing the distance_context.txt file.
    data : pd.DataFrame
        AMRFinderPlus results for a single strain.

    Returns
    -------
    str
        Formatted string showing gene co-location patterns.
        Genes within 8000 bp are joined by ';'.
        Genes on different contigs or > 8000 bp apart are joined by ' || '.

    Notes
    -----
    Only AMR genes are analyzed (virulence genes are excluded).
    Writes pairwise distances to distance_context.txt for reference.
    The 8000 bp threshold represents typical mobile genetic element sizes.
    """
    amrfinderplus_version = find_amrfinderplus_version()
    if amrfinderplus_version == '3':
        gene_symbol_key = 'Gene symbol'
    else:
        gene_symbol_key = 'Element symbol'

    d = []
    data_AMR = data[~data['Class'].isin( list(set(get_virulence_extended())| set(get_virulence())))]
    fi = open(outdir+'/distance_context.txt', 'a', encoding='utf-8')
    for contigs in data_AMR['Contig id'].value_counts().keys() :
        table_contigs  = data_AMR[data_AMR['Contig id'] == contigs]

        if len(table_contigs) == 1 :
            d.append(table_contigs[gene_symbol_key].value_counts().keys()[0])
        else:
            t = table_contigs[gene_symbol_key].iloc[0]
            for i in range(0,len(table_contigs)-1):
                dis = int(table_contigs['Start'].iloc[i+1]) - int(table_contigs['Stop'].iloc[i])
                fi.write(table_contigs[gene_symbol_key].iloc[i]+'\t'+table_contigs[gene_symbol_key].iloc[i+1]+'\t'+str(abs(dis))+'\n')
                # 8000 bp threshold: genes within this distance are considered co-located
                # (potentially on the same mobile genetic element or gene cluster)
                if abs(dis) <=  8000 :
                    t +=  ";" + table_contigs[gene_symbol_key].iloc[i+1]
                else :
                    t +=  " || " + table_contigs[gene_symbol_key].iloc[i+1]
            d.append(t)
    fi.close()
    return " || ".join(d)


def determine_biovar(results_dict: dict, species: str = None) -> str:
    """
    Determine the biovar of C. diphtheriae based on spuA and narG gene presence.

    Biovars are historically determined by biochemical tests, but can be
    predicted genomically based on the presence of key marker genes:
    - spuA: Spermidine/putrescine utilization (positive in gravis)
    - narG: Nitrate reductase catalytic subunit (nitrate reduction capability)

    Parameters
    ----------
    results_dict : dict
        Dictionary containing gene detection results. Should include
        keys for resistance/virulence results with gene names.
    species : str, optional
        Species name. Biovar is only determined for C. diphtheriae.

    Returns
    -------
    str
        Biovar classification:
        - 'gravis': spuA positive AND nitrate reductase positive
        - 'mitis': spuA negative AND nitrate reductase negative
        - 'belfantii': spuA negative AND nitrate reductase positive
        - 'intermediate': spuA positive AND nitrate reductase negative (rare)
        - 'NA': Not applicable (non-C. diphtheriae species)
        - 'unknown': Cannot be determined (missing gene data)

    Notes
    -----
    Traditional biovar determination:
    - gravis: ferments starch, reduces nitrate
    - mitis: does not ferment starch, does not reduce nitrate
    - belfantii: does not ferment starch, reduces nitrate (formerly intermedius)

    The spuA gene cluster is associated with starch fermentation capability,
    while narG is essential for nitrate reduction.
    """
    # Biovar is only applicable to C. diphtheriae
    if species and 'diphtheriae' not in species.lower():
        return 'NA'

    # Look for spuA and narG in the results
    spuA_present = False
    narG_present = False

    # Check various possible locations for gene results
    for key, value in results_dict.items():
        if value is None or value == '-' or value == '':
            continue

        value_str = str(value).lower()

        # Check for spuA (may be in 'spuA' column or within a gene list)
        if 'spua' in key.lower() or 'spua-cluster' in key.lower():
            if value != '-' and value != '' and 'spua' in value_str:
                spuA_present = True
        elif 'spua' in value_str and key not in ['species', 'ST']:
            spuA_present = True

        # Check for narG (may be in 'narG' column or within a gene list)
        if 'narg' in key.lower() or 'narijhk' in key.lower():
            if value != '-' and value != '' and ('narg' in value_str or 'nari' in value_str or 'narj' in value_str or 'nark' in value_str or 'narh' in value_str):
                narG_present = True
        elif 'narg' in value_str and key not in ['species', 'ST']:
            narG_present = True

    # Determine biovar based on gene presence
    if spuA_present and narG_present:
        return 'gravis'
    elif not spuA_present and not narG_present:
        return 'mitis'
    elif not spuA_present and narG_present:
        return 'belfantii'
    elif spuA_present and not narG_present:
        return 'intermediate'

    return 'unknown'
