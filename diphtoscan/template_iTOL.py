"""
Copyright 2023 Melanie Hennart (melanie.hennart@pasteur.fr)
Copyright 2023 Martin Rethoret Pasty (martin.rethoret-pasty@pasteur.fr)
https://gitlab.pasteur.fr/BEBP

This file is part of diphtOscan. diphtOscan is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. diphtOscan is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with diphtOscan. If
not, see <http://www.gnu.org/licenses/>.
"""

import pandas as pd 


def get_BINARY_header() -> str:
    """
    Return the iTOL binary dataset header template.

    Returns
    -------
    str
        Header template for iTOL binary presence/absence datasets
        with placeholders for Title, Shapes, Labels, and Colors.
    """
    header_BINARY = """DATASET_BINARY
SEPARATOR COMMA
DATASET_LABEL,Title
COLOR,#ff0000
FIELD_SHAPES,Shapes_Binary
FIELD_LABELS,Labels_Binary
FIELD_COLORS,Colors_Binary
#=================================================================#
DATA
"""
    return header_BINARY


def get_STRIP_header() -> str:
    """
    Return the iTOL color strip dataset header template.

    Returns
    -------
    str
        Header template for iTOL color strip datasets
        with placeholder for Title.
    """
    header_STRIP = """DATASET_COLORSTRIP
SEPARATOR COMMA
DATASET_LABEL,Title
COLOR,#ff0000
SHOW_LABELS,1
LABEL_SHIFT,10
BORDER_WIDTH,1
BORDER_COLOR,#ffffff
COMPLETE_BORDER,1

#=================================================================#
DATA
"""
    return header_STRIP 


def get_TOX_header() -> str:
    """
    Return the iTOL toxin-specific binary dataset header.

    Returns
    -------
    str
        Pre-configured header for toxin visualization with
        separate columns for intact and truncated toxin.
    """
    header_TOX = """DATASET_BINARY
SEPARATOR COMMA
DATASET_LABEL,toxin
COLOR,#ff0000
FIELD_SHAPES,3,3
FIELD_LABELS,toxin,toxin truncated
FIELD_COLORS,#cc0000,#ee6500
SYMBOL_SPACING,-27
#=================================================================#
DATA
"""
    return header_TOX


def writeTemplateBinary(outdir: str, file: pd.DataFrame, column: str,
                        values: list, colors: list, symbols: list) -> pd.DataFrame:
    """
    Write an iTOL binary presence/absence dataset file.

    Parameters
    ----------
    outdir : str
        Output directory path.
    file : pd.DataFrame
        Results DataFrame with strain indices.
    column : str
        Column name from results to visualize.
    values : list
        Gene/feature names to check for presence.
    colors : list
        Hex color codes for each value.
    symbols : list
        iTOL shape codes for each value (e.g., '1', '2', '3').

    Returns
    -------
    pd.DataFrame
        Binary presence matrix (1/0) for the analyzed column.
    """
    f = open(outdir+"/"+column.replace('/','_')+".txt", 'w', encoding='utf-8')
    header_BINARY = get_BINARY_header()
    header = header_BINARY.replace("Title", column)
    header = header.replace("Shapes_Binary", ','.join(symbols))
    header = header.replace("Labels_Binary", ','.join(values))
    header = header.replace("Colors_Binary", ','.join(colors))
    f.write(header)
    data = pd.DataFrame(index = file.index, columns=values)
    for gene in values:
        data[gene] = file[column].apply(lambda x : "1" if (gene in x) else "0" )
    for strain in data.index : 
        line = strain+','+",".join(data.loc[strain].values) + "\n"
        f.write(line)  
    f.close()
    return data


def writeTemplateTOX(outdir: str, file: pd.DataFrame, column: str) -> None:
    """
    Write an iTOL toxin visualization dataset file.

    Creates a binary dataset with two columns: intact toxin and
    truncated toxin.

    Parameters
    ----------
    outdir : str
        Output directory path.
    file : pd.DataFrame
        Results DataFrame with strain indices and toxin column.
    column : str
        Column name containing toxin detection results.
    """
    f = open(outdir+'/'+column.replace('/','_')+".txt", 'w', encoding='utf-8')
    f.write(get_TOX_header())
    for strain in file.index :
        line = [strain, '-1','-1']
        if "tox" in file[column][strain] :
            if "-" in file[column][strain] : 
               line[2] = '1' 
            else : 
               line[1] = '1'                             
        line = ",".join(line) + "\n"
        f.write(line)  
    f.close()
    return 


def writeTemplateStrip(outdir: str, file: pd.DataFrame, column: str,
                       list_familiesRes: dict) -> None:
    """
    Write an iTOL color strip dataset file for AMR families.

    Parameters
    ----------
    outdir : str
        Output directory path.
    file : pd.DataFrame
        Results DataFrame with strain indices and AMR columns.
    column : str
        AMR family column name to visualize.
    list_familiesRes : dict
        Dictionary mapping family names to [absent_color, present_color].
    """
    f = open(outdir+'/'+column+".txt", 'w', encoding='utf-8')
    header_STRIP = get_STRIP_header()
    header = header_STRIP.replace("Title", column)
    f.write(header)
    for strain in file.index : 
        if file[column][strain] == "-" : 
            line = strain+","+list_familiesRes[column][0]+","+file[column][strain]+"\n"
        else :
            line = strain+","+list_familiesRes[column][1]+","+file[column][strain]+"\n"
        f.write(line)                
    f.close()
    return 

def spuA(results: pd.DataFrame, arguments) -> None:
    """
    Generate iTOL visualization file for spuA gene presence.

    Parameters
    ----------
    results : pd.DataFrame
        Analysis results DataFrame.
    arguments : argparse.Namespace
        Parsed arguments with outdir attribute.
    """
    if "spuA" in results.columns:
        SpuA_CLUSTER = ["spuA"]
        SpuA_CLUSTER_color = ['#002b00']
        SpuA_CLUSTER_symbol = ["2"]
        writeTemplateBinary(arguments.outdir , results, "spuA", SpuA_CLUSTER, SpuA_CLUSTER_color, SpuA_CLUSTER_symbol)

        
def narG(results: pd.DataFrame, arguments) -> None:
    """
    Generate iTOL visualization file for narG gene presence.

    Parameters
    ----------
    results : pd.DataFrame
        Analysis results DataFrame.
    arguments : argparse.Namespace
        Parsed arguments with outdir attribute.
    """
    if "narG" in results.columns:
        narIJHGK = ["narG"]
        narIJHGK_color = ['#f1c40f']
        narIJHGK_symbol = ["2"]
        writeTemplateBinary(arguments.outdir, results, "narG", narIJHGK, narIJHGK_color, narIJHGK_symbol)


def toxin(results: pd.DataFrame, arguments) -> None:
    """
    Generate iTOL visualization file for toxin gene status.

    Parameters
    ----------
    results : pd.DataFrame
        Analysis results DataFrame.
    arguments : argparse.Namespace
        Parsed arguments with outdir attribute.
    """
    if "TOXIN" in results.columns:
        writeTemplateTOX(arguments.outdir, results, 'TOXIN')

list_familiesRes ={'AMINOGLYCOSIDE' : ['#a6cee3', '#1f78b4'],
                   'MACROLIDE' : ['#b2df8a', '#33a02c'],
                   'PHENICOL' : ['#fb9a99', '#e31a1c'],
                   'SULFONAMIDE' : ['#fdbf6f', '#ff7f00'],
                   'TETRACYCLINE' : ['#cab2d6', "#6a3d9a"],
                   'TRIMETHOPRIM' : ['#ffff99', '#b15928'],
                   'QUATERNARY AMMONIUM' : ["#e0eaf4", "#3c6498"],
                   'BETA-LACTAM' : ["#da74da", "#9e3c8b"],
                   'QUINOLONE' : ["#b0c665", "#6c7b38"],
                   'RIFAMYCIN' : ["#bd924f", "#926114"]
                   }

def amr_families(results: pd.DataFrame, arguments) -> None:
    """
    Generate iTOL visualization files for all AMR families.

    Creates color strip files for each antimicrobial resistance family
    present in the results.

    Parameters
    ----------
    results : pd.DataFrame
        Analysis results DataFrame.
    arguments : argparse.Namespace
        Parsed arguments with outdir attribute.
    """
    for family in list_familiesRes : 
        if family in results.columns:
            writeTemplateStrip (arguments.outdir, results, family, list_familiesRes)