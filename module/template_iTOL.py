"""
Copyright 2022 Melanie Hennart (melanie.hennart@pasteur.fr)
https://gitlab.pasteur.fr/BEBP

This file is part of diphtOscan. diphtOscan is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. diphtOscan is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with diphtOscan. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import pandas as pd 



def get_BINARY_header():   
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


def get_STRIP_header():   
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


def get_TOX_header():  
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


def writeTemplateBinary (outdir, file, column, values, colors, symbols):
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




def writeTemplateTOX (outdir, file, column):
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


def writeTemplateStrip (outdir, file, column, list_familiesRes):
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


       




     


