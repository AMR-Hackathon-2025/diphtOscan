"""
Copyright 2020 Kat Holt
Copyright 2020 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

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

import os


def get_species_results(contigs:str, folder:str, threads:str) -> dict:
    species, species_hit_strength = get_corynebacterium_species(contigs, folder, threads)
    return {'species': species,
            'species_match': species_hit_strength}


def get_corynebacterium_species(contigs:str, folder:str, threads:str) -> tuple:
    f = os.popen('mash dist '+folder+'/species_mash_sketches.msh -p '+ threads + ' ' + contigs)

    best_species = None
    best_distance = 1.0

    for line in f:
        line_parts = line.split('\t')
        reference = line_parts[0]
        if len(line_parts) < 4:
            continue
        species = reference.split('/')[1]
        distance = float(line_parts[2])

        # Fix up the species name formatting a bit.

        species = species.replace('C.', 'C. ')
        species = species.split('_')[0]
        species = species.replace('pseudotub', 'pseudotuberculosis')
        
        if distance < best_distance:
            best_distance = distance
            best_species = species

    f.close()

    if best_distance <= 0.05:
        return best_species, 'strong'
    elif best_distance <= 0.1:
        return best_species, 'weak'
    else:
        return 'unknown', ''


def is_cd_complex(results:dict) -> bool:
    """
    Returns True if the species call is in the Cd-complex, otherwise false.
    """
    assert 'species' in results
    species = results['species']
    if species.startswith('C. diphtheriae'):
        return True
    if species.startswith('C. belfantii'):
        return True
    if species.startswith('C. rouxii'):
        return True
    if species.startswith('C. ulcerans'):
        return True
    if species.startswith('C. pseudotuberculosis'):
        return True
    if species.startswith('C. ramonii'):
        return True
    return False
