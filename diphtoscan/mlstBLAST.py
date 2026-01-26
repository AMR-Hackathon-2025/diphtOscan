#!/usr/bin/env python3
"""
Reports best match
  * :   best matching allele is not precise match
 -nLV : best matching ST is n-locus variant

If an annotation column is provided (such as clonal complex) in the final column of the profiles
file, this annotation will be reported in column 2 of the output table.

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

import collections
import os
import re

from typing import List
from .blastn import run_blastn, BlastHit
from .truncation import truncation_check


def mlst_blast(seqs: str,
               database: str,
               info_arg: str,
               assemblies: list,
               min_cov: float,
               min_ident: float,
               max_missing: int,
               check_for_truncation: bool = False,
               report_incomplete: bool = False,
               allow_multiple: bool = False,
               min_gene_count: int = None,
               unknown_group_name: str = None,
               min_spurious_cov: float = None,
               min_spurious_ident: float = None
               ) -> tuple:
    """
    Perform MLST typing using BLAST-based allele matching.

    Identifies sequence type by matching assembly sequences against
    an allele database and looking up the ST profile.

    Parameters
    ----------
    seqs : str
        Path to the combined allele FASTA database.
    database : str
        Path to the ST profile definitions file.
    info_arg : str
        'yes' to include info column from profiles, 'no' otherwise.
    assemblies : list
        List containing assembly file path (uses first element).
    min_cov : float
        Minimum coverage threshold for allele matches (0-100).
    min_ident : float
        Minimum identity threshold for allele matches (0-100).
    max_missing : int
        Maximum number of missing/mismatched loci to attempt ST calling.
    check_for_truncation : bool, optional
        Whether to check for truncated genes (default: False).
    report_incomplete : bool, optional
        Whether to add '(incomplete)' to info for missing loci (default: False).
    allow_multiple : bool, optional
        Whether to allow multiple ST calls from different contigs (default: False).
    min_gene_count : int, optional
        Minimum genes required for unknown group assignment (default: None).
    unknown_group_name : str, optional
        Name to assign when minimum gene count met but no ST (default: None).
    min_spurious_cov : float, optional
        Coverage threshold for tracking spurious hits (default: None).
    min_spurious_ident : float, optional
        Identity threshold for tracking spurious hits (default: None).

    Returns
    -------
    tuple
        (final_call, final_alleles, final_info, spurious_hits) where:
        - final_call: ST number (with -nLV suffix if variant), '0' if no call
        - final_alleles: List of allele calls per locus (with * for inexact)
        - final_info: Associated info from profile database
        - spurious_hits: List of low-quality hits (if min_spurious_* set)

    Notes
    -----
    ST calling requires at least half of the loci (len(header)/2) to have
    exact allele matches. Imprecise matches are marked with '*'.
    """
    st_names, alleles_to_st, st_to_info, header = load_st_database(database, info_arg)

    # Required exact matches = half of total loci (rounded down)
    # For 7-locus MLST: requires at least 3 exact allele matches
    # This threshold balances confidence in ST assignment with tolerance
    # for sequencing errors or novel alleles
    required_exact_matches = int(len(header) / 2)

    contigs = assemblies[0]
    _, filename = os.path.split(contigs)
    name, _ = os.path.splitext(filename)

    if min_spurious_cov is not None:
        hits = run_blastn(seqs, contigs, min_spurious_cov, min_spurious_ident)
        num_hits_before = len(hits)
        spurious_hits = [h for h in hits
                         if h.ref_cov * 100 < min_cov or h.pcid * 100 < min_ident]
        hits = [h for h in hits if h.ref_cov * 100 >= min_cov and h.pcid * 100 >= min_ident]
        assert len(hits) + len(spurious_hits) == num_hits_before
    else:
        hits = run_blastn(seqs, contigs, min_cov, min_ident)
        spurious_hits = None

    final_call = ''
    final_alleles = [''] * len(header)
    final_info = ''

    # See if we have any loci for which there are multiple genes.
    hits_one_per_locus = keep_only_one_hit_per_locus(hits)
    any_multiple_hits = len(hits_one_per_locus) < len(hits)

    # If we have multiple hits and allow multiple MLST calls, then we cluster the hits by contig
    # before continuing.
    if any_multiple_hits and allow_multiple:
        hit_groups = cluster_hits_by_contig(hits)

    # If we aren't considering multiple MLST calls, then all the hits go into a single cluster.
    else:
        hit_groups = [hits_one_per_locus]

    for hit_group in hit_groups:
        call, alleles, info = \
            call_one_st(hit_group, header, check_for_truncation, max_missing, alleles_to_st,
                        required_exact_matches, info_arg, st_to_info, report_incomplete,
                        min_gene_count, unknown_group_name)
        final_call = add_to_string(final_call, call)
        final_alleles = add_to_strings(final_alleles, alleles)
        final_info = add_to_string(final_info, info)

    if spurious_hits is not None:
        spurious_hits = process_spurious_hits(spurious_hits)

    return final_call, final_alleles, final_info, spurious_hits


def call_one_st(hits: List[BlastHit],
                header: List[str],
                check_for_truncation: bool,
                max_missing: int,
                alleles_to_st: dict,
                required_exact_matches: int,
                info_arg: str,
                st_to_info: dict,
                report_incomplete: bool,
                min_gene_count: int,
                unknown_group_name: str) -> tuple:
    """
    Call a single sequence type from a set of BLAST hits.

    Parameters
    ----------
    hits : List[BlastHit]
        BLAST hits for one contig cluster or all hits.
    header : List[str]
        List of locus names in the MLST scheme.
    check_for_truncation : bool
        Whether to check for amino acid level truncations.
    max_missing : int
        Maximum mismatching loci to still attempt ST lookup.
    alleles_to_st : dict
        Mapping from allele combination string to ST number.
    required_exact_matches : int
        Minimum exact allele matches required for ST assignment.
    info_arg : str
        'yes' to include info from profiles, 'no' otherwise.
    st_to_info : dict
        Mapping from ST to additional info (e.g., clonal complex).
    report_incomplete : bool
        Whether to mark incomplete calls in info string.
    min_gene_count : int
        Minimum genes for unknown group assignment.
    unknown_group_name : str
        Group name when enough genes but unknown ST.

    Returns
    -------
    tuple
        (bst, best_st_annotated, info_final) where:
        - bst: ST number (with -nLV suffix), '0' if no call
        - best_st_annotated: List of allele numbers with annotations
        - info_final: Info string from profile database
    """
    best_alleles = get_best_allele_per_locus(hits, check_for_truncation)

    best_st = []
    best_st_annotated = []

    mismatch_loci_including_snps, missing_loci = 0, 0

    for locus in header:
        if locus in best_alleles:
            allele = best_alleles[locus]

            # Remove * (inexact match) and truncation percentages
            allele_number = allele.replace('*', '')
            allele_number = re.sub(r'-\d+%', '', allele_number)

            if '*' in allele:
                mismatch_loci_including_snps += 1
            best_st.append(allele_number)
            best_st_annotated.append(allele)  # will still have character if imperfect match
        else:
            best_st.append('-')
            best_st_annotated.append('-')
            mismatch_loci_including_snps += 1
            missing_loci += 1

    # assign ST
    bst = ','.join(best_st)

    if mismatch_loci_including_snps <= max_missing:
        # only report ST if enough loci are precise matches
        if bst in alleles_to_st:
            # note may have mismatching alleles due to SNPs, this will be recorded in
            # mismatch_loci_including_snps
            bst = alleles_to_st[bst]
        else:
            # determine closest ST
            bst, _, mismatch_loci_including_snps = \
                get_closest_locus_variant(best_st, best_st_annotated, alleles_to_st)
    else:
        bst = '0'

    exact_matches = len(best_st) - mismatch_loci_including_snps
    if exact_matches < required_exact_matches:
        bst = '0'

    # pull info column
    info_final = ''
    if info_arg == 'yes' and bst in st_to_info:
        info_final = st_to_info[bst]
        if report_incomplete and missing_loci > 0:
            info_final += ' (incomplete)'

    if mismatch_loci_including_snps > 0 and bst != '0':
        bst += '-' + str(mismatch_loci_including_snps) + 'LV'

    if min_gene_count is not None and unknown_group_name is not None:
        if info_final == '':
            if sum(0 if x == '-' else 1 for x in best_st_annotated) >= min_gene_count:
                info_final = unknown_group_name
                bst = '0'
            else:
                info_final = '-'

    return bst, best_st_annotated, info_final


def get_allele_and_locus(hit:BlastHit) -> tuple:
    """
    Parses the allele name and locus name from the hit's gene ID.
    """
    if '__' in hit.gene_id:  # srst2 formatted file
        gene_id_components = hit.gene_id.split('__')
        locus = gene_id_components[1]
        allele = gene_id_components[2]
    else:
        allele = hit.gene_id
        locus = hit.gene_id.split('_')[0]
    return allele, locus


def keep_only_one_hit_per_locus(hits:List[BlastHit]) -> List[BlastHit]:
    hits_per_locus = collections.defaultdict(list)
    for hit in hits:
        _, locus = get_allele_and_locus(hit)
        hits_per_locus[locus].append(hit)
    kept_hits = []
    for locus, locus_hits in hits_per_locus.items():
        locus_hits = sorted(locus_hits, key=lambda x: x.score, reverse=True)  # sort best to worst
        best_hit = locus_hits[0]
        kept_hits.append(best_hit)
    return kept_hits


def get_best_allele_per_locus(hits:List[BlastHit], check_for_truncation:bool) -> dict:
    best_scores = {}   # key = locus, value = BLAST score for best allele encountered so far
    best_alleles = {}  # key = locus, value = best allele (* if imprecise match)

    for hit in hits:
        allele, locus = get_allele_and_locus(hit)
        if hit.pcid < 100.00 or hit.alignment_length < hit.ref_length:
            allele += '*'  # inexact match
        if check_for_truncation:
            allele += truncation_check(hit)[0]
        # store best match for each one locus
        if locus in best_scores:
            if hit.score > best_scores[locus]:    # update
                best_scores[locus] = hit.score
                best_alleles[locus] = allele.split('_')[1]  # store number only
        else:  # initialise
            best_scores[locus] = hit.score
            best_alleles[locus] = allele.split('_')[1]  # store number only
    return best_alleles


def process_spurious_hits(hits:List[BlastHit]) -> list:
    hit_strings = []
    for hit in hits:
        allele, locus = get_allele_and_locus(hit)
        if hit.pcid < 100.00 or hit.alignment_length < hit.ref_length:
            allele += '*'  # inexact match
        allele += truncation_check(hit)[0]
        hit_strings.append(allele)
    return hit_strings


def load_st_database(database:str, info_arg:str) -> tuple:
    st_names = []
    alleles_to_st = {}  # key = concatenated string of alleles, value = st
    st_to_info = {}  # key = st, value = info relating to this ST, eg clonal group
    header = []
    with open(database, 'r') as f:
        for line in f:
            fields = line.rstrip().split('\t')
            if len(header) == 0:
                header = fields
                header.pop(0)  # remove st label
                if info_arg == 'yes':
                    header.pop()  # remove info label
            else:
                st = fields.pop(0)
                if info_arg == 'yes':
                    info = fields.pop()
                else:
                    info = ''
                st_names.append(st)
                alleles_to_st[','.join(fields)] = st
                if info_arg == 'yes':
                    st_to_info[st] = info
    return st_names, alleles_to_st, st_to_info, header


def get_closest_locus_variant(query:List[str], annotated_query:List[str], sts:dict) -> tuple:
    annotated_query = list(annotated_query)  # copy the list so we don't change the original
    closest = []
    closest_alleles = {}   # key = st, value = list
    min_dist = len(query)  # number mismatching loci, ignoring SNPs

    for index, item in enumerate(query):
        if item == '-':
            query[index] = '0'

    # get distance from closest ST, ignoring SNPs
    for st in sts:
        d = sum(map(lambda x, y: bool(int(x)-int(y)), st.split(','), query))
        if d == min_dist:
            closest.append(int(sts[st]))
            closest_alleles[sts[st]] = st
        elif d < min_dist:
            # reset
            closest = [int(sts[st])]
            closest_alleles[sts[st]] = st
            min_dist = d  # distance from closest ST, ignoring SNPs

    closest_st = str(min(closest))

    for index, item in enumerate(annotated_query):
        annotated_query[index] = re.sub(r'-\d+%', '', item)
        if item == '-' or '*' in item:
            annotated_query[index] = '0'

    # get distance from closest ST, including SNPs
    min_dist_incl_snps = sum(map(lambda x, y: bool(int(x)-int(y)),
                                 closest_alleles[closest_st].split(','), annotated_query))

    return closest_st, min_dist, min_dist_incl_snps


def add_to_string(existing_str:str, new_str:str) -> str:
    if existing_str == '':
        return new_str
    elif new_str == '':
        return existing_str
    else:
        return existing_str + ',' + new_str


def add_to_strings(existing_strs:str, new_strs:str) -> List[str]:
    assert len(existing_strs) == len(new_strs)
    return [add_to_string(existing_str, new_str)
            for existing_str, new_str in zip(existing_strs, new_strs)]


def cluster_hits_by_contig(hits:List[BlastHit]):
    hits_by_contig = collections.defaultdict(list)
    for h in hits:
        hits_by_contig[h.contig_name].append(h)
    return [x for x in hits_by_contig.values()]
