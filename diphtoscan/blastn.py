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

import os
import subprocess

from typing import List
from .misc import reverse_complement


class BlastHit(object):
    """
    Represents a single BLAST nucleotide alignment hit.

    Parses a tab-separated BLAST output line (outfmt 6) and provides
    structured access to alignment properties.

    Parameters
    ----------
    line : str
        Tab-separated BLAST output line with fields:
        sacc, pident, slen, length, bitscore, qseq, sstrand, sstart,
        send, qacc, qstart, qend, qframe

    Attributes
    ----------
    gene_id : str
        Subject accession (reference gene identifier).
    pcid : float
        Percent identity of the alignment.
    ref_length : int
        Length of the reference sequence.
    alignment_length : int
        Length of the alignment.
    score : float
        BLAST bit score.
    hit_seq : str
        Query sequence from the alignment.
    strand : str
        Strand orientation ('plus' or 'minus').
    ref_start : int
        Start position on reference.
    ref_end : int
        End position on reference.
    contig_name : str
        Query contig accession.
    contig_start : int
        Start position on query contig.
    contig_end : int
        End position on query contig.
    frame : str
        Reading frame on query.
    ref_cov : float
        Fraction of reference covered by alignment (0-1).
    """
    def __init__(self, line):
        fields = line.rstrip().split('\t')
        self.gene_id = fields[0]                # sacc
        self.pcid = float(fields[1])            # pident
        self.ref_length = int(fields[2])        # slen
        self.alignment_length = int(fields[3])  # length
        self.score = float(fields[4])           # score
        self.hit_seq = fields[5]                # qseq
        self.strand = fields[6]                 # sstrand
        self.ref_start = int(fields[7])         # sstart
        self.ref_end = int(fields[8])           # send
        self.contig_name = fields[9]            # qacc
        self.contig_start = int(fields[10])     # qstart
        self.contig_end = int(fields[11])       # qend
        self.frame = fields[12]                 # qframe

        if self.strand == 'plus':
            assert self.ref_end >= self.ref_start
            self.ref_hit_len = self.ref_end - self.ref_start + 1
        else:
            assert self.strand == 'minus'
            assert self.ref_start >= self.ref_end
            self.ref_hit_len = self.ref_start - self.ref_end + 1
        assert self.contig_end >= self.contig_start

        self.ref_cov = self.ref_hit_len / self.ref_length

    def get_seq_start_end_pos_strand(self) -> tuple:
        """
        Extract the nucleotide sequence with position and strand normalization.

        Removes gaps from the alignment sequence and handles strand
        orientation to return sequence in reference orientation.

        Returns
        -------
        tuple
            (nucleotide_sequence, start_position, end_position)
            where positions are oriented relative to the reference.
        """
        # BLAST gives the aligned sequence, so we might need to remove dashes if there are
        # deletions relative to the reference.
        nucl_seq = self.hit_seq.replace('-', '')

        # BLAST also returns the contig's sequence so we might need to flip to the ref strand.
        if self.strand == 'minus':
            return reverse_complement(nucl_seq), self.ref_end, self.ref_start
        else:
            return nucl_seq, self.ref_start, self.ref_end


def run_blastn(db: str, query: str, min_cov: float, min_ident: float) -> List[BlastHit]:
    """
    Execute a BLASTN search and return filtered hits.

    Parameters
    ----------
    db : str
        Path to the FASTA file to use as BLAST database.
        Database files (.nin, etc.) are created if needed.
    query : str
        Path to the query assembly FASTA file.
    min_cov : float
        Minimum reference coverage threshold (0-100).
    min_ident : float
        Minimum percent identity threshold (0-100).

    Returns
    -------
    List[BlastHit]
        List of BlastHit objects passing coverage and identity filters,
        with redundant overlapping hits removed.

    Notes
    -----
    BLAST parameters used:
    - E-value threshold: 1E-20
    - Word size: 32
    - Dust filter: disabled
    - Max target sequences: 10000
    """
    build_blast_database_if_needed(db)

    cmd = 'blastn -task blastn -db {} -query {}'.format(db, query)
    cmd += " -outfmt '6 sacc pident slen length bitscore qseq sstrand sstart send" \
           " qacc qstart qend qframe'"
    # E-value 1E-20: Stringent threshold to ensure statistically significant hits
    # word_size 32: Optimized for finding highly similar sequences (allele matching)
    cmd += ' -dust no -evalue 1E-20 -word_size 32 -max_target_seqs 10000'
    cmd += ' -perc_identity {}'.format(min_ident)

    # TODO: switch this over to subprocess
    blast_hits = []
    f = os.popen(cmd)
    for line in f:
        blast_hits.append(BlastHit(line))
    f.close()

    # Toss out low identity and low coverage hits.
    if min_ident is not None:
        blast_hits = [h for h in blast_hits if h.pcid * 100 >= min_ident]
    if min_cov is not None:
        blast_hits = [h for h in blast_hits if h.ref_cov * 100 >= min_cov]

    # Clean up redundant hits so we only have one hit for each part of the genome.
    blast_hits = cull_redundant_hits(blast_hits)

    # Get rid of any hits that start with 'delete_'. This is used in cases where there can be
    # confounding homologous hits. E.g. searching for rmpA2 alleles but not wanting to be misled by
    # homologous rmpA alleles.
    blast_hits = [h for h in blast_hits if not h.gene_id.startswith('delete_')]
    return blast_hits


def cull_redundant_hits(blast_hits: list) -> List[BlastHit]:
    """
    Remove overlapping BLAST hits, keeping only the best for each region.

    Implements custom culling logic similar to BLAST's -culling_limit 1
    but with more sophisticated overlap detection.

    Parameters
    ----------
    blast_hits : list
        List of BlastHit objects to filter.

    Returns
    -------
    List[BlastHit]
        Filtered list with redundant overlapping hits removed.
        Hits are sorted by quality (pcid * score * ref_cov).
    """
    # Sort the hits from best to worst. Hit quality is defined as the product of gene coverage,
    # identity and score.
    blast_hits = sorted(blast_hits, key=lambda x: (1/(x.pcid * x.score * x.ref_cov), x.gene_id))

    filtered_blast_hits = []

    for h in blast_hits:
        if not overlapping(h, filtered_blast_hits):
            filtered_blast_hits.append(h)

    return filtered_blast_hits


def overlapping(hit: BlastHit, existing_hits: list) -> bool:
    """
    Check if a hit overlaps significantly with any existing hits.

    Parameters
    ----------
    hit : BlastHit
        The hit to check for overlaps.
    existing_hits : list
        List of previously accepted BlastHit objects.

    Returns
    -------
    bool
        True if hit overlaps with any existing hit in the same
        contig, strand, and frame.
    """
    # Only consider hits in the same reading frame.
    existing_hits = [h for h in existing_hits if
                     h.strand == hit.strand and h.frame == hit.frame and
                     h.contig_name == hit.contig_name]

    for existing_hit in existing_hits:
        if hits_overlap(hit, existing_hit):
            return True

    return False


def hits_overlap(a: BlastHit, b: BlastHit) -> bool:
    """
    Determine if two BLAST hits have significant positional overlap.

    Parameters
    ----------
    a : BlastHit
        First hit to compare.
    b : BlastHit
        Second hit to compare.

    Returns
    -------
    bool
        True if hits overlap by more than 50 bp on the query contig.
    """
    if a.contig_start <= b.contig_end and b.contig_start <= a.contig_end:  # There is some overlap
        allowed_overlap = 50
        overlap_size = len(range(max(a.contig_start, b.contig_start),
                                 min(a.contig_end, b.contig_end) + 1))
        return overlap_size > allowed_overlap
    else:
        return False


def build_blast_database_if_needed(seqs: str) -> None:
    """
    Create BLAST nucleotide database if it doesn't exist.

    Parameters
    ----------
    seqs : str
        Path to the FASTA file to index.

    Notes
    -----
    Creates .nin, .nsq, .nhr files alongside the FASTA file.
    """
    if not os.path.exists(seqs + '.nin'):
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call('makeblastdb -dbtype nucl -in ' + seqs, stdout=devnull,
                                  shell=True)
