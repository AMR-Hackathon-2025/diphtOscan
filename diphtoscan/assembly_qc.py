#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Assembly quality control metrics for diphtOscan.

This module calculates standard assembly quality metrics from FASTA files,
providing quality validation for input genomes before downstream analyses.
"""

from typing import Dict, List, Tuple


def parse_fasta_sequences(fasta_path: str) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file and return sequences with their headers.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    List[Tuple[str, str]]
        List of (header, sequence) tuples.
    """
    sequences = []
    current_header = None
    current_seq_parts = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq_parts)))
                current_header = line[1:].split()[0]
                current_seq_parts = []
            else:
                current_seq_parts.append(line.upper())

        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq_parts)))

    return sequences


def calculate_n50_l50(contig_lengths: List[int]) -> Tuple[int, int]:
    """
    Calculate N50 and L50 metrics from contig lengths.

    N50 is the sequence length of the shortest contig at 50% of the total
    genome length. L50 is the number of contigs whose summed length is N50.

    Parameters
    ----------
    contig_lengths : List[int]
        List of contig lengths in base pairs.

    Returns
    -------
    Tuple[int, int]
        (N50, L50) tuple.
    """
    if not contig_lengths:
        return 0, 0

    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    half_length = total_length / 2

    cumulative = 0
    for i, length in enumerate(sorted_lengths, 1):
        cumulative += length
        if cumulative >= half_length:
            return length, i

    return sorted_lengths[-1], len(sorted_lengths)


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage of a sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence string.

    Returns
    -------
    float
        GC content as a percentage (0-100).
    """
    if not sequence:
        return 0.0

    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    return (gc_count / total) * 100 if total > 0 else 0.0


def calculate_n_content(sequence: str) -> float:
    """
    Calculate percentage of ambiguous bases (N) in a sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence string.

    Returns
    -------
    float
        N content as a percentage (0-100).
    """
    if not sequence:
        return 0.0

    n_count = sequence.count('N')
    total = len(sequence)
    return (n_count / total) * 100 if total > 0 else 0.0


def calculate_assembly_stats(fasta_path: str) -> Dict[str, any]:
    """
    Calculate comprehensive assembly quality metrics from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the input FASTA file.

    Returns
    -------
    Dict[str, any]
        Dictionary containing:
        - total_length (int): Total assembly length in base pairs
        - num_contigs (int): Number of contigs
        - n50 (int): N50 value in base pairs
        - l50 (int): L50 value (number of contigs)
        - largest_contig (int): Length of the largest contig
        - gc_percent (float): GC content percentage
        - n_percent (float): Percentage of ambiguous bases (N)

    Raises
    ------
    FileNotFoundError
        If the FASTA file does not exist.
    ValueError
        If the FASTA file is empty or malformed.
    """
    sequences = parse_fasta_sequences(fasta_path)

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

    contig_lengths = [len(seq) for _, seq in sequences]
    combined_sequence = ''.join(seq for _, seq in sequences)

    n50, l50 = calculate_n50_l50(contig_lengths)

    return {
        'total_length': sum(contig_lengths),
        'num_contigs': len(contig_lengths),
        'n50': n50,
        'l50': l50,
        'largest_contig': max(contig_lengths),
        'gc_percent': round(calculate_gc_content(combined_sequence), 2),
        'n_percent': round(calculate_n_content(combined_sequence), 2)
    }


def get_qc_header() -> List[str]:
    """
    Return the list of QC metric column names.

    Returns
    -------
    List[str]
        List of column names for QC metrics.
    """
    return [
        'total_length',
        'num_contigs',
        'n50',
        'l50',
        'largest_contig',
        'gc_percent',
        'n_percent'
    ]


def format_qc_results(stats: Dict[str, any]) -> Dict[str, str]:
    """
    Format QC statistics for output, prefixing keys with 'qc_'.

    Parameters
    ----------
    stats : Dict[str, any]
        Raw QC statistics from calculate_assembly_stats.

    Returns
    -------
    Dict[str, str]
        Formatted dictionary with 'qc_' prefixed keys.
    """
    return {f'qc_{k}': str(v) for k, v in stats.items()}
