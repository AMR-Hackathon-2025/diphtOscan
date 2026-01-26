"""Tests for diphtoscan/mlstBLAST.py - MLST typing functions."""

import pytest
from unittest.mock import MagicMock, patch

from diphtoscan.mlstBLAST import (
    get_allele_and_locus,
    keep_only_one_hit_per_locus,
    get_best_allele_per_locus,
    load_st_database,
    get_closest_locus_variant,
    add_to_string,
    add_to_strings,
    cluster_hits_by_contig,
)
from diphtoscan.blastn import BlastHit


class TestGetAlleleAndLocus:
    """Tests for get_allele_and_locus function."""

    def test_srst2_format(self, mock_srst2_formatted_line):
        """Test parsing srst2 formatted gene ID."""
        hit = BlastHit(mock_srst2_formatted_line)
        allele, locus = get_allele_and_locus(hit)

        assert locus == 'atpA'
        assert allele == 'atpA_1'

    def test_simple_format(self, blast_hit_line_plus_strand):
        """Test parsing simple format gene ID (locus_allele)."""
        hit = BlastHit(blast_hit_line_plus_strand)
        allele, locus = get_allele_and_locus(hit)

        assert locus == 'atpA'
        assert allele == 'atpA_1'

    def test_multiple_underscores_simple_format(self):
        """Test parsing gene ID with multiple underscores."""
        line = "rpoB_variant_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1"
        hit = BlastHit(line)
        allele, locus = get_allele_and_locus(hit)

        assert locus == 'rpoB'
        assert allele == 'rpoB_variant_1'


class TestKeepOnlyOneHitPerLocus:
    """Tests for keep_only_one_hit_per_locus function."""

    def test_single_hit_per_locus(self):
        """Test that single hits per locus are all kept."""
        lines = [
            "atpA_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1",
            "dnaE_1\t100.0\t1200\t1200\t2200.0\tATG\tplus\t1\t1200\tcontig_1\t2000\t3199\t1",
        ]
        hits = [BlastHit(line) for line in lines]

        result = keep_only_one_hit_per_locus(hits)

        assert len(result) == 2
        loci = {get_allele_and_locus(h)[1] for h in result}
        assert loci == {'atpA', 'dnaE'}

    def test_multiple_hits_per_locus_keeps_best(self):
        """Test that only the best hit per locus is kept."""
        lines = [
            "atpA_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1",
            "atpA_2\t99.0\t1500\t1500\t2700.0\tATG\tplus\t1\t1500\tcontig_2\t100\t1599\t1",  # Lower score
        ]
        hits = [BlastHit(line) for line in lines]

        result = keep_only_one_hit_per_locus(hits)

        assert len(result) == 1
        assert result[0].gene_id == 'atpA_1'  # Best score

    def test_empty_hits_list(self):
        """Test handling of empty hits list."""
        result = keep_only_one_hit_per_locus([])
        assert result == []


class TestGetBestAllelePerLocus:
    """Tests for get_best_allele_per_locus function."""

    def test_exact_match(self):
        """Test exact match returns allele number without marker."""
        line = "atpA_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1"
        hits = [BlastHit(line)]

        result = get_best_allele_per_locus(hits, check_for_truncation=False)

        assert 'atpA' in result
        assert result['atpA'] == '1'  # No * marker

    def test_inexact_match_marked_with_asterisk(self):
        """Test inexact match (identity < 100%) is marked with asterisk."""
        # pcid is 99.5, not 100
        line = "atpA_1\t99.5\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1"
        hits = [BlastHit(line)]

        result = get_best_allele_per_locus(hits, check_for_truncation=False)

        assert 'atpA' in result
        assert '*' in result['atpA']

    def test_partial_alignment_marked_with_asterisk(self):
        """Test partial alignment (length < ref_length) is marked with asterisk."""
        # alignment_length is 1400, ref_length is 1500
        line = "atpA_1\t100.0\t1500\t1400\t2600.0\tATG\tplus\t1\t1400\tcontig_1\t100\t1499\t1"
        hits = [BlastHit(line)]

        result = get_best_allele_per_locus(hits, check_for_truncation=False)

        assert 'atpA' in result
        assert '*' in result['atpA']

    def test_best_score_selected(self):
        """Test that the hit with the best score is selected."""
        lines = [
            "atpA_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1",
            "atpA_2\t100.0\t1500\t1500\t2900.0\tATG\tplus\t1\t1500\tcontig_2\t100\t1599\t1",  # Higher score
        ]
        hits = [BlastHit(line) for line in lines]

        result = get_best_allele_per_locus(hits, check_for_truncation=False)

        assert 'atpA' in result
        assert result['atpA'] == '2'  # Best score


class TestLoadStDatabase:
    """Tests for load_st_database function."""

    def test_load_without_info(self, sample_st_profiles_file):
        """Test loading ST database without info column."""
        st_names, alleles_to_st, st_to_info, header = load_st_database(
            sample_st_profiles_file, 'no'
        )

        assert len(header) == 7  # atpA, dnaE, dnaK, fusA, leuA, odhA, rpoB
        assert '1' in st_names
        assert '2' in st_names
        assert '1,1,1,1,1,1,1' in alleles_to_st
        assert alleles_to_st['1,1,1,1,1,1,1'] == '1'
        assert st_to_info == {}

    def test_load_with_info(self, sample_st_profiles_with_info_file):
        """Test loading ST database with info column."""
        st_names, alleles_to_st, st_to_info, header = load_st_database(
            sample_st_profiles_with_info_file, 'yes'
        )

        assert len(header) == 7  # Info column removed from header
        assert '1' in st_to_info
        assert st_to_info['1'] == 'CC1'

    def test_alleles_to_st_mapping(self, sample_st_profiles_file):
        """Test that alleles correctly map to ST."""
        st_names, alleles_to_st, _, header = load_st_database(
            sample_st_profiles_file, 'no'
        )

        # ST1 has all 1s
        assert alleles_to_st['1,1,1,1,1,1,1'] == '1'
        # ST2 has dnaE=2
        assert alleles_to_st['1,2,1,1,1,1,1'] == '2'


class TestGetClosestLocusVariant:
    """Tests for get_closest_locus_variant function."""

    def test_exact_match(self):
        """Test exact match returns ST with 0 distance."""
        query = ['1', '1', '1', '1', '1', '1', '1']
        annotated_query = ['1', '1', '1', '1', '1', '1', '1']
        sts = {'1,1,1,1,1,1,1': '1', '2,2,2,2,2,2,2': '2'}

        closest_st, min_dist, min_dist_incl_snps = get_closest_locus_variant(
            query, annotated_query, sts
        )

        assert closest_st == '1'
        assert min_dist == 0
        assert min_dist_incl_snps == 0

    def test_single_locus_variant(self):
        """Test single locus variant detection."""
        query = ['1', '2', '1', '1', '1', '1', '1']  # dnaE differs
        annotated_query = ['1', '2', '1', '1', '1', '1', '1']
        sts = {'1,1,1,1,1,1,1': '1', '2,2,2,2,2,2,2': '2'}

        closest_st, min_dist, min_dist_incl_snps = get_closest_locus_variant(
            query, annotated_query, sts
        )

        assert closest_st == '1'  # Closest to ST1
        assert min_dist == 1  # 1 locus different

    def test_missing_allele_handled(self):
        """Test that missing alleles ('-') are handled."""
        query = ['-', '1', '1', '1', '1', '1', '1']  # atpA missing
        annotated_query = ['-', '1', '1', '1', '1', '1', '1']
        sts = {'1,1,1,1,1,1,1': '1'}

        # Missing alleles are converted to '0'
        closest_st, min_dist, min_dist_incl_snps = get_closest_locus_variant(
            query, annotated_query, sts
        )

        assert closest_st == '1'
        assert min_dist >= 1  # At least one mismatch due to missing allele


class TestAddToString:
    """Tests for add_to_string function."""

    def test_empty_existing(self):
        """Test adding to empty string."""
        assert add_to_string('', 'new') == 'new'

    def test_empty_new(self):
        """Test adding empty string."""
        assert add_to_string('existing', '') == 'existing'

    def test_both_non_empty(self):
        """Test adding two non-empty strings."""
        assert add_to_string('existing', 'new') == 'existing,new'

    def test_both_empty(self):
        """Test adding two empty strings."""
        assert add_to_string('', '') == ''


class TestAddToStrings:
    """Tests for add_to_strings function."""

    def test_combine_lists(self):
        """Test combining two lists of strings."""
        existing = ['a', 'b', 'c']
        new = ['1', '2', '3']

        result = add_to_strings(existing, new)

        assert result == ['a,1', 'b,2', 'c,3']

    def test_with_empty_strings(self):
        """Test combining lists with empty strings."""
        existing = ['a', '', 'c']
        new = ['', '2', '']

        result = add_to_strings(existing, new)

        assert result == ['a', '2', 'c']

    def test_empty_lists(self):
        """Test combining empty lists."""
        result = add_to_strings([], [])
        assert result == []


class TestClusterHitsByContig:
    """Tests for cluster_hits_by_contig function."""

    def test_single_contig(self):
        """Test clustering hits from single contig."""
        lines = [
            "atpA_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1",
            "dnaE_1\t100.0\t1200\t1200\t2200.0\tATG\tplus\t1\t1200\tcontig_1\t2000\t3199\t1",
        ]
        hits = [BlastHit(line) for line in lines]

        result = cluster_hits_by_contig(hits)

        assert len(result) == 1
        assert len(result[0]) == 2

    def test_multiple_contigs(self):
        """Test clustering hits from multiple contigs."""
        lines = [
            "atpA_1\t100.0\t1500\t1500\t2800.0\tATG\tplus\t1\t1500\tcontig_1\t100\t1599\t1",
            "dnaE_1\t100.0\t1200\t1200\t2200.0\tATG\tplus\t1\t1200\tcontig_2\t100\t1299\t1",
        ]
        hits = [BlastHit(line) for line in lines]

        result = cluster_hits_by_contig(hits)

        assert len(result) == 2
        contig_names = [group[0].contig_name for group in result]
        assert set(contig_names) == {'contig_1', 'contig_2'}

    def test_empty_hits(self):
        """Test clustering empty hits list."""
        result = cluster_hits_by_contig([])
        assert result == []
