"""Shared fixtures for diphtOscan tests."""

import os
import pytest
from pathlib import Path
from unittest.mock import MagicMock
from argparse import Namespace

import pandas as pd


# Path to fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir():
    """Return the path to the fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture(scope="session")
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """>contig_1
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCG
GGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTT
TTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAG
>contig_2
GCTGGCGCTAATCAGGTTCCATGAGGATGAAATCATCGATGGTGTGGTGGCGGCGGCGA
TGAAAGCCCTGAAAGAGGGCTGTAATATCATCGGTAACCTGCTGCCGACCATTAAAGAG
"""


@pytest.fixture
def sample_fasta_file(tmp_path, sample_fasta_content):
    """Create a temporary FASTA file."""
    fasta_file = tmp_path / "sample.fasta"
    fasta_file.write_text(sample_fasta_content)
    return str(fasta_file)


@pytest.fixture(scope="session")
def blast_hit_line_plus_strand():
    """Sample BLAST output line for plus strand hit."""
    # Fields: sacc, pident, slen, length, score, qseq, sstrand, sstart, send, qacc, qstart, qend, qframe
    return "atpA_1\t100.0\t1500\t1500\t2800.0\tATGAAACGCATTAGCACC\tplus\t1\t1500\tcontig_1\t100\t1599\t1"


@pytest.fixture(scope="session")
def blast_hit_line_minus_strand():
    """Sample BLAST output line for minus strand hit."""
    return "dnaE_2\t99.5\t1200\t1200\t2200.0\tGGTGCTAATGCGTTTCAT\tminus\t1200\t1\tcontig_2\t500\t1699\t-1"


@pytest.fixture
def mock_blast_hit(blast_hit_line_plus_strand):
    """Create a mock BlastHit object."""
    from diphtoscan.blastn import BlastHit
    return BlastHit(blast_hit_line_plus_strand)


@pytest.fixture
def mock_blast_hit_minus(blast_hit_line_minus_strand):
    """Create a mock BlastHit object for minus strand."""
    from diphtoscan.blastn import BlastHit
    return BlastHit(blast_hit_line_minus_strand)


@pytest.fixture(scope="session")
def sample_st_profiles_content():
    """Sample ST profiles content for MLST testing."""
    return """ST\tatpA\tdnaE\tdnaK\tfusA\tleuA\todhA\trpoB
1\t1\t1\t1\t1\t1\t1\t1
2\t1\t2\t1\t1\t1\t1\t1
3\t2\t2\t2\t2\t2\t2\t2
5\t1\t1\t1\t1\t2\t1\t1
"""


@pytest.fixture
def sample_st_profiles_file(tmp_path, sample_st_profiles_content):
    """Create a temporary ST profiles file."""
    profiles_file = tmp_path / "st_profiles.txt"
    profiles_file.write_text(sample_st_profiles_content)
    return str(profiles_file)


@pytest.fixture(scope="session")
def sample_st_profiles_with_info_content():
    """Sample ST profiles content with info column."""
    return """ST\tatpA\tdnaE\tdnaK\tfusA\tleuA\todhA\trpoB\tclonal_complex
1\t1\t1\t1\t1\t1\t1\t1\tCC1
2\t1\t2\t1\t1\t1\t1\t1\tCC1
3\t2\t2\t2\t2\t2\t2\t2\tCC3
"""


@pytest.fixture
def sample_st_profiles_with_info_file(tmp_path, sample_st_profiles_with_info_content):
    """Create a temporary ST profiles file with info column."""
    profiles_file = tmp_path / "st_profiles_info.txt"
    profiles_file.write_text(sample_st_profiles_with_info_content)
    return str(profiles_file)


@pytest.fixture
def sample_amrfinder_output():
    """Sample AMRFinder output DataFrame."""
    return pd.DataFrame({
        'Name': ['strain1', 'strain1', 'strain2'],
        'Protein identifier': ['prot1', 'prot2', 'prot3'],
        'Contig id': ['contig_1', 'contig_1', 'contig_2'],
        'Start': ['100', '500', '200'],
        'Stop': ['400', '800', '600'],
        'Strand': ['+', '+', '-'],
        'Element symbol': ['tetA', 'ermC', 'tetA'],
        'Gene symbol': ['tetA', 'ermC', 'tetA'],
        'Element name': ['tetA gene', 'ermC gene', 'tetA gene'],
        'Scope': ['core', 'core', 'core'],
        'Element type': ['AMR', 'AMR', 'AMR'],
        'Element subtype': ['AMR', 'AMR', 'AMR'],
        'Class': ['TETRACYCLINE', 'MACROLIDE', 'TETRACYCLINE'],
        'Subclass': ['TETRACYCLINE', 'MACROLIDE', 'TETRACYCLINE'],
        'Method': ['EXACTX', 'BLASTX', 'PARTIALX'],
        'Target length': ['300', '300', '400'],
        'Reference sequence length': ['100', '100', '133'],
        '% Coverage of reference sequence': ['100.00', '95.00', '75.00'],
        '% Coverage of reference': ['100.00', '95.00', '75.00'],
        '% Identity to reference sequence': ['100.00', '99.00', '98.00'],
        'Alignment length': ['300', '285', '300'],
        'Accession of closest sequence': ['NG_1234', 'NG_5678', 'NG_1234'],
        'Name of closest sequence': ['tetA', 'ermC', 'tetA'],
        'HMM id': ['HMM001', 'HMM002', 'HMM001'],
        'HMM description': ['desc1', 'desc2', 'desc1'],
        'File': ['/path/to/file1.fasta', '/path/to/file1.fasta', '/path/to/file2.fasta'],
    })


@pytest.fixture
def mock_args(tmp_path):
    """Create mock command line arguments."""
    return Namespace(
        assemblies=[str(tmp_path / "sample.fasta")],
        outdir=str(tmp_path / "output"),
        min_identity=80.0,
        min_coverage=50.0,
        threads=4,
        mlst=True,
        tox=True,
        resistance_virulence=True,
        extend_genotyping=False,
        integron=False,
        tree=False,
        update=False,
        overwrite=False,
        path=str(tmp_path),
        extract=False,
    )


@pytest.fixture(scope="session")
def mock_mash_output():
    """Mock Mash distance output lines."""
    return [
        "Corynebacterium/C.diphtheriae_ref.fasta\tsample.fasta\t0.02\t1.5e-50\t1000/1000\n",
        "Corynebacterium/C.ulcerans_ref.fasta\tsample.fasta\t0.08\t1.0e-30\t800/1000\n",
        "Corynebacterium/C.pseudotub_ref.fasta\tsample.fasta\t0.15\t1.0e-10\t500/1000\n",
    ]


@pytest.fixture(scope="session")
def mock_blast_output_lines():
    """Mock BLAST output lines for testing."""
    return [
        "atpA_1\t100.0\t1500\t1500\t2800.0\tATGAAACGCATTAGCACC\tplus\t1\t1500\tcontig_1\t100\t1599\t1",
        "dnaE_1\t100.0\t1200\t1200\t2200.0\tGCTGGCGCTAATCAGGTT\tplus\t1\t1200\tcontig_1\t1700\t2899\t1",
        "dnaK_1\t99.5\t1800\t1800\t3300.0\tGAAGTTCGGCGGTACATCA\tplus\t1\t1800\tcontig_2\t100\t1899\t1",
    ]


@pytest.fixture(scope="session")
def mock_srst2_formatted_line():
    """Mock BLAST output line with srst2 format gene ID."""
    return "CD__atpA__atpA_1__1\t100.0\t1500\t1500\t2800.0\tATGAAACGCATTAGCACC\tplus\t1\t1500\tcontig_1\t100\t1599\t1"


@pytest.fixture(scope="session")
def sample_mlst_seqs_content():
    """Sample MLST allele sequences content."""
    return """>atpA_1
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCG
GGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTT
>atpA_2
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCG
GGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTA
>dnaE_1
GCTGGCGCTAATCAGGTTCCATGAGGATGAAATCATCGATGGTGTGGTGGCGGCGGCGA
TGAAAGCCCTGAAAGAGGGCTGTAATATCATCGGTAACCTGCTGCCGACCATTAAAGAG
"""


@pytest.fixture
def sample_mlst_seqs_file(tmp_path, sample_mlst_seqs_content):
    """Create a temporary MLST sequences file."""
    seqs_file = tmp_path / "mlst_seqs.fas"
    seqs_file.write_text(sample_mlst_seqs_content)
    return str(seqs_file)


@pytest.fixture
def sample_contig_file(tmp_path):
    """Create a sample contig file for testing is_contig_edge."""
    content = """>contig_1
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCG
GGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTT
TTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAG
"""
    contig_file = tmp_path / "contigs.fasta"
    contig_file.write_text(content)
    return str(contig_file)


@pytest.fixture
def mock_itol_results_df():
    """Mock results DataFrame for iTOL template testing."""
    return pd.DataFrame({
        'spuA': ['spuA', '-', 'spuA'],
        'narG': ['narG', 'narG', '-'],
        'TOXIN': ['tox', 'tox-50%', '-'],
        'AMINOGLYCOSIDE': ['aac(6)', '-', 'aph(3)'],
        'MACROLIDE': ['-', 'ermC', '-'],
    }, index=['strain1', 'strain2', 'strain3'])


@pytest.fixture
def clean_temp_dir(tmp_path):
    """Provide a clean temporary directory and ensure cleanup."""
    yield tmp_path
    # Cleanup is automatic with tmp_path fixture


class MockResponse:
    """Mock requests.Response object."""

    def __init__(self, status_code=200, json_data=None, text=""):
        self.status_code = status_code
        self._json_data = json_data
        self.text = text

    def json(self):
        return self._json_data

    def iter_lines(self):
        for line in self.text.split('\n'):
            yield line.encode('utf-8')

    def iter_content(self, chunk_size=8192):
        yield self.text.encode('utf-8')


@pytest.fixture
def mock_response_200():
    """Create a mock successful response."""
    return MockResponse(status_code=200, json_data={'loci': []})


@pytest.fixture
def mock_response_404():
    """Create a mock 404 response."""
    return MockResponse(status_code=404)


@pytest.fixture
def mock_bigsdb_loci_response():
    """Mock BIGSdb API response for loci."""
    return MockResponse(
        status_code=200,
        json_data={
            'loci': [
                'https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/loci/atpA',
                'https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/loci/dnaE',
            ]
        }
    )


@pytest.fixture
def mock_bigsdb_locus_response():
    """Mock BIGSdb API response for a single locus."""
    return MockResponse(
        status_code=200,
        json_data={
            'id': 'atpA',
            'alleles_fasta': 'https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/loci/atpA/alleles_fasta'
        }
    )
