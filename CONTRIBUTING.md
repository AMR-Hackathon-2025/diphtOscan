# Contributing to diphtOscan

Thank you for your interest in contributing to diphtOscan! This document provides guidelines for contributing to the project.

## Development Environment Setup

### Prerequisites

- Python 3.8 or higher
- Conda (recommended for dependency management)
- Git

### Setting Up the Development Environment

1. **Clone the repository**
   ```bash
   git clone https://github.com/AMR-Hackathon-2025/diphtOscan.git
   cd diphtOscan
   ```

2. **Create and activate the conda environment**
   ```bash
   conda env create -f environment.yml
   conda activate diphtoscan
   ```

3. **Install the package in development mode**
   ```bash
   pip install -e . --no-deps
   ```

4. **Initialize the database**
   ```bash
   diphtoscan -u
   ```

5. **Verify installation**
   ```bash
   diphtoscan --version
   diphtoscan --help
   ```

## Code Style

### Python Style Guidelines

We follow [PEP 8](https://pep8.org/) with the following specifications:

- **Line length**: Maximum 100 characters
- **Indentation**: 4 spaces (no tabs)
- **Imports**: Group in order - standard library, third-party, local
- **Naming conventions**:
  - Functions and variables: `snake_case`
  - Classes: `PascalCase`
  - Constants: `UPPER_SNAKE_CASE`

### Type Hints

Use type hints for function signatures:

```python
def process_assembly(assembly_path: str, threads: int = 4) -> dict:
    ...
```

### Linting

Before submitting code, check for style issues:

```bash
# Install development tools
pip install flake8 black

# Check style
flake8 diphtoscan/

# Auto-format code
black diphtoscan/
```

## Docstring Format

We use [NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html):

```python
def function_name(param1: str, param2: int, optional_param: float = 1.0) -> dict:
    """
    Brief one-line description of the function.

    Longer description if needed, explaining the function's purpose,
    algorithm, or important details.

    Parameters
    ----------
    param1 : str
        Description of param1.
    param2 : int
        Description of param2.
    optional_param : float, optional
        Description of optional_param (default: 1.0).

    Returns
    -------
    dict
        Description of return value.
        Keys include:
        - 'key1': Description
        - 'key2': Description

    Raises
    ------
    ValueError
        When param1 is empty.
    FileNotFoundError
        When the specified file doesn't exist.

    Examples
    --------
    >>> result = function_name("input", 42)
    >>> print(result['key1'])
    'output'

    Notes
    -----
    Additional notes about implementation details or caveats.

    See Also
    --------
    related_function : Description of relationship.
    """
```

### Class Docstrings

```python
class BlastHit:
    """
    Represents a single BLAST hit result.

    This class parses and stores information from a BLAST tabular output line,
    providing convenient access to hit properties.

    Parameters
    ----------
    line : str
        Tab-separated BLAST output line.

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

    Examples
    --------
    >>> hit = BlastHit(blast_output_line)
    >>> print(f"Gene: {hit.gene_id}, Identity: {hit.pcid}%")
    """
```

## Git Workflow

### Branch Naming

- Feature branches: `feature/description-of-feature`
- Bug fixes: `fix/description-of-bug`
- Documentation: `docs/description-of-change`
- Refactoring: `refactor/description-of-change`

### Commit Messages

Follow the [Conventional Commits](https://www.conventionalcommits.org/) format:

```
<type>(<scope>): <short description>

<optional longer description>

<optional footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, no logic change)
- `refactor`: Code refactoring
- `test`: Adding or modifying tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(mlst): add support for novel allele detection

fix(species): correct Mash distance threshold for weak matches

docs(readme): update installation instructions for conda

refactor(utils): simplify genomic context calculation
```

### Pull Request Process

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Write clean, documented code
   - Add or update tests as needed
   - Update documentation if necessary

3. **Test your changes**
   ```bash
   # Run the tool with test data
   diphtoscan -a test_data/*.fasta --mlst --resistance_virulence -o test_output
   ```

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "feat(scope): description of changes"
   ```

5. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request**
   - Provide a clear description of changes
   - Reference any related issues
   - Ensure all checks pass

### PR Review Checklist

Before submitting, ensure:
- [ ] Code follows PEP 8 style guidelines
- [ ] All functions have docstrings
- [ ] Type hints are added for new functions
- [ ] Tests pass with the changes
- [ ] Documentation is updated if needed
- [ ] CHANGELOG.md is updated for significant changes
- [ ] No debug print statements are left in code

## Testing

diphtOscan includes a comprehensive test suite with over 250 tests covering all modules.

### Test Structure

```
tests/
├── conftest.py              # Shared fixtures for all tests
├── fixtures/                # Test data files
├── unit/                    # Unit tests for individual modules
│   ├── test_assembly_qc.py
│   ├── test_blastn.py
│   ├── test_download_alleles_st.py
│   ├── test_jolytree_generation.py
│   ├── test_misc.py
│   ├── test_mlstBLAST.py
│   ├── test_species.py
│   ├── test_template_iTOL.py
│   ├── test_truncation.py
│   ├── test_updating_database.py
│   └── test_utils.py
└── integration/             # Integration tests
    └── test_cli.py
```

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=diphtoscan --cov-report=term-missing

# Run only unit tests
pytest tests/unit/ -v

# Run only integration tests
pytest tests/integration/ -v

# Run tests for a specific module
pytest tests/unit/test_blastn.py -v

# Run tests matching a pattern
pytest tests/ -k "test_blast" -v

# Run tests in parallel (requires pytest-xdist)
pytest tests/ -n auto
```

### Using Test Fixtures

Shared fixtures are defined in `tests/conftest.py`. Common fixtures include:

| Fixture | Description |
|---------|-------------|
| `sample_fasta_file` | Temporary FASTA file with sample sequences |
| `sample_fasta_content` | Raw FASTA content string |
| `mock_blast_hit` | Mock BlastHit object for testing |
| `sample_st_profiles_file` | Temporary ST profiles file |
| `sample_amrfinder_output` | Sample AMRFinder DataFrame |
| `mock_args` | Mock command-line arguments |
| `tmp_path` | pytest built-in for temporary directories |

Example usage in tests:

```python
def test_my_function(sample_fasta_file, tmp_path):
    """Test my function with a sample FASTA file."""
    result = my_function(sample_fasta_file)
    assert result is not None
```

### Adding New Tests

When adding tests for new features:

1. **Unit tests**: Add to `tests/unit/test_<module>.py`
2. **Integration tests**: Add to `tests/integration/test_cli.py`
3. **Fixtures**: Add shared fixtures to `tests/conftest.py`

Follow this pattern for test functions:

```python
def test_function_name_describes_behavior():
    """Docstring explaining what is being tested."""
    # Arrange
    input_data = prepare_test_data()

    # Act
    result = function_under_test(input_data)

    # Assert
    assert result == expected_value
```

### Manual Testing

For manual testing with the actual tool:

```bash
# Basic test with subcommand format
diphtoscan species -a test_assembly.fasta -o test_output

# Full analysis
diphtoscan all -a test_assembly.fasta -st -t -res_vir -o test_output

# Test JSON output
diphtoscan all -a test_assembly.fasta -st --format json -o test_output

# Verify output files exist
ls -la test_output/
```

### Test Data

When adding new features, include appropriate test data if possible. Test assemblies should be:
- Small (to keep repository size manageable)
- Representative of the feature being tested
- Free of any sensitive or confidential data

## Adding New Features

### Database Extensions

To add new genes to the custom database:

1. Add protein sequences to `data/resistance/Corynebacterium_diphtheriae/AMRProt_Cd`
2. Add family annotations to `data/resistance/Corynebacterium_diphtheriae/fam_Cd.tab`
3. Update the `node_class` dictionary in `updating_database.py` if needed
4. Document new genes in `docs/DATABASE.md`

### New Analysis Modules

When adding new analysis capabilities:

1. Create a new module in `diphtoscan/`
2. Add appropriate docstrings and type hints
3. Integrate with `cli.py`
4. Add command-line arguments as needed
5. Update documentation

## Reporting Issues

### Bug Reports

Include:
- diphtOscan version (`diphtoscan --version`)
- Python version (`python --version`)
- Operating system
- Complete error message and traceback
- Minimal example to reproduce the issue
- Input file details (if relevant, without sensitive data)

### Feature Requests

Include:
- Clear description of the proposed feature
- Use case and motivation
- Example of expected behavior
- Any relevant references (papers, tools)

## Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Focus on the code, not the person
- Help others learn and grow

## Questions?

- Open an issue for questions about the codebase
- See the documentation in `docs/` for usage questions
- Contact the maintainers for other inquiries

Thank you for contributing to diphtOscan!
