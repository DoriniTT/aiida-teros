# Tests for AiiDA-TEROS

This directory contains tests for the AiiDA-TEROS package.

## Running Tests

To run the tests, you need to have pytest installed. You can then run:

```bash
cd /path/to/aiida-teros
pytest
```

Or to run with coverage:

```bash
pytest --cov=aiida_teros
```

## Test Data

The `test_data` directory contains sample input files needed for testing, including:
- Sample configuration files
- Sample structure files (VASP POSCAR format)
- Other test-specific data

## Adding New Tests

When adding new tests:
1. Create a new test file with the prefix `test_` (e.g., `test_mymodule.py`)
2. Add test functions with the prefix `test_` (e.g., `def test_myfunction():`)
3. If needed, add test data to the `test_data` directory
4. Document any specific requirements for running the test