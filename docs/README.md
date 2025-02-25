# AiiDA-TEROS Documentation

This directory contains documentation for the AiiDA-TEROS package.

## Structure

- `images/` - Contains diagrams, screenshots and other images used in the documentation
- `source/` - Source files for the documentation
- `examples/` - Example files and tutorials

## Building the Documentation

To build the documentation, you need to have Sphinx installed. You can install it along with the required theme and extensions using:

```bash
pip install sphinx sphinx-rtd-theme myst-parser
```

Then, to build the documentation:

```bash
cd docs
make html
```

The built documentation will be available in the `_build/html` directory.

## Contributing to the Documentation

When contributing to the documentation:
1. Write clear and concise text
2. Include examples where appropriate
3. Use proper formatting (headers, lists, code blocks, etc.)
4. Add references to other sections of the documentation or external resources as needed