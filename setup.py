#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="aiida-teros",
    version="0.1.0",
    description="AiiDA plugin for Thermodynamics of Oxide Surfaces",
    author="Dr. Thiago Trevizam Dorini",
    author_email="thiagotd@unicamp.br",
    url="https://github.com/DoriniTT/aiida-teros",
    classifiers=[
        "Programming Language :: Python",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "aiida-core>=1.0.0",
        "aiida-vasp>=1.0.0",
        "pymatgen>=2022.0.0",
        "ase>=3.20.0",
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "pint>=0.17",
        "pyyaml>=5.1",
        "tabulate>=0.8.7",
        "click>=7.0",
        "jsonschema>=3.2.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
            "pre-commit>=2.7.0",
            "black>=20.8b1",
            "flake8>=3.8.0",
        ],
        "docs": [
            "sphinx>=3.0.0",
            "sphinx-rtd-theme>=0.5.0",
            "myst-parser>=0.13.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "aiida-teros = aiida_teros.cli:cli",
        ],
    },
    include_package_data=True,
    packages=find_packages(),
    zip_safe=False,
)