#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example script for running TEROS workflow for Ag2O binary oxide system.

This script demonstrates how to set up and run the TEROS workflow for
calculating the surface thermodynamics of Ag2O. It can be run in two ways:
1. Using the command-line interface (recommended)
2. Using the Python API directly

Both methods are shown below.
"""

import os
import sys
import subprocess
from pathlib import Path

def run_cli_method():
    """Run the workflow using the AiiDA-TEROS CLI."""
    
    print("=" * 80)
    print("RUNNING TEROS WORKFLOW FOR Ag2O USING CLI METHOD")
    print("=" * 80)
    
    # Get the current script directory
    current_dir = Path(__file__).parent.absolute()
    
    # Path to the configuration file
    config_file = current_dir / "config.yaml"
    
    # Check if the config file exists
    if not config_file.exists():
        print(f"Error: Configuration file not found at {config_file}")
        sys.exit(1)
    
    # Run the TEROS workflow using the CLI
    print(f"Running: aiida-teros run {config_file}")
    try:
        result = subprocess.run(
            ["aiida-teros", "run", str(config_file)],
            check=True,
            text=True,
            capture_output=True
        )
        print(result.stdout)
        if result.stderr:
            print(f"Warnings/Errors: {result.stderr}")
    except subprocess.CalledProcessError as e:
        print(f"Error running aiida-teros CLI: {e}")
        if e.stdout:
            print(f"Output: {e.stdout}")
        if e.stderr:
            print(f"Error: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: aiida-teros command not found. Make sure AiiDA-TEROS is installed properly.")
        sys.exit(1)
    
    print("\nWorkflow submitted successfully using CLI method!")

def run_api_method():
    """Run the workflow using the AiiDA-TEROS Python API."""
    
    print("=" * 80)
    print("RUNNING TEROS WORKFLOW FOR Ag2O USING PYTHON API METHOD")
    print("=" * 80)
    
    import yaml
    from ase.io import read
    from aiida import load_profile
    from aiida.engine import submit
    from aiida.orm import (
        load_code,
        StructureData,
        Dict,
        Float,
        Str,
        List,
        Int
    )
    from aiida_teros.workflows.teros import TerosWorkChain
    from aiida_teros.utils.io import load_config
    
    # Get the current script directory
    current_dir = Path(__file__).parent.absolute()
    
    # Path to the configuration file
    config_file = current_dir / "config.yaml"
    
    # Load AiiDA profile
    load_profile()
    
    # Load configuration
    config = load_config(str(config_file))
    
    # Load bulk structure
    print(f"Loading bulk structure from {config['bulk_structure_path']}...")
    bulk_structure_path = current_dir / config['bulk_structure_path']
    bulk_structure = StructureData(ase=read(bulk_structure_path))
    
    # Load bulk metal structure
    print(f"Loading bulk metal structure from {config['bulk_metal_structure_path']}...")
    bulk_metal_path = current_dir / config['bulk_metal_structure_path']
    bulk_metal = StructureData(ase=read(bulk_metal_path))
    
    # Load VASP code
    print(f"Loading VASP code {config['code_label']}...")
    try:
        code = load_code(config['code_label'])
    except Exception as e:
        print(f"Error: Could not load code '{config['code_label']}': {e}")
        print("Please make sure the code is set up in AiiDA with 'verdi code setup'.")
        sys.exit(1)
    
    # Prepare inputs for the workflow
    inputs = {
        'code': code,
        'bulk_structure': bulk_structure,
        'bulk_metal': bulk_metal,
        'incar_parameters_bulk': config['incar_parameters_bulk'],
        'incar_parameters_bulk_metal': config['incar_parameters_bulk_metal'],
        'incar_parameters_slab': config['incar_parameters_slab'],
        'kpoints_precision': Float(config['workflow_settings']['kpoints_precision']),
        'potential_mapping': Dict(dict=config['potential_mapping']),
        'potential_family': Str(config['potential_family']),
        'parser_settings': Dict(dict=config['parser_settings']),
        'computer_options': Dict(dict=config['computer_options']),
        'miller_indices': List(list=config['slab_parameters']['miller_indices']),
        'min_slab_thickness': Float(config['slab_parameters']['min_slab_thickness']),
        'HF_bulk': Float(config['thermodynamic_parameters']['hf_bulk']),
        'total_energy_first_element': Float(config['total_energies']['total_energy_first_element']),
        'total_energy_o2': Float(config['total_energies']['total_energy_o2']),
        'precision_phase_diagram': Int(500),  # Default value
    }
    
    # Add vacuum if specified
    if 'vacuum' in config['slab_parameters']:
        inputs['vacuum'] = Float(config['slab_parameters']['vacuum'])
    
    # Add path to graphs if specified
    if 'path_to_graphs' in config:
        inputs['path_to_graphs'] = Str(config['path_to_graphs'])
    else:
        inputs['path_to_graphs'] = Str(str(current_dir))
    
    # Submit the workflow
    print("Submitting TEROS workflow...")
    process = submit(TerosWorkChain, **inputs)
    
    # Save the process ID to a file
    with open(current_dir / 'pks.txt', 'a') as f:
        f.write(f"{process.pk}\n")
    
    print(f"Workflow submitted with PK: {process.pk}")
    print(f"Process ID saved to pks.txt")
    print("You can monitor the progress with 'verdi process list' or 'verdi process show {process.pk}'")
    
    print("\nWorkflow submitted successfully using Python API method!")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--api":
        # Run using Python API
        run_api_method()
    else:
        # Run using CLI (default)
        run_cli_method()
    
    print("\nNOTE: To monitor the calculation progress, run:")
    print("    verdi process list")
    print("    verdi process show <PK>")