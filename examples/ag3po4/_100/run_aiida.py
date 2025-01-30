#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
run_aiida.py

This script sets up and submits the AiiDATEROSWorkChain for automating surface thermodynamics
calculations using AiiDA and VASP. It ensures the AiiDA daemon is running, prepares all necessary
inputs, and submits the workflow for execution.

Usage:
    python run_aiida.py
"""

import os, sys, time, yaml
import importlib
from ase.io import read
from aiida.engine import submit
from aiida import load_profile
from aiida.orm import (
    load_code,
    StructureData,
    Dict,
    Float,
    Str,
    List
)
# ================================================
# Configuration Section
# ================================================

# Load configuration from config.yaml
CONFIG_FILE = 'config.yaml'

def load_config(config_file):
    """
    Load configuration from a YAML file.

    :param config_file: Path to the YAML configuration file.
    :return: Configuration dictionary.
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f'Configuration file not found: {config_file}')
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    print(f'Configuration loaded from {config_file}.')
    return config

# Load configuration
config = load_config(CONFIG_FILE)

# Retrieve the module path from the configuration dictionary
module_path = config['module_path']
# Import the module dynamically using the module path
module = importlib.import_module(module_path)
# Get the AiiDATEROSWorkChain class from the imported module
AiiDATEROSWorkChain = getattr(module, "AiiDATEROSWorkChain")

# Extract configuration parameters
BULK_STRUCTURE_PATH = config['bulk_structure_path']
BULK_METAL_STRUCTURE_PATH = config['bulk_metal_structure_path']
POTENTIAL_FAMILY = config['potential_family']
CODE_LABEL = config['code_label']
HF_BULK = config['thermodynamic_parameters']['hf_bulk']
INCAR_PARAMETERS_BULK = config['incar_parameters_bulk']
INCAR_PARAMETERS_BULK_METAL = config['incar_parameters_bulk_metal']
INCAR_PARAMETERS_O2 = config['incar_parameters_o2']
INCAR_PARAMETERS_SLAB = config['incar_parameters_slab']
WORKFLOW_SETTINGS = config['workflow_settings']
POTENTIAL_MAPPING = config['potential_mapping']
PARSER_SETTINGS = config['parser_settings']
COMPUTER_OPTIONS = config['computer_options']
SLAB_PARAMETERS = config['slab_parameters']

# Check if 'terminations' exist in the config.yaml
if 'terminations' in config:
    dict_terminations = config['terminations']
    TERMINATIONS = {struc: StructureData(ase=read(termination)) for struc, termination in dict_terminations.items()} # Load terminations
else:
    TERMINATIONS = None

# ================================================
# Helper Functions
# ================================================

def start_and_restart_daemon(max_retries=5, delay=3):
    """
    Ensure the AiiDA daemon is running by starting and restarting it if necessary.

    :param max_retries: Maximum number of retries for starting/restarting the daemon.
    :param delay: Delay in seconds between retries.
    """
    from aiida.engine import get_daemon_client

    client = get_daemon_client()

    # Check if the daemon is already running
    if client.is_daemon_running:
        print('Daemon is already running.')
    else:
        # Attempt to start the daemon
        print('Daemon is not running. Attempting to start the daemon.')
        client.start_daemon()

        # Wait and check if the daemon started successfully
        for attempt in range(max_retries):
            time.sleep(delay)
            if client.is_daemon_running:
                print('Daemon started successfully.')
                break
            else:
                print(f'Daemon not running after attempt {attempt + 1}. Retrying...')
                client.start_daemon()
        else:
            raise RuntimeError('Failed to start the daemon after multiple attempts.')

    # Restart the daemon
    print('Restarting the daemon.')
    client.restart_daemon()

    # Wait and check if the daemon restarted successfully
    for attempt in range(max_retries):
        time.sleep(delay)
        if client.is_daemon_running:
            print('Daemon restarted successfully.')
            break
        else:
            print(f'Daemon not running after restart attempt {attempt + 1}. Retrying restart...')
            client.restart_daemon()
    else:
        raise RuntimeError('Failed to restart the daemon after multiple attempts.')

def validate_file_exists(file_path):
    """
    Validate that a file exists at the given path.

    :param file_path: Path to the file.
    :raises FileNotFoundError: If the file does not exist.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f'Bulk structure file not found: {file_path}')

def load_bulk_structure(file_path):
    """
    Load the bulk structure from a VASP POSCAR file and create a StructureData node.

    :param file_path: Path to the POSCAR file.
    :return: AiiDA StructureData node.
    """
    try:
        ase_structure = read(file_path)
        structure = StructureData(ase=ase_structure)
        print(f'Loaded bulk structure from {file_path}.')
        return structure
    except Exception as e:
        raise RuntimeError(f'Failed to read bulk structure file: {e}')

def load_vasp_code(code_label):
    """
    Load the VASP code from AiiDA using the given label.

    :param code_label: Label of the VASP code configured in AiiDA.
    :return: AiiDA Code node.
    """
    try:
        code = load_code(code_label)
        print(f'Loaded VASP code: {code_label}.')
        return code
    except Exception as e:
        raise RuntimeError(f'Failed to load VASP code "{code_label}": {e}')

# ================================================
# Main Execution
# ================================================

def main():
    # Load the AiiDA profile
    try:
        load_profile()
        print('AiiDA profile loaded successfully.')
    except Exception as e:
        sys.exit(f'Failed to load AiiDA profile: {e}')

    # Start and restart the AiiDA daemon
    try:
        start_and_restart_daemon()
    except RuntimeError as e:
        sys.exit(e)

    # Validate the bulk structure file
    try:
        validate_file_exists(BULK_STRUCTURE_PATH)
    except FileNotFoundError as e:
        sys.exit(e)

    # Load the bulk structure
    bulk_structure = load_bulk_structure(BULK_STRUCTURE_PATH)
    bulk_metal_structure = load_bulk_structure(BULK_METAL_STRUCTURE_PATH)

    # Create AiiDA data nodes for inputs
    inputs = {
        'code': load_vasp_code(CODE_LABEL),
        'bulk_structure': bulk_structure,
        'bulk_metal': bulk_metal_structure,
        'incar_parameters_bulk': INCAR_PARAMETERS_BULK,
        'incar_parameters_bulk_metal': INCAR_PARAMETERS_BULK_METAL,
        'incar_parameters_o2': INCAR_PARAMETERS_O2,
        'incar_parameters_slab': INCAR_PARAMETERS_SLAB,
        'kpoints_precision': Float(WORKFLOW_SETTINGS['kpoints_precision']),
        'potential_mapping': Dict(dict=POTENTIAL_MAPPING),
        'potential_family': Str(POTENTIAL_FAMILY),
        'parser_settings': Dict(dict=PARSER_SETTINGS),
        'computer_options': Dict(dict=COMPUTER_OPTIONS),
        'miller_indices': List(list=SLAB_PARAMETERS['miller_indices']),
        'min_slab_thickness': Float(SLAB_PARAMETERS['min_slab_thickness']),
        'HF_bulk': Float(HF_BULK),
    }

    # Add terminations if they exist
    if TERMINATIONS is not None:
        inputs['terminations'] = TERMINATIONS

    # Submit the WorkChain
    try:
        print('Submitting AiiDATEROSWorkChain...')
        future = submit(AiiDATEROSWorkChain, **inputs)
        print(f'Submitted AiiDATEROSWorkChain with PK {future.pk}')

        PKS_FILE = 'pks.txt'  # File to store submitted WorkChain PKs
        # Save the PK of the WorkChain for future reference
        with open(PKS_FILE, 'a') as f:
            f.write(f'{future.pk}\n')
        print(f'WorkChain PK {future.pk} saved to {PKS_FILE}.')

    except Exception as e:
        sys.exit(f'Failed to submit AiiDATEROSWorkChain: {e}')

if __name__ == '__main__':
    main()