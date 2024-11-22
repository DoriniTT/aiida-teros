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

import os
import sys
import time
from ase.io import read
from aiida.engine import submit
from aiida import load_profile
from aiida.orm import (
    load_code,
    StructureData,
    Dict,
    Float,
    Int,
    Str,
    List
)
from examples.ag2moo4.AiiDA_teros import AiiDATEROS

# ================================================
# Configuration Section
# ================================================

# Paths and File Names
BULK_STRUCTURE_PATH = 'bulk_ag3po4.vasp'
POTENTIAL_FAMILY = 'your_potential_family_here'  # Example: 'PBE', 'GGA', etc.
CODE_LABEL = 'your_code@here'  # Label of the configured VASP code in AiiDA
PKS_FILE = 'pks.txt'  # File to store submitted WorkChain PKs

# Thermodynamic Parameters
DIVIDE_TO_GET_MIN_COMPOSITION = 2  # For Ag3PO4
HF_BULK = -10.979718545  # Heat of formation for bulk
TOTAL_ENERGY_FIRST_ELEMENT = -2.8289  # Total energy for the first element (e.g., Ag)
TOTAL_ENERGY_O2 = -9.82  # Total energy for O2 molecule

# INCAR Parameters for Bulk and Slab Relaxations
INCAR_PARAMETERS = {
    'bulk': {
        'ISMEAR': 0,
        'SIGMA': 0.01,
        'ENCUT': 550,
        'NCORE': 2,
        'ISPIN': 1,
        'LREAL': 'Auto',
        'PREC': 'Accurate',
        'NELM': 60,
        'NELMIN': 6,
        'EDIFF': 1e-5,
        'LWAVE': True,
        'LORBIT': 11,
        'IVDW': 12,
    },
    'slab': {
        'ISMEAR': 0,
        'SIGMA': 0.01,
        'ENCUT': 550,
        'NCORE': 2,
        'ISPIN': 1,
        'LREAL': 'Auto',
        'PREC': 'Accurate',
        'NELM': 60,
        'NELMIN': 6,
        'EDIFF': 1e-5,
        'LWAVE': True,
        'LORBIT': 11,
        'IVDW': 12,
    }
}

# Workflow Settings
WORKFLOW_SETTINGS = {
    'force_cutoff': 0.01,        # Force convergence criterion in eV/Å
    'steps': 1000,               # Maximum number of relaxation steps
    'kpoints_precision': 0.3,    # K-points mesh density
    'phase_diagram_precision': 500,  # Precision for phase diagram calculations
    'path_to_graphs': os.getcwd(),    # Directory to save generated plots
}

# Potential Mapping
POTENTIAL_MAPPING = {
    'Ag': 'Ag',
    'O': 'O',
    'P': 'P',
}

# Parser Settings
PARSER_SETTINGS = {
    'add_energies': True,
    'add_trajectory': True,
    'add_forces': True,
    'add_structure': True,
    'add_kpoints': True,
}

# Computer Options
COMPUTER_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 40
    },
    'queue_name': 'par40',
}

# Slab Generation Parameters
SLAB_PARAMETERS = {
    'miller_indices': [1, 1, 0],  # Example: [1, 1, 0]
    'min_slab_thickness': 10.0,   # Minimum slab thickness in Å
    'vacuum': 15.0,                # Vacuum spacing in Å
}

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

    # Create AiiDA data nodes for inputs
    inputs = {
        'code': load_vasp_code(CODE_LABEL),
        'bulk_structure': bulk_structure,
        'incar_parameters_bulk': Dict(dict=INCAR_PARAMETERS['bulk']),
        'incar_parameters_slabs': Dict(dict=INCAR_PARAMETERS['slab']),
        'force_cutoff': Float(WORKFLOW_SETTINGS['force_cutoff']),
        'steps': Int(WORKFLOW_SETTINGS['steps']),
        'kpoints_precision': Float(WORKFLOW_SETTINGS['kpoints_precision']),
        'potential_mapping': Dict(dict=POTENTIAL_MAPPING),
        'potential_family': Str(POTENTIAL_FAMILY),
        'parser_settings': Dict(dict=PARSER_SETTINGS),
        'computer_options': Dict(dict=COMPUTER_OPTIONS),
        'miller_indices': List(list=SLAB_PARAMETERS['miller_indices']),
        'min_slab_thickness': Float(SLAB_PARAMETERS['min_slab_thickness']),
        'vacuum': Float(SLAB_PARAMETERS['vacuum']),
        'divide_to_get_minimal_bulk_composition': Int(DIVIDE_TO_GET_MIN_COMPOSITION),
        'precision_phase_diagram': Int(WORKFLOW_SETTINGS['phase_diagram_precision']),
        'HF_bulk': Float(HF_BULK),
        'total_energy_first_element': Float(TOTAL_ENERGY_FIRST_ELEMENT),
        'total_energy_o2': Float(TOTAL_ENERGY_O2),
        'path_to_graphs': Str(WORKFLOW_SETTINGS['path_to_graphs']),
    }

    # Submit the WorkChain
    try:
        print('Submitting AiiDATEROSWorkChain...')
        future = submit(AiiDATEROSWorkChain, **inputs)
        print(f'Submitted AiiDATEROSWorkChain with PK {future.pk}')

        # Save the PK of the WorkChain for future reference
        with open(PKS_FILE, 'a') as f:
            f.write(f'{future.pk}\n')
        print(f'WorkChain PK {future.pk} saved to {PKS_FILE}.')

    except Exception as e:
        sys.exit(f'Failed to submit AiiDATEROSWorkChain: {e}')

if __name__ == '__main__':
    main()
