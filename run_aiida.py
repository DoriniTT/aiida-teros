#!/usr/bin/env python
import os
import time
from ase.io import read
from aiida.engine import submit
from aiida import load_profile
from aiida.orm import load_code, StructureData, Dict, Float, Int, Str, List
from workchains.thermo.complete_aiida_thermo.AiiDA_complete_thermo import MasterThermoWorkChain

# Load the AiiDA profile
load_profile()

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

    # Now restart the daemon
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

# Start and restart the daemon
start_and_restart_daemon()

# Path to the bulk structure file
bulk_structure_path = f'/home/thiagotd/git/unicamp_posdoc/calculos/workchains/thermo/complete_aiida_thermo/bulk_ag3po4.vasp'

# Validate that the bulk structure file exists
if not os.path.exists(bulk_structure_path):
    raise FileNotFoundError(f'Bulk structure file not found: {bulk_structure_path}')

# Read the bulk structure and create a StructureData node
bulk_structure = StructureData(ase=read(bulk_structure_path))
divide_to_get_minimal_bulk_composition = Int(2) # Ag3PO4
hf_bulk = Float(-10.979718545)
total_energy_first_element = Float(-2.8289)
total_energy_o2 = Float(-9.82)

# Define INCAR parameters for bulk relaxation
incar_parameters_bulk = {'incar': {
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
}}

# Define INCAR parameters for slab relaxation
incar_parameters_slabs = {'incar': {
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
}}

# Define relaxation settings
force_cutoff = Float(0.01)          # Force convergence criterion
steps = Int(1000)                    # Maximum number of relaxation steps
kpoints_precision = Float(0.3)      # K-points mesh density

# Define potential mapping
potential_mapping = Dict(dict={
    'Ag': 'Ag',
    'O': 'O',
    'P': 'P',
})

# Define potential family
potential_family = Str('PBE')  # Example: 'PBE', 'GGA', etc.

# Define parser settings
parser_settings = Dict(dict={
    'add_energies': True,
    'add_trajectory': True,
    'add_forces': True,
    'add_structure': True,
    'add_kpoints': True,
})

# Define computer options
computer_options = Dict(dict={
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 40
    },
    'queue_name': 'par40',
})

# Define slab generation parameters
miller_indices = List(list=[1, 1, 0])  # Example: [1, 1, 0]
min_slab_thickness = Float(10.0)       # Minimum slab thickness in Angstroms
vacuum = Float(15.0)                    # Vacuum spacing in Angstroms

# Load the VASP code
code_label = 'VASPVTST-6.4.1@bohr-vtst'
code = load_code(code_label)

phase_diagram_precision = Int(500)
path_to_graphs = Str(os.getcwd())

# Prepare inputs
thermo_inputs = {
    'code': code,
    'bulk_structure': bulk_structure,
    'incar_parameters_bulk': incar_parameters_bulk,
    'incar_parameters_slabs': incar_parameters_slabs,
    'force_cutoff': force_cutoff,
    'steps': steps,
    'kpoints_precision': kpoints_precision,
    'potential_mapping': potential_mapping,
    'potential_family': potential_family,
    'parser_settings': parser_settings,
    'computer_options': computer_options,
    'miller_indices': miller_indices,
    'min_slab_thickness': min_slab_thickness,
    'vacuum': vacuum,
    'divide_to_get_minimal_bulk_composition': divide_to_get_minimal_bulk_composition,
    'precision_phase_diagram': phase_diagram_precision,
    'HF_bulk': hf_bulk,
    'total_energy_first_element': total_energy_first_element,
    'total_energy_o2': total_energy_o2,
    'path_to_graphs': path_to_graphs,
}

# Submit the WorkChain
future = submit(MasterThermoWorkChain, **thermo_inputs)

# Save the PK of the WorkChain
with open('pks.txt', 'a') as f:
    f.write(f'{future.pk}\n')

print(f'Submitted MasterThermoWorkChain with PK {future.pk}')