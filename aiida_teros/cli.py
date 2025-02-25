#!/usr/bin/env python
"""
Command-line interface for AiiDA-TEROS.

This module provides a command-line interface for running AiiDA-TEROS workflows.
"""

import os
import sys
import time
import click
import yaml
import jsonschema
from ase.io import read
from aiida import load_profile
from aiida.engine import submit
from aiida.orm import load_code, StructureData, Dict, Float, Str, List, Int
from aiida.common.exceptions import NotExistent

from aiida_teros.workflows.teros import TerosWorkChain
from aiida_teros.schemas.config_schema import CONFIG_SCHEMA
from aiida_teros.utils.io import load_config

def ensure_daemon_running():
    """
    Ensure the AiiDA daemon is running.
    
    Returns:
        bool: True if daemon is running, False otherwise.
    """
    from aiida.engine import get_daemon_client
    
    client = get_daemon_client()
    
    if client.is_daemon_running:
        click.echo("AiiDA daemon is running.")
        return True
    else:
        click.echo("AiiDA daemon is not running. Please start it with 'verdi daemon start'.")
        return False

def validate_file_exists(ctx, param, value):
    """
    Validate that a file exists.
    
    Args:
        ctx: Click context.
        param: Click parameter.
        value: Parameter value.
        
    Returns:
        str: Value if file exists.
        
    Raises:
        click.BadParameter: If file doesn't exist.
    """
    if not os.path.exists(value):
        raise click.BadParameter(f"File {value} does not exist.")
    return value

@click.group()
def cli():
    """AiiDA-TEROS command line interface."""
    pass

@cli.command()
@click.argument('config_file', type=click.Path(exists=True), callback=validate_file_exists)
@click.option('--profile', '-p', help="AiiDA profile to use.")
@click.option('--daemon-timeout', type=int, default=30, 
              help="Timeout in seconds for daemon check.")
def run(config_file, profile, daemon_timeout):
    """
    Run the TEROS workflow using the provided configuration file.
    """
    # Load AiiDA profile
    if profile:
        load_profile(profile)
    else:
        try:
            load_profile()
        except Exception as e:
            click.echo(f"Error loading AiiDA profile: {e}")
            sys.exit(1)
    
    # Check if daemon is running
    if not ensure_daemon_running():
        sys.exit(1)
    
    try:
        # Load and validate configuration
        click.echo(f"Loading configuration from {config_file}...")
        config = load_config(config_file)
        
        # Load bulk structure
        click.echo(f"Loading bulk structure from {config['bulk_structure_path']}...")
        bulk_structure = StructureData(ase=read(config['bulk_structure_path']))
        
        # Load bulk metal structure
        click.echo(f"Loading bulk metal structure from {config['bulk_metal_structure_path']}...")
        bulk_metal = StructureData(ase=read(config['bulk_metal_structure_path']))
        
        # Load terminations if specified
        terminations = None
        if 'terminations' in config:
            click.echo("Loading user-specified terminations...")
            terminations = {
                name: StructureData(ase=read(path))
                for name, path in config['terminations'].items()
            }
        
        # Load VASP code
        click.echo(f"Loading VASP code {config['code_label']}...")
        try:
            code = load_code(config['code_label'])
        except NotExistent:
            click.echo(f"Error: Code {config['code_label']} not found in AiiDA database.")
            click.echo("Please set up the code with 'verdi code setup'.")
            sys.exit(1)
        
        # Prepare inputs for the workflow
        inputs = {
            'code': code,
            'bulk_structure': bulk_structure,
            'bulk_metal': bulk_metal,
            'incar_parameters_bulk': config['incar_parameters_bulk'],
            'incar_parameters_bulk_metal': config['incar_parameters_bulk_metal'],
            'incar_parameters_o2': config['incar_parameters_o2'],
            'incar_parameters_slab': config['incar_parameters_slab'],
            'kpoints_precision': Float(config['workflow_settings']['kpoints_precision']),
            'potential_mapping': Dict(dict=config['potential_mapping']),
            'potential_family': Str(config['potential_family']),
            'parser_settings': Dict(dict=config['parser_settings']),
            'computer_options': Dict(dict=config['computer_options']),
            'miller_indices': List(list=config['slab_parameters']['miller_indices']),
            'min_slab_thickness': Float(config['slab_parameters']['min_slab_thickness']),
            'HF_bulk': Float(config['thermodynamic_parameters']['hf_bulk']),
            'precision_phase_diagram': Int(500),  # Default value
        }
        
        # Add vacuum if specified
        if 'vacuum' in config['slab_parameters']:
            inputs['vacuum'] = Float(config['slab_parameters']['vacuum'])
        
        # Add terminations if specified
        if terminations:
            inputs['terminations'] = terminations
        
        # Add path to graphs if specified
        if 'path_to_graphs' in config:
            inputs['path_to_graphs'] = Str(config['path_to_graphs'])
        else:
            inputs['path_to_graphs'] = Str(os.getcwd())
            
        # Add backward compatibility for total energies
        if 'total_energies' in config:
            click.echo("Note: 'total_energies' section is deprecated and will be ignored. O2 and metal energies are now calculated automatically.")
            inputs['total_energy_first_element'] = Float(config['total_energies'].get('total_energy_first_element', 0.0))
            inputs['total_energy_o2'] = Float(config['total_energies'].get('total_energy_o2', 0.0))
        
        # Submit the workflow
        click.echo("Submitting TEROS workflow...")
        process = submit(TerosWorkChain, **inputs)
        
        # Save the process ID to a file
        with open('pks.txt', 'a') as f:
            f.write(f"{process.pk}\n")
        
        click.echo(f"Workflow submitted with PK: {process.pk}")
        click.echo(f"Process ID saved to pks.txt")
        click.echo("You can monitor the progress with 'verdi process list' or 'verdi process show {process.pk}'")
        
    except Exception as e:
        click.echo(f"Error: {e}")
        sys.exit(1)

@cli.command()
@click.argument('config_file', type=click.Path(exists=False))
def validate(config_file):
    """
    Validate a configuration file without running the workflow.
    """
    try:
        # Load configuration
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Validate against schema
        jsonschema.validate(instance=config, schema=CONFIG_SCHEMA)
        
        click.echo(f"Configuration file {config_file} is valid.")
        
        # Check if all file paths exist
        files_to_check = [
            ('bulk_structure_path', "Bulk structure"),
            ('bulk_metal_structure_path', "Bulk metal structure"),
        ]
        
        # Check if all required INCAR parameters are present
        required_incar_params = [
            'incar_parameters_bulk',
            'incar_parameters_bulk_metal',
            'incar_parameters_o2',
            'incar_parameters_slab'
        ]
        
        for param in required_incar_params:
            if param not in config:
                click.echo(f"Warning: Required INCAR parameter set '{param}' not found in configuration.")
                all_files_exist = False
        
        # Add termination files if specified
        if 'terminations' in config:
            for name, path in config['terminations'].items():
                files_to_check.append((f"terminations.{name}", f"Termination {name}"))
        
        all_files_exist = True
        for path_key, description in files_to_check:
            path = None
            if '.' in path_key:
                # Handle nested keys (e.g., terminations.slab1)
                parent_key, child_key = path_key.split('.')
                path = config[parent_key][child_key]
            else:
                path = config[path_key]
            
            if not os.path.exists(path):
                click.echo(f"Warning: {description} file not found: {path}")
                all_files_exist = False
        
        if all_files_exist:
            click.echo("All referenced files exist.")
        
        # Check if code exists in AiiDA database
        try:
            load_profile()
            code = load_code(config['code_label'])
            click.echo(f"VASP code '{config['code_label']}' exists in the AiiDA database.")
        except Exception as e:
            click.echo(f"Warning: Could not verify VASP code '{config['code_label']}': {e}")
        
    except FileNotFoundError:
        click.echo(f"Error: Configuration file {config_file} not found.")
        sys.exit(1)
    except yaml.YAMLError as e:
        click.echo(f"Error: Invalid YAML in configuration file: {e}")
        sys.exit(1)
    except jsonschema.exceptions.ValidationError as e:
        click.echo(f"Error: Configuration file validation failed: {e.message}")
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    cli()
