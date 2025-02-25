"""
I/O utility functions for AiiDA-TEROS.

This module provides functions for file I/O operations, including reading
configuration files and handling VASP files.
"""

import os
import yaml
import jsonschema
import shutil
from aiida.orm import FolderData
from aiida_teros.schemas.config_schema import CONFIG_SCHEMA

def load_config(config_file):
    """
    Load configuration from a YAML file and validate against schema.
    
    Parameters:
        config_file (str): Path to the YAML configuration file.
        
    Returns:
        dict: Configuration dictionary.
        
    Raises:
        FileNotFoundError: If the configuration file does not exist.
        jsonschema.ValidationError: If the configuration is invalid.
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Validate configuration against schema
    jsonschema.validate(instance=config, schema=CONFIG_SCHEMA)
    
    return config

def ensure_directory(directory_path):
    """
    Ensure that a directory exists; create it if it doesn't.
    
    Parameters:
        directory_path (str): Path to the directory.
        
    Returns:
        str: Path to the directory.
    """
    os.makedirs(directory_path, exist_ok=True)
    return directory_path

def save_latex_table(data, headers, filename, title=""):
    """
    Save data as a LaTeX table.
    
    Parameters:
        data (list): List of rows for the table.
        headers (list): List of column headers.
        filename (str): Path to save the LaTeX table.
        title (str): Title for the table section.
    """
    from tabulate import tabulate
    
    # Generate the LaTeX table
    table_latex = tabulate(data, headers=headers, tablefmt="latex")
    
    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    # Save the LaTeX table to a file
    with open(filename, 'w') as f:
        f.write(r"\documentclass{article}\n")
        f.write(r"\usepackage{amsmath}\n")
        f.write(r"\usepackage{geometry}\n")
        f.write(r"\geometry{a4paper, margin=1in}\n")
        f.write(r"\begin{document}\n")
        if title:
            f.write(f"\\section*{{{title}}}\n")
        f.write(table_latex)
        f.write(r"\end{document}\n")

def extract_file_from_folderdata(folderdata, filename, destination=None):
    """
    Extract a file from a FolderData node and optionally save it to a destination.
    
    Parameters:
        folderdata (FolderData): AiiDA FolderData node.
        filename (str): Name of the file to extract.
        destination (str, optional): Destination path to save the file.
                                   If None, the file is saved in the current directory.
                                   
    Returns:
        str: Path to the extracted file.
    """
    if destination is None:
        destination = filename
    
    with folderdata.open(filename, mode='rb') as source:
        with open(destination, mode='wb') as dest:
            shutil.copyfileobj(source, dest)
    
    return destination

def clean_extracted_files(filenames):
    """
    Clean up extracted files.
    
    Parameters:
        filenames (list): List of filenames to remove.
    """
    for filename in filenames:
        if os.path.exists(filename):
            os.remove(filename)