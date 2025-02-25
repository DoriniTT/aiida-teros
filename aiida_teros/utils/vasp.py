"""
VASP utility functions for AiiDA-TEROS.

This module provides functions for configuring VASP calculations.
"""

from aiida.orm import KpointsData, Dict, Str, Bool
from aiida.plugins import WorkflowFactory

VaspWorkflow = WorkflowFactory('vasp.vasp')

def get_vasp_builder(
    code, structure, incar_parameters, potential_family, potential_mapping,
    kpoint_density=0.3, parser_settings=None, computer_options=None,
    label="vasp_calculation", description="", clean_workdir=False, slab=False
):
    """
    Create a builder for a VASP calculation.
    
    Parameters:
        code (Code): VASP code.
        structure (StructureData): Structure to calculate.
        incar_parameters (dict): INCAR parameters.
        potential_family (str): Potential family for VASP.
        potential_mapping (dict): Mapping of elements to potentials.
        kpoint_density (float): K-points mesh density.
        parser_settings (dict): Settings for the parser.
        computer_options (dict): Computational resources and options.
        label (str): Label for the calculation.
        description (str): Description for the calculation.
        clean_workdir (bool): Whether to clean the working directory after calculation.
        slab (bool): Whether the structure is a slab (for k-point adjustment).
        
    Returns:
        ProcessBuilder: Builder for the VASP calculation.
    """
    builder = VaspWorkflow.get_builder()
    
    # Set structure
    builder.structure = structure
    
    # Set INCAR parameters
    builder.parameters = Dict(dict=incar_parameters)
    
    # Set k-points
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(kpoint_density)
    
    # Adjust the k-point mesh for the slab: set the third direction to 1 (for slab calculations)
    if slab:
        kpoints_list = kpoints.get_kpoints_mesh()[0]
        kpoints_list[2] = 1
        kpoints.set_kpoints_mesh(kpoints_list)
    
    builder.kpoints = kpoints
    
    # Set code
    builder.code = code
    
    # Set computational options
    if computer_options:
        builder.options = Dict(dict=computer_options)
    
    # Set metadata
    builder.metadata.label = label
    builder.metadata.description = description
    
    # Set settings (parser settings)
    if parser_settings:
        builder.settings = Dict(dict=parser_settings)
    
    builder.clean_workdir = Bool(clean_workdir)
    
    # Set potential family and mapping
    builder.potential_family = Str(potential_family)
    builder.potential_mapping = Dict(dict=potential_mapping)
    
    return builder