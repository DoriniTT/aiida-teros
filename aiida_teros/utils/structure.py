"""
Structure utility functions for AiiDA-TEROS.

This module provides functions for manipulating structures, generating slabs,
and analyzing structural properties.
"""

import numpy as np
from collections import Counter
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from ase.build import molecule
from aiida.orm import StructureData

def generate_slabs(bulk_structure, miller_indices, min_slab_thickness, vacuum_thickness, 
                  center_slab=True, lll_reduce=True):
    """
    Generate slab structures from a bulk structure.
    
    Parameters:
        bulk_structure (StructureData): AiiDA StructureData of the bulk.
        miller_indices (tuple): Miller indices for the surface orientation (e.g., (1, 1, 0)).
        min_slab_thickness (float): Minimum slab thickness in Angstroms.
        vacuum_thickness (float): Vacuum thickness in Angstroms.
        center_slab (bool): Center the slab in the simulation cell.
        lll_reduce (bool): Apply LLL reduction to find smaller cell.
        
    Returns:
        list: List of StructureData objects representing the generated slabs.
    """
    # Convert AiiDA StructureData to pymatgen Structure
    bulk_pymatgen = bulk_structure.get_pymatgen()
    
    # Generate slabs
    slab_generator = SlabGenerator(
        bulk_pymatgen,
        miller_indices,
        min_slab_thickness,
        vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab
    )
    
    slabs = slab_generator.get_slabs(symmetrize=True)
    
    # Convert slabs to AiiDA StructureData
    aiida_slabs = []
    for slab in slabs:
        orthogonal_slab = slab.get_orthogonal_c_slab()
        ase_atoms = AseAtomsAdaptor().get_atoms(orthogonal_slab)
        aiida_slabs.append(StructureData(ase=ase_atoms))
    
    return aiida_slabs

def get_minimal_composition_factor(structure):
    """
    Calculate the greatest common divisor of the element counts
    to determine the minimal composition factor.
    
    Parameters:
        structure (StructureData): AiiDA StructureData.
        
    Returns:
        int: The GCD of element counts.
    """
    ase_structure = structure.get_ase()
    element_counts = Counter(atom.symbol for atom in ase_structure)
    
    # Find the greatest common divisor for the number of atoms of each element
    gcd = np.gcd.reduce(list(element_counts.values()))
    
    return gcd

def is_binary_oxide(structure):
    """
    Check if the structure is a binary oxide.
    
    Parameters:
        structure (StructureData): AiiDA StructureData.
        
    Returns:
        bool: True if the structure is a binary oxide, False if it's a ternary oxide.
    """
    unique_elements = len(set(structure.get_ase().get_chemical_symbols()))
    
    if unique_elements == 2:
        return True
    elif unique_elements == 3:
        return False
    else:
        raise ValueError("Structure must contain two (binary) or three (ternary) elements.")

def get_element_counts(structure):
    """
    Get element counts in the structure.
    
    Parameters:
        structure (StructureData): AiiDA StructureData.
        
    Returns:
        dict: Dictionary with element symbols as keys and counts as values.
    """
    ase_structure = structure.get_ase()
    element_counts = Counter(atom.symbol for atom in ase_structure)
    return dict(element_counts)

def get_surface_area(structure):
    """
    Calculate the surface area of a slab structure.
    
    Parameters:
        structure (StructureData): AiiDA StructureData of a slab.
        
    Returns:
        float: Surface area in Å².
    """
    cell = structure.get_ase().get_cell()
    return np.linalg.norm(np.cross(cell[0], cell[1]))

def create_o2_molecule(box_size=(13.0, 14.0, 15.0)):
    """
    Create an O2 molecule in a box.
    
    Parameters:
        box_size (tuple): Size of the box in Angstroms (x, y, z).
        
    Returns:
        StructureData: AiiDA StructureData with O2 molecule.
    """
    o2_molecule = molecule('O2')
    o2_molecule.set_cell(box_size)
    o2_molecule.center()  # Center the molecule in the cell
    return StructureData(ase=o2_molecule)