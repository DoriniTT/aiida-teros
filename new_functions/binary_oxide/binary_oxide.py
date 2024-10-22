import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read
from aiida.orm import StructureData
from aiida import load_profile
from collections import Counter
load_profile()

def surface_total_energy(slab_structure, bulk_structure, E_slab, E_bulk, mu_O, gcd_value):
    """
    Calculate the surface total energy according to the given equation.

    Parameters:
    slab_structure (StructureData): AiiDA StructureData of the slab.
    bulk_structure (StructureData): AiiDA StructureData of the bulk.
    E_slab (float): DFT total energy of the slab.
    E_bulk (float): DFT total energy of the bulk.
    mu_O (float): Chemical potential of oxygen.
    gcd_value (int): Value by which the number of atoms needs to be divided to get the minimal stoichiometry.

    Returns:
    float: Surface Gibbs free energy.
    """
    # Convert StructureData to ASE Atoms
    slab_structure = slab_structure.get_ase()
    bulk_structure = bulk_structure.get_ase()

    element_symbols = list(set(atom.symbol for atom in slab_structure if atom.symbol != 'O'))
    if len(element_symbols) != 1:
        raise ValueError("Structure must contain exactly one non-oxygen element.")
    element = element_symbols[0]

    N_element_slab = len([atom for atom in slab_structure if atom.symbol == element])
    N_O_slab = len([atom for atom in slab_structure if atom.symbol == 'O'])

    if N_element_slab == 0 or N_O_slab == 0:
        raise ValueError("Number of non-oxygen or oxygen atoms cannot be zero.")
    
    lengths = slab_structure.get_cell_lengths_and_angles()
    if lengths[0] <= 0 or lengths[1] <= 0:
        raise ValueError("Cell dimensions must be positive values.")
    
    A = lengths[0] * lengths[1]  # Area is the product of the x and y coordinates of the cell

    if A == 0:
        raise ValueError("Surface area cannot be zero.")
    
    N_element_bulk = len([atom for atom in bulk_structure if atom.symbol == element])
    N_O_bulk = len([atom for atom in bulk_structure if atom.symbol == 'O'])
    
    if N_element_bulk == 0:
        raise ValueError("Bulk structure must contain the element.")
    
    E_bulk_per_fu = E_bulk / gcd_value  # Energy per formula unit
    
    gamma = (E_slab - N_element_slab * E_bulk_per_fu + (2 * N_element_slab - N_O_slab) * mu_O) / (2 * A)
    if not np.isfinite(gamma):
        raise ValueError("Calculated surface Gibbs free energy is not finite. Check input values.")
    
    return gamma

def chemical_potential_limits(delta_Hf, E_total_O2):
    """
    Calculate the limits for the chemical potential of oxygen.

    Parameters:
    delta_Hf (float): Formation enthalpy at 0 K and 0 Pa.
    E_total_O2 (float): Total energy of O2.

    Returns:
    tuple: Lower and upper limits for the chemical potential of oxygen.
    """
    if delta_Hf >= 0 or E_total_O2 >= 0:
        raise ValueError("Formation enthalpy and total energy of O2 must be negative values.")
    
    lower_limit = 0.5 * delta_Hf
    upper_limit = 0.5 * E_total_O2
    return lower_limit, upper_limit

def plot_surface_gibbs_free_energy_all(slabs, bulk_structure, E_slabs, E_bulk, lower_limit, upper_limit, gcd_value):
    """
    Plot the surface Gibbs free energy for all slab structures as a function of the chemical potential of oxygen.

    Parameters:
    slabs (list of StructureData): List of AiiDA StructureData of the slab terminations.
    bulk_structure (StructureData): AiiDA StructureData of the bulk.
    E_slabs (list of float): List of DFT total energies of the slabs.
    E_bulk (float): DFT total energy of the bulk.
    lower_limit (float): Lower limit for the chemical potential of oxygen.
    upper_limit (float): Upper limit for the chemical potential of oxygen.
    gcd_value (int): Value by which the number of atoms needs to be divided to get the minimal stoichiometry.
    """
    mu_O_values = np.linspace(lower_limit, upper_limit, 100)
    plt.figure()

    for i, slab_structure in enumerate(slabs):
        gamma_values = []
        E_slab = E_slabs[i]
        for mu_O in mu_O_values:
            try:
                gamma = surface_total_energy(slab_structure, bulk_structure, E_slab, E_bulk, mu_O, gcd_value)
                gamma_values.append(gamma)
            except ValueError as e:
                gamma_values.append(np.nan)
                print(f"Warning: {e}")

        plt.plot(mu_O_values, gamma_values, label=f'Slab Termination {i+1}')

    plt.xlabel('Chemical Potential of O (eV)')
    plt.ylabel('Surface Gibbs Free Energy (eV/Angstrom^2)')
    plt.title('Surface Gibbs Free Energy vs Chemical Potential of O for Different Terminations')
    plt.legend()
    plt.grid(True)
    plt.show()

def check_minimal_composition(bulk_structure):
    """
    Check if the bulk structure has the minimal composition.

    Parameters:
    bulk_structure (StructureData): AiiDA StructureData of the bulk.

    Returns:
    int: Value by which the number of atoms needs to be divided to get the minimal stoichiometry.
    """
    # Convert StructureData to ASE Atoms
    bulk_structure = bulk_structure.get_ase()
    element_counts = Counter(atom.symbol for atom in bulk_structure)

    # Find the greatest common divisor for the number of atoms of each element
    gcd = np.gcd.reduce(list(element_counts.values()))

    return gcd

def main():
    # Load bulk structure from VASP file
    bulk_structure = StructureData(ase=read("ag2o_bulk.vasp"))
    
    # Load multiple slab structures representing different surface terminations
    slab_files = ["ag2o_slab_termination1.vasp", "ag2o_slab_termination2.vasp", "ag2o_slab_termination3.vasp"]
    slabs = [StructureData(ase=read(slab_file)) for slab_file in slab_files]
    
    # Example values
    E_slabs = [-1000.0, -1200.0, -1100.0]  # DFT total energies of the slabs, in eV
    E_bulk = -500.0  # DFT total energy of the bulk, in eV
    delta_Hf = -4.0  # Formation enthalpy at 0 K and 0 Pa, in eV
    E_total_O2 = -8.0  # Total energy of O2, in eV

    # Calculate limits for the chemical potential of oxygen
    lower_limit, upper_limit = chemical_potential_limits(delta_Hf, E_total_O2)
    print("Chemical potential limits for O:", lower_limit, "to", upper_limit, "eV")

    # Check minimal composition of the bulk structure
    gcd_value = check_minimal_composition(bulk_structure)
    print("Value to divide the number of atoms by to get minimal stoichiometry:", gcd_value)

    # Plot surface Gibbs free energy as a function of the chemical potential of O for all slab terminations
    plot_surface_gibbs_free_energy_all(slabs, bulk_structure, E_slabs, E_bulk, lower_limit, upper_limit, gcd_value)

if __name__ == "__main__":
    main()
