"""
Thermodynamic calculations for AiiDA-TEROS.

This module provides functions for calculating surface energies
and phase diagrams for oxide surfaces.
"""

import numpy as np
from collections import Counter
import pint
from .constants import ureg

def calculate_binary_oxide_limits(bulk_energy, bulk_structure, metal_energy, hf_bulk):
    """
    Calculate the chemical potential limits for a binary oxide.
    
    Parameters:
        bulk_energy (float): Energy of the bulk oxide in eV.
        bulk_structure (StructureData): Bulk oxide structure.
        metal_energy (float): Energy of bulk metal per atom in eV.
        hf_bulk (float): Heat of formation of the bulk oxide in eV.
        
    Returns:
        tuple: (lower_limit, upper_limit) for the oxygen chemical potential.
    """
    element_counts = Counter(atom.symbol for atom in bulk_structure.get_ase())
    gcd_value = np.gcd.reduce(list(element_counts.values()))
    
    # Identify oxygen and metal counts
    y = element_counts['O'] / gcd_value  # Oxygen stoichiometry
    metal_symbols = [symbol for symbol in element_counts.keys() if symbol != 'O']
    if len(metal_symbols) != 1:
        raise ValueError("Binary oxide should have exactly one non-oxygen element.")
    
    metal_symbol = metal_symbols[0]
    x = element_counts[metal_symbol] / gcd_value  # Metal stoichiometry
    
    # Calculate limits
    lower_limit = (1/y) * (bulk_energy - x * metal_energy)
    upper_limit = lower_limit + (1/y) * hf_bulk
    
    return lower_limit, upper_limit

def calculate_binary_surface_energy(
    slab_energy, slab_structure, bulk_energy, bulk_structure, 
    metal_energy, oxygen_chemical_potential
):
    """
    Calculate the surface energy for a binary oxide slab.
    
    Parameters:
        slab_energy (float): Energy of the slab in eV.
        slab_structure (StructureData): Slab structure.
        bulk_energy (float): Energy of the bulk in eV.
        bulk_structure (StructureData): Bulk structure.
        metal_energy (float): Energy of bulk metal per atom in eV.
        oxygen_chemical_potential (float): Chemical potential of oxygen in eV.
        
    Returns:
        float: Surface energy in J/m².
    """
    # Get element counts in bulk structure
    bulk_ase = bulk_structure.get_ase()
    element_counts_bulk = Counter(atom.symbol for atom in bulk_ase)
    gcd_value = np.gcd.reduce(list(element_counts_bulk.values()))
    
    # Identify oxygen and metal in bulk
    metal_symbols = [symbol for symbol in element_counts_bulk.keys() if symbol != 'O']
    if len(metal_symbols) != 1:
        raise ValueError("Binary oxide should have exactly one non-oxygen element.")
    
    metal_symbol = metal_symbols[0]
    x = element_counts_bulk[metal_symbol] / gcd_value  # Metal stoichiometry in bulk
    y = element_counts_bulk['O'] / gcd_value  # Oxygen stoichiometry in bulk
    
    # Get element counts in slab
    slab_ase = slab_structure.get_ase()
    N_metal_slab = len([atom for atom in slab_ase if atom.symbol == metal_symbol])
    N_O_slab = len([atom for atom in slab_ase if atom.symbol == 'O'])
    
    # Calculate surface area
    cell = slab_ase.get_cell()
    A = np.linalg.norm(np.cross(cell[0], cell[1]))  # Area in Å²
    
    # Energy per formula unit
    E_bulk_per_fu = bulk_energy / gcd_value
    
    # Surface energy calculation: Eqn 12 from Reuter and Scheffler, PRB 65, 035406 (2001)
    gamma_ev = (1 / (2 * A)) * (slab_energy - (N_metal_slab / x) * E_bulk_per_fu 
                               + ((y / x) * N_metal_slab - N_O_slab) * oxygen_chemical_potential)
    
    # Convert to J/m²
    gamma_jm2 = gamma_ev * 1.602176634e-19 * 1e20  # eV/Å² to J/m²
    
    return gamma_jm2

def calculate_ternary_surface_energy(
    slab_energy, slab_structure, bulk_energy, bulk_structure,
    metal_energy, delta_mu_metal, delta_mu_oxygen, area=None
):
    """
    Calculate the surface energy for a ternary oxide slab.
    
    Parameters:
        slab_energy (float): Energy of the slab in eV.
        slab_structure (StructureData): Slab structure.
        bulk_energy (float): Energy of the bulk in eV.
        bulk_structure (StructureData): Bulk structure.
        metal_energy (float): Energy of first metal element per atom in eV.
        delta_mu_metal (float): Change in chemical potential of first metal in eV.
        delta_mu_oxygen (float): Change in chemical potential of oxygen in eV.
        area (float, optional): Surface area in Å². If None, calculated from structure.
        
    Returns:
        float: Surface energy in J/m².
    """
    # Get bulk structure composition and stoichiometry
    bulk_ase = bulk_structure.get_ase()
    element_counts_bulk = Counter(atom.symbol for atom in bulk_ase)
    unique_elements = list(element_counts_bulk.keys())
    
    if len(unique_elements) != 3:
        raise ValueError("Ternary oxide should have exactly three elements.")
    
    if 'O' not in unique_elements:
        raise ValueError("Ternary oxide must contain oxygen.")
    
    unique_elements.remove('O')
    first_element, second_element = unique_elements
    
    # Get minimal composition factor
    gcd_value = np.gcd.reduce(list(element_counts_bulk.values()))
    n_first = element_counts_bulk[first_element] / gcd_value
    n_second = element_counts_bulk[second_element] / gcd_value
    n_oxygen = element_counts_bulk['O'] / gcd_value
    
    # Get element counts in slab
    slab_ase = slab_structure.get_ase()
    element_counts_slab = Counter(atom.symbol for atom in slab_ase)
    n_first_slab = element_counts_slab.get(first_element, 0)
    n_second_slab = element_counts_slab.get(second_element, 0)
    n_oxygen_slab = element_counts_slab.get('O', 0)
    
    # Calculate excess atoms (Δ values)
    Delta_Me_A = n_first_slab - n_first * n_second_slab / n_second
    Delta_Me_O = n_oxygen_slab - n_oxygen * n_second_slab / n_second
    
    # Calculate Ψ term (energy term in surface energy formula)
    energy_o2 = oxygen_chemical_potential * ureg.eV
    psi = (slab_energy * ureg.eV - 
           n_second_slab * bulk_energy * ureg.eV / n_second - 
           metal_energy * ureg.eV * Delta_Me_A - 
           0.5 * Delta_Me_O * energy_o2)
    
    # Calculate surface area if not provided
    if area is None:
        cell = slab_ase.get_cell()
        area = np.linalg.norm(np.cross(cell[0], cell[1]))  # Area in Å²
    
    # Calculate surface energy
    a = (cell[0][0]**2 + cell[0][1]**2 + cell[0][2]**2)**0.5 * ureg.angstrom
    b = (cell[1][0]**2 + cell[1][1]**2 + cell[1][2]**2)**0.5 * ureg.angstrom
    
    gamma = (1 / (2 * a * b)) * (psi - 
                                Delta_Me_A * delta_mu_metal * ureg.eV - 
                                Delta_Me_O * delta_mu_oxygen * ureg.eV)
    
    # Convert to J/m²
    gamma_jm2 = gamma.to('J/m^2').magnitude
    
    return gamma_jm2

def generate_phase_diagram_data(terminations, delta_mu_a_range, delta_mu_o_range, precision=500):
    """
    Generate phase diagram data for a ternary oxide.
    
    Parameters:
        terminations (dict): Dictionary with termination names as keys and functions to calculate
                           surface energy as values.
        delta_mu_a_range (tuple): Range of chemical potential for element A (min, max).
        delta_mu_o_range (tuple): Range of chemical potential for oxygen (min, max).
        precision (int): Number of points in each direction.
        
    Returns:
        tuple: (delta_mu_a, delta_mu_o, phase_indices, gamma_min)
               delta_mu_a: 1D array of chemical potential values for element A.
               delta_mu_o: 1D array of chemical potential values for oxygen.
               phase_indices: 2D array of phase indices (which termination is stable).
               gamma_min: 2D array of minimum surface energies.
    """
    delta_mu_a = np.linspace(delta_mu_a_range[0], delta_mu_a_range[1], precision)
    delta_mu_o = np.linspace(delta_mu_o_range[0], delta_mu_o_range[1], precision)
    
    # Initialize arrays
    gamma_min = np.full((precision, precision), np.inf)
    phase_indices = np.zeros((precision, precision), dtype=int)
    
    # Calculate surface energies for each termination
    for term_idx, (term_name, energy_func) in enumerate(terminations.items(), start=1):
        for i, mu_a in enumerate(delta_mu_a):
            for j, mu_o in enumerate(delta_mu_o):
                gamma = energy_func(mu_a, mu_o)
                
                # Update minimum surface energy and phase index
                if gamma < gamma_min[i, j]:
                    gamma_min[i, j] = gamma
                    phase_indices[i, j] = term_idx
    
    return delta_mu_a, delta_mu_o, phase_indices, gamma_min

def calculate_temperature_pressure_map(temperature_range, pressure_range, energies=None, mu_ex=0.0):
    """
    Calculate oxygen chemical potential as a function of temperature and pressure.
    
    Parameters:
        temperature_range (list): List of temperatures in K.
        pressure_range (list): List of logarithmic pressures (ln(p/p⁰)).
        energies (list, optional): List of O₂ free energies at different temperatures.
                                 If None, a simplified model is used.
        mu_ex (float): Excess chemical potential term in eV.
        
    Returns:
        tuple: (temperatures, log_pressures, chemical_potentials)
               temperatures: 1D array of temperatures.
               log_pressures: 1D array of logarithmic pressures.
               chemical_potentials: 2D array of chemical potentials.
    """
    kB = 1.38064852e-23 * ureg["J/K"]  # Boltzmann constant
    mu_ex = mu_ex * ureg.eV  # Excess chemical potential
    
    # Simple model for O₂ free energy if not provided
    if energies is None:
        energies = [-0.15, -0.341, -0.548, -0.765, -0.991, -1.222, -1.460, -1.702, -1.949, -2.199, -2.453, -2.711]
        energies = [e * ureg.eV for e in energies]
        # Ensure temperature range and energies have the same length
        if len(temperature_range) != len(energies):
            raise ValueError("Temperature range and energies must have the same length when using provided energies.")
    
    # Convert temperature range to Quantity objects
    temperatures = [t * ureg.K for t in temperature_range]
    log_pressures = pressure_range
    
    # Calculate chemical potentials for each temperature and pressure
    chemical_potentials = np.zeros((len(temperatures), len(log_pressures)))
    
    for i, (temp, energy) in enumerate(zip(temperatures, energies)):
        for j, log_p in enumerate(log_pressures):
            p = 10**log_p
            mu = (0.5 * (energy + kB * temp * np.log(p)) + mu_ex).to('eV').magnitude
            chemical_potentials[i, j] = mu
    
    return temperature_range, log_pressures, chemical_potentials