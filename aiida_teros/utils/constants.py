"""
Physical constants and unit conversion for AiiDA-TEROS.

This module provides physical constants and unit conversion factors
for thermodynamic calculations.
"""

import pint

# Initialize unit registry
ureg = pint.UnitRegistry()

# Constants
KB = 1.38064852e-23 * ureg.joule / ureg.kelvin  # Boltzmann constant
EV_TO_JOULE = 1.602176634e-19  # Conversion from eV to J
JOULE_TO_EV = 1.0 / EV_TO_JOULE  # Conversion from J to eV
ANGSTROM_TO_METER = 1e-10  # Conversion from Angstrom to meter
METER_TO_ANGSTROM = 1e10  # Conversion from meter to Angstrom

# Energy conversion
def ev_to_jm2(energy_ev, area_angstrom2):
    """
    Convert energy from eV/Å² to J/m².
    
    Parameters:
        energy_ev (float): Energy in eV.
        area_angstrom2 (float): Area in Å².
        
    Returns:
        float: Energy in J/m².
    """
    return energy_ev * EV_TO_JOULE / (area_angstrom2 * ANGSTROM_TO_METER**2)

def jm2_to_ev(energy_jm2, area_angstrom2):
    """
    Convert energy from J/m² to eV/Å².
    
    Parameters:
        energy_jm2 (float): Energy in J/m².
        area_angstrom2 (float): Area in Å².
        
    Returns:
        float: Energy in eV/Å².
    """
    return energy_jm2 * JOULE_TO_EV * (area_angstrom2 * ANGSTROM_TO_METER**2)

# Temperature dependence of oxygen chemical potential
def oxygen_chemical_potential(temperature, pressure=1.0, reference_energy=0.0):
    """
    Calculate the chemical potential of oxygen as a function of temperature and pressure.
    
    Parameters:
        temperature (float): Temperature in K.
        pressure (float): Pressure in atm.
        reference_energy (float): Reference energy in eV.
        
    Returns:
        float: Chemical potential in eV.
    """
    # This is a simplified model; for a more accurate model, you need to include
    # the contribution from vibrational modes, etc.
    t_kelvin = temperature * ureg.kelvin
    p_atm = pressure * ureg.atm
    
    # Standard free energy change with temperature (simplified model)
    # Values from NIST-JANAF tables could be used for more accuracy
    delta_g = reference_energy
    
    # Pressure contribution to chemical potential (ideal gas approximation)
    mu = delta_g + (KB * t_kelvin * ureg.log(p_atm / ureg.atm)).to('eV')
    
    return mu.magnitude