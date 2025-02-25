"""
Plotting package for AiiDA-TEROS.

This package provides functions for creating plots from
thermodynamic calculation results.
"""

from .binary import plot_binary_surface_energy
from .ternary import plot_ternary_surface_energy, plot_phase_diagram

__all__ = [
    'plot_binary_surface_energy',
    'plot_ternary_surface_energy',
    'plot_phase_diagram'
]