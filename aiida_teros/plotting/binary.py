"""
Binary oxide plotting functions for AiiDA-TEROS.

This module provides functions for creating plots from
binary oxide thermodynamic calculation results.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_binary_surface_energy(
    gammas, mu_o_values, output_path,
    figsize=(10, 6), title="Surface Energy vs O Chemical Potential"
):
    """
    Plot surface energies as a function of oxygen chemical potential for different terminations.
    
    Parameters:
        gammas (dict): Dictionary with termination names as keys and dictionaries containing
                      'gamma_lower' and 'gamma_upper' as values.
        mu_o_values (list): List of oxygen chemical potential values (lower and upper limits).
        output_path (str): Path to save the plot.
        figsize (tuple): Figure size in inches (width, height).
        title (str): Plot title.
        
    Returns:
        tuple: (fig, ax) Matplotlib figure and axis objects.
    """
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create color cycle for different terminations
    colors = plt.cm.tab10(np.linspace(0, 1, len(gammas)))
    
    # Plot lines for each termination
    for (term_name, term_data), color in zip(gammas.items(), colors):
        # Create points for the line
        x_points = mu_o_values
        y_points = [term_data['gamma_lower'], term_data['gamma_upper']]
        
        # Plot the line
        ax.plot(x_points, y_points, '-', label=term_name, color=color, linewidth=2)
        
        # Add points to show exact values
        ax.plot(x_points, y_points, 'o', color=color, markersize=6)
    
    # Customize the plot
    ax.set_xlabel('O Chemical Potential (eV)', fontsize=12)
    ax.set_ylabel('Surface Energy (J/m²)', fontsize=12)
    ax.set_title(title, fontsize=14, pad=15)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Add legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Adjust layout to prevent cutting off the legend
    plt.tight_layout()
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save the figure
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    
    return fig, ax