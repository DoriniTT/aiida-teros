"""
Ternary oxide plotting functions for AiiDA-TEROS.

This module provides functions for creating plots from
ternary oxide thermodynamic calculation results.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

def plot_ternary_surface_energy(
    results_dict, output_path,
    figsize=(8, 6), dpi=200, y_padding=0.2
):
    """
    Plot surface energies as a function of oxygen chemical potential for different terminations.
    
    Parameters:
        results_dict (dict): Dictionary with termination labels as keys and dictionaries containing
                          'delta_mu' and 'gamma_delta_o' as values.
        output_path (str): Path to save the plot.
        figsize (tuple): Figure size in inches (width, height).
        dpi (int): Figure resolution.
        y_padding (float): Padding for y-axis limits.
        
    Returns:
        tuple: (fig, ax) Matplotlib figure and axis objects.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Get unique termination labels
    labels = sorted(results_dict.keys())
    
    # Create colors for terminations
    colors = sns.color_palette("husl", min(max(len(results_dict), 8), 20))
    
    lines = []
    
    # Extract delta_mu (should be the same for all terminations)
    delta_mu = np.array(results_dict[labels[0]]['delta_mu'])
    
    # Plot each termination
    for n, (label, termination_data) in enumerate(sorted(results_dict.items()), start=1):
        line, = ax.plot(
            delta_mu, 
            termination_data['gamma_delta_o'], 
            label=f'ST-{n}', 
            linewidth=2.5, 
            color=colors[(n - 1) % len(colors)]
        )
        lines.append(line)
    
    # Set axis limits based on min and max values
    all_gamma_values = [value for termination_data in results_dict.values() 
                        for value in termination_data['gamma_delta_o']]
    y_min = min(all_gamma_values) - y_padding
    y_max = max(all_gamma_values) + y_padding
    
    # Font sizes
    fontsize_ticks = 18
    
    # Set labels and ticks
    ax.set_xlabel(r'$\Delta \mu_{\mathrm{O}}$ (eV)', fontsize=fontsize_ticks + 2)
    ax.set_ylabel(r'$\gamma$ (J/m$^2$)', fontsize=fontsize_ticks + 2)
    
    # X-axis ticks
    primary_xticks = np.linspace(min(delta_mu), 0, 5)
    ax.set_xticks(primary_xticks)
    ax.set_xticklabels(primary_xticks, fontsize=fontsize_ticks - 2)
    
    # Y-axis ticks
    yticks = np.round(np.linspace(y_min, y_max, 8), 0)
    ax.set_yticks(yticks)
    ax.set_yticklabels([int(tick) for tick in yticks], fontsize=fontsize_ticks - 2)
    
    # Set axis limits
    ax.set_xlim(-2.0, 0)
    ax.set_ylim(y_min, y_max)
    
    # Add first secondary x-axis (temperature)
    ax2 = ax.twiny()
    ax2.spines['top'].set_position(('outward', 10))
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.linspace(-4, 0, 9))
    T = [1100, 1000, 900, 800, 700, 600, 500, 400, 300]
    ax2.set_xticklabels(T, fontsize=fontsize_ticks - 5)
    ax2.set_xlabel('Temperature (K) @ 1 atm', fontsize=fontsize_ticks - 5)
    
    # Add second secondary x-axis (pressure)
    ax3 = ax.twiny()
    ax3.spines['top'].set_position(('outward', 60))
    ax3.set_xlim(ax.get_xlim())
    ax3.set_xticks(np.linspace(-4, 0, 12))
    ax3.set_xticklabels([None, None, r'$10^{-27}$', r'$10^{-24}$', r'$10^{-21}$', r'$10^{-18}$', 
                         r'$10^{-15}$', r'$10^{-12}$', r'$10^{-9}$', r'$10^{-6}$', 
                         r'$10^{-3}$', r'$1$'], fontsize=fontsize_ticks - 5)
    ax3.set_xlabel('Pressure (atm) @ 300 K', fontsize=fontsize_ticks - 5)
    
    # Create the legend outside of the plot area
    ax.legend(lines, [f'ST-{i}' for i in range(1, len(lines) + 1)], 
              loc='center left', bbox_to_anchor=(1, 0.5), 
              fontsize=13, borderaxespad=0.)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight")
    
    return fig, ax

def plot_phase_diagram(
    results_dict, unique_symbols, lim_delta_ag, lim_delta_o,
    output_path, figsize=(13, 6), dpi=200
):
    """
    Plot the phase diagram for a ternary oxide.
    
    Parameters:
        results_dict (dict): Dictionary with termination labels as keys and dictionaries containing
                          'delta_ag', 'delta_o', and 'gamma_delta_ag_delta_o' as values.
        unique_symbols (list): List of unique element symbols in the structure.
        lim_delta_ag (float): Upper limit for delta_ag (first metal).
        lim_delta_o (float): Upper limit for delta_o (oxygen).
        output_path (str): Path to save the plot.
        figsize (tuple): Figure size in inches (width, height).
        dpi (int): Figure resolution.
        
    Returns:
        tuple: (fig, ax) Matplotlib figure and axis objects.
    """
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=figsize, dpi=dpi,
        gridspec_kw={'wspace': 0, 'width_ratios': [2, 1]}
    )
    
    # Sort labels by termination number
    labels = sorted(results_dict.keys(), key=lambda x: int(x.replace('relax_slab_', '')))
    
    # Get arrays from first termination (should be the same for all)
    delta_ag = results_dict[labels[0]]['delta_ag']
    delta_o = results_dict[labels[0]]['delta_o']
    
    # Initialize arrays for phase diagram
    gamma_min = np.full((len(delta_ag), len(delta_o)), np.inf)
    phase = np.zeros((len(delta_ag), len(delta_o)), dtype=int)
    
    # Calculate minimum energy at each point
    for label in labels:
        gamma = np.array(results_dict[label]['gamma_delta_ag_delta_o'])
        mask = gamma < gamma_min
        gamma_min[mask] = gamma[mask]
        phase[mask] = int(label.replace('relax_slab_', ''))
    
    # Create a custom colormap from the seaborn husl palette
    husl_palette = sns.color_palette("husl", len(labels))
    cmap = LinearSegmentedColormap.from_list("custom_husl", husl_palette, N=len(labels))
    
    # Plot the phase diagram
    contour = ax1.contourf(delta_ag, delta_o, phase.T, levels=len(labels), cmap=cmap)
    
    # Configure the left plot
    ax1.spines['top'].set_visible(True)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(True)
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    ax1.yaxis.set_ticks_position('left')
    ax1.set_ylabel(r'$\Delta \mu_{\mathrm{O}}$ (eV)', fontsize=16)
    ax1.set_xlabel(rf'$\Delta \mu_{{\mathrm{{{unique_symbols[0]}}}}}$ (eV)', fontsize=16, labelpad=10)
    ax1.set_ylim(-2, 0)
    ax1.set_xlim(lim_delta_ag, 0)
    ax1.yaxis.set_visible(False)
    
    # Uniform tick and label sizes
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax1.tick_params(axis='both', which='minor', labelsize=14)
    
    # Plot Delta Mu vs Pressure on the right (ax2)
    log_pressures = np.linspace(-30, 10, 10)
    T = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
    G = [-0.15, -0.341, -0.548, -0.765, -0.991, -1.222, -1.460, -1.702, -1.949, -2.199, -2.453, -2.711]
    mu_ex = -0.03
    
    # Calculate and plot temperature-pressure relationships
    colors = sns.color_palette("husl", 12)
    for color, t, g in zip(colors, T, G):
        delta_mu = [
            0.5 * (g + 1.38064852e-23 * t * 8.617333262e-5 * np.log(10**p)) + mu_ex
            for p in log_pressures
        ]
        ax2.plot(log_pressures, delta_mu, label=f'T = {t} K', linewidth=2.5, color=color)
    
    # Configure the right plot
    ax2.set_xlabel(r'ln(p/p$^0$)', fontsize=16, labelpad=10)
    ax2.set_ylabel(r'$\Delta \mu_{\mathrm{O}}$ (eV)', fontsize=16)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='minor', labelsize=14)
    ax2.set_ylim(-2, 0)
    ax2.set_xlim(-30, 5)
    xticks = [-25, -20, -15, -10, -5, 0, 5]
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticks, fontsize=16)
    
    # Create a mapping from phase number to color based on label order
    label_to_color = {
        int(label.replace('relax_slab_', '')): husl_palette[i]
        for i, label in enumerate(labels)
    }
    
    # Create a legend using the mapping
    unique_phases = np.unique(phase)
    legend_patches = [
        mpatches.Patch(color=label_to_color[phase_val], label=f'ST-{phase_val}')
        for phase_val in unique_phases
    ]
    ax1.legend(handles=legend_patches, loc='lower left', bbox_to_anchor=(0.05, 0.05))
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Adjust layout and save the figure
    plt.tight_layout()
    fig.subplots_adjust(wspace=0)
    fig.savefig(output_path)
    
    return fig, ax1, ax2