# AiiDA-TEROS: Automating Surface Thermodynamics

Author: Dr. Thiago Trevizam Dorini (thiagotd@unicamp.br)

Contact this email for support, bug reports, or information.

## Introduction

Welcome to the **AiiDA-TEROS** repository!

**TEROS** stands for **Thermodynamics of Oxide Surfaces**, and this project is dedicated to automating the study of surface thermodynamics for oxide materials using the powerful combination of [AiiDA](https://www.aiida.net/) and [VASP](https://www.vasp.at/).

Understanding the thermodynamic properties of oxide surfaces is crucial in fields such as catalysis, corrosion science, and materials engineering. Oxide surfaces often exhibit complex behaviors due to their interactions with the environment, especially in the presence of oxygen and under varying temperature and pressure conditions. Traditional computational studies of these systems can be time-consuming and labor-intensive, requiring significant manual setup and data handling.

**AiiDA-TEROS** addresses these challenges by providing a comprehensive and fully automated workflow that streamlines the entire process—from pre-processing and calculations to post-processing and analysis. By leveraging AiiDA's workflow management capabilities and VASP's robust computational engine, researchers can focus on interpreting results rather than managing computational tasks.

### Why Use AiiDA-TEROS?

- **Efficiency**: Automate repetitive and complex tasks involved in surface thermodynamics calculations, saving time and reducing the possibility of human error.
- **Scalability**: Easily extend the workflow to study various ternary oxide materials and explore multiple surface orientations and terminations.
- **Reproducibility**: Ensure consistent computational procedures, which is essential for validating results and comparing studies across different systems.
- **Customization**: Tailor calculation parameters to fit specific research needs while maintaining the integrity of the automated workflow.

### Key Features

- **Automated Surface Generation**: Generate symmetric surface terminations from any bulk structure and crystallographic orientation automatically.
- **Efficient Relaxation Calculations**: Perform relaxation calculations on all generated slabs using VASP with predefined settings, ensuring consistency and accuracy.
- **Surface Thermodynamics Analysis**: Compute surface Gibbs free energies and construct surface phase diagrams using _ab initio_ atomistic thermodynamics methods.
- **Stable Surface Identification**: Automatically identify and output the most thermodynamically stable surface terminations based on comprehensive energy analyses.
- **Integration with AiiDA**: Utilize AiiDA's powerful data management and workflow automation features to seamlessly track computations and results.

By automating these complex tasks, **AiiDA-TEROS** accelerates the research process, enabling you to explore new materials and surface phenomena more efficiently. Whether you are investigating catalytic surfaces, studying corrosion resistance, or designing novel materials, this tool provides a robust foundation for your computational thermodynamics studies.

---

This README will guide you through the prerequisites, installation, and usage of AiiDA-TEROS. Let's get started!

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Preparing Input Files](#1-preparing-input-files)
  - [2. Running the Workflow](#2-running-the-workflow)
  - [3. Monitoring and Retrieving Results](#3-monitoring-and-retrieving-results)
- [Code Structure](#code-structure)
- [Output Explanation](#output-explanation)
- [Acknowledgements](#acknowledgements)

---

## Prerequisites

Before using AiiDA-TEROS, ensure you have the following software and packages installed:

- **Python 3.7+**
- **AiiDA** (v1.0 or higher)
- **VASP** (and a valid license)
- **AiiDA-VASP** plugin
- **ASE** (Atomic Simulation Environment)
- **pymatgen** (Python Materials Genomics)
- **NumPy**, **Matplotlib**, **Seaborn**
- **pint** (for unit handling)
- **Additional Python packages**: `jsonschema`, `tabulate`, `pyyaml`, `click`

### Note on VASP

VASP is proprietary software and requires a valid license to use. Ensure that you have access to VASP and that it is correctly set up with AiiDA on your computing resources.

---

## Installation

### Option 1: Install from Source

1. **Clone the Repository**

   ```bash
   git clone git@github.com:DoriniTT/aiida-teros.git
   cd aiida-teros
   ```

2. **Install the Package**

   ```bash
   pip install -e .
   ```

   This will install the package in development mode, allowing you to modify the code and see the changes immediately.

### Option 2: Set Up a Python Virtual Environment

It is recommended to set up a Python virtual environment to keep your dependencies isolated and avoid conflicts with other Python projects.

1. **Create and Activate a Virtual Environment**

   ```bash
   python3 -m venv aiida_teros_env
   source aiida_teros_env/bin/activate
   ```

2. **Install AiiDA-TEROS and its Dependencies**

   ```bash
   pip install -e /path/to/aiida-teros
   ```

3. **Configure AiiDA and AiiDA-VASP**

   Ensure that your computer is properly configured for using AiiDA and AiiDA-VASP. Please refer to the official [AiiDA](https://aiida.readthedocs.io/) and [AiiDA-VASP documentation](https://aiida-vasp.readthedocs.io/) for configuration instructions.

---

## Usage

### 1. Preparing Input Files

AiiDA-TEROS requires a configuration file in YAML format to define all input parameters for the workflow. Here's how to prepare the necessary files:

#### Configuration File

Create a `config.yaml` file with the following sections:

```yaml
# Paths and File Names
bulk_structure_path: "path/to/bulk.vasp"          # Path to bulk oxide structure
bulk_metal_structure_path: "path/to/metal.vasp"   # Path to bulk metal structure
potential_family: "PBE"                           # Potential family (e.g., 'PBE')
code_label: "vasp@computer"                       # Label of the VASP code in AiiDA

# Thermodynamic Parameters
thermodynamic_parameters:
  hf_bulk: -1.884                                 # Heat of formation for bulk in eV

# Total Energies
total_energies:
  total_energy_first_element: -2.8289             # Energy of the first element (e.g., Ag)
  total_energy_o2: -9.82                          # Energy for O2 molecule

# INCAR Parameters for Bulk Relaxations
incar_parameters_bulk:
  incar:
    ISMEAR: 0
    SIGMA: 0.01
    # ... other INCAR parameters

# INCAR Parameters for Bulk Metal Relaxations
incar_parameters_bulk_metal:
  incar:
    ISMEAR: 0
    SIGMA: 0.01
    # ... other INCAR parameters

# INCAR Parameters for Slab Relaxations
incar_parameters_slab:
  incar:
    ISMEAR: 0
    SIGMA: 0.01
    # ... other INCAR parameters

# Workflow Settings
workflow_settings:
  kpoints_precision: 0.3                          # K-points mesh density

# Potential Mapping
potential_mapping:
  Ag: "Ag"
  O: "O"
  # ... mappings for other elements

# Parser Settings
parser_settings:
  parser_settings:
    add_energies: true
    add_trajectory: true
    # ... other parser settings

# Computer Options
computer_options:
  resources:
    num_machines: 1
    num_cores_per_machine: 40
  queue_name: "queue_name"

# Slab Generation Parameters
slab_parameters:
  miller_indices: [1, 1, 0]                       # Miller indices
  min_slab_thickness: 10.0                        # Minimum slab thickness in Å
  vacuum: 15.0                                    # Vacuum thickness in Å

# Optional: Specific Terminations
# terminations:
#   slab1: "path/to/termination1.vasp"
#   slab2: "path/to/termination2.vasp"
```

#### Structure Files

Prepare the following structure files in VASP POSCAR format:

1. **Bulk Oxide Structure**: The bulk structure of the oxide material.
2. **Bulk Metal Structure**: The bulk structure of the metal component.
3. **Optional Termination Structures**: Specific surface terminations, if you want to use predefined slabs instead of automatically generated ones.

### 2. Running the Workflow

AiiDA-TEROS provides a command-line interface for running the workflow. There are two ways to run the workflow:

#### Using the Command-Line Interface (Recommended)

```bash
# Run the workflow
aiida-teros run /path/to/config.yaml

# Validate a configuration file without running the workflow
aiida-teros validate /path/to/config.yaml
```

#### Using the Python Script (Backward Compatibility)

For backward compatibility, you can also use the `run_aiida.py` script:

```bash
# Ensure that the AiiDA daemon is running
verdi daemon start

# Run the workflow
cd /path/to/work/directory  # Directory containing config.yaml
python run_aiida.py
```

### 3. Monitoring and Retrieving Results

#### Monitor the Workflow

Use AiiDA's command-line tools to monitor the progress:

```bash
# List running processes
verdi process list

# Show details of a specific process
verdi process show <PK>
```

Replace `<PK>` with the process ID printed after submitting the workflow.

#### Retrieve Outputs

The workflow generates several outputs:

- **Relaxed Structures**: Bulk and slab structures after relaxation.
- **Thermodynamic Data**: Surface energies, chemical potentials, and phase diagrams.
- **LaTeX Tables**: Tables with parameters used for calculating surface energies.
- **Plots**: Surface energy plots and phase diagrams saved in the `thermo_results/` directory.

## Code Structure

AiiDA-TEROS follows a modular architecture to improve maintainability and extensibility:

```
aiida_teros/
├── __init__.py             # Package initialization
├── cli.py                  # Command-line interface
├── core/                   # Core components
├── schemas/                # JSON schemas for configuration validation
│   └── config_schema.py    # Schema for config.yaml
├── utils/                  # Utility functions
│   ├── constants.py        # Physical constants and unit conversion
│   ├── io.py               # I/O operations and file handling
│   ├── structure.py        # Structure manipulation and analysis
│   ├── thermodynamics.py   # Thermodynamic calculations
│   └── vasp.py             # VASP calculation setup
├── plotting/               # Plotting functions
│   ├── binary.py           # Binary oxide plotting
│   └── ternary.py          # Ternary oxide plotting
└── workflows/              # AiiDA workflows
    └── teros.py            # Main TEROS workflow
```

## Output Explanation

The **AiiDA-TEROS** workflow generates several key outputs organized into categories. These outputs include relaxed structures and associated data, essential for further thermodynamic analysis. Below is an overview:

### 1. Bulk Relaxation

Outputs from the relaxation of the bulk structure:

- **structure** (`StructureData`): The relaxed bulk structure.
- **trajectory** (`TrajectoryData`): Contains positions, velocities, and forces over the relaxation steps.
- **energies** (`ArrayData`): Total energies at each step of the relaxation.
- **forces** (`ArrayData`): Forces acting on each atom during relaxation.
- **kpoints** (`KpointsData`): The k-point mesh used in the calculation.
- **misc** (`Dict`): Miscellaneous information related to the relaxation, such as convergence parameters.
- **remote_folder** (`RemoteData`): Directory where the calculations were performed.
- **retrieved** (`FolderData`): Retrieved output files from the calculations.

### 2. Slab Relaxations

Multiple symmetric slab terminations are generated from the bulk structure based on the specified Miller indices. For each slab termination, the outputs include:

- **structure** (`StructureData`): The relaxed slab structure.
- **trajectory** (`TrajectoryData`): Relaxation process details.
- **energies** (`ArrayData`): Total energies at each relaxation step.
- **forces** (`ArrayData`): Forces acting on each atom.
- **kpoints** (`KpointsData`): The k-point mesh used.
- **misc** (`Dict`): Additional information from the relaxation process.
- **remote_folder** (`RemoteData`): Directory where calculations were performed.
- **retrieved** (`FolderData`): Retrieved calculation files.

### 3. Stable Structures

After relaxing all slab terminations, the workflow identifies the most thermodynamically stable surfaces.

### 4. Generated Plots and Tables

#### For Binary Oxides (e.g., Ag₂O)

- **Surface Energy Plot**: Shows surface Gibbs free energy (γ) as a function of oxygen chemical potential (Δμ_O).
- **LaTeX Table**: Contains parameters used for calculating surface energies.

Output files:
- `thermo_results/binary/surface_energy_plot.pdf`
- `thermo_results/binary/termination_parameters_table.tex`

#### For Ternary Oxides (e.g., Ag₂MoO₄)

- **Surface Energy Plot**: Shows surface Gibbs free energy vs. oxygen chemical potential.
- **Phase Diagram**: Illustrates stable terminations as a function of chemical potentials.
- **LaTeX Table**: Contains parameters for ternary oxide surface energy calculations.

Output files:
- `thermo_results/ternary/surface_free_energies.pdf`
- `thermo_results/ternary/surface_phase_diagram.pdf`
- `thermo_results/ternary/ternary_parameters_table.tex`

## Acknowledgements

We gratefully acknowledge the financial support from the Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP) under project CEPID, multiusuário, grant numbers 2013/07296-2, 2016/23891-6, 2017/26105-4, and 2023/03447-8.

---

*For any questions or issues, please contact Dr. Thiago Trevizam Dorini at thiagotd@unicamp.br.*