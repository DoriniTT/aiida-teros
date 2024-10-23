# AiiDA-TEROS: Automating Surface Thermodynamics

## Introduction

Welcome to the **AiiDA-TEROS** repository!

**TEROS** stands for **Thermodynamics of Oxide Surfaces**, and this project is dedicated to automating the study of surface thermodynamics for oxide materials using the powerful combination of [AiiDA](https://www.aiida.net/) and [VASP](https://www.vasp.at/).

Understanding the thermodynamic properties of oxide surfaces is crucial in fields such as catalysis, corrosion science, and materials engineering. Oxide surfaces often exhibit complex behaviors due to their interactions with the environment, especially in the presence of oxygen and under varying temperature and pressure conditions. Traditional computational studies of these systems can be time-consuming and labor-intensive, requiring significant manual setup and data handling.

**AiiDA-TEROS** addresses these challenges by providing a comprehensive and fully automated workflow that streamlines the entire process—from pre-processing and calculations to post-processing and analysis. By leveraging AiiDA's workflow management capabilities and VASP's robust computational engine, researchers can focus on interpreting results rather than managing computational tasks.

### Why Use AiiDA-TEROS?

- **Efficiency**: Automate repetitive and complex tasks involved in surface thermodynamics calculations, saving time and reducing the possibility of human error.
- **Scalability**: Easily extend the workflow to study various ternary oxides materials (for the moment!), and explore multiple surface orientations and terminations.
- **Reproducibility**: Ensure consistent computational procedures, which is essential for validating results and comparing studies across different systems.
- **Customization**: Tailor calculation parameters to fit specific research needs while maintaining the integrity of the automated workflow.

### Key Features

- **Automated Surface Generation**: Generate symmetric surface terminations from any given bulk structure and crystallographic orientation without manual intervention.
- **Efficient Relaxation Calculations**: Perform relaxation calculations on all generated slabs using VASP with predefined settings, ensuring consistency and accuracy.
- **Surface Thermodynamics Analysis**: Compute surface Gibbs free energies and construct surface phase diagrams using _ab initio_ atomistic thermodynamics methods.
- **Stable Surface Identification**: Automatically identify and output the most thermodynamically stable surface terminations based on comprehensive energy analyses.
- **Integration with AiiDA**: Utilize AiiDA's powerful data management and workflow automation features to keep track of computations and results seamlessly.

By automating these complex tasks, **AiiDA-TEROS** accelerates the research process, enabling you to explore new materials and surface phenomena more efficiently. Whether you are investigating catalytic surfaces, studying corrosion resistance, or designing novel materials, this tool provides a robust foundation for your computational thermodynamics studies.

---

This README will guide you through the prerequisites, installation, and usage of AiiDA-TEROS. Let's get started!

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Setting Up the Environment](#1-setting-up-the-environment)
  - [2. Preparing Input Files](#2-preparing-input-files)
  - [3. Running the WorkChain](#3-running-the-workchain)
  - [4. Monitoring and Retrieving Results](#4-monitoring-and-retrieving-results)
- [Output Explanation](#output-explanation)

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
- **Additional Python packages**: `subprocess`, `shutil`, `os`, `time`

### Note on VASP

VASP is a proprietary software and requires a valid license to use. Ensure that you have access to VASP and that it is correctly set up with AiiDA on your computing resources.

---

## Installation

1. **Clone the Repository**

   ```bash
   git clone git@github.com:DoriniTT/aiida-teros.git
   cd aiida-teros
   ```
2. **Set Up a Python Virtual Environment (Optional but Recommended)**

It is recommended to set up a Python virtual environment to keep your dependencies isolated and avoid conflicts with other Python projects. You can create and activate a virtual environment by running the following commands:

  ```bash
  python3 -m venv aiida_teros
  source aiida_teros/bin/activate
  ```
3. **Install Python Dependencies**

After setting up the virtual environment, install the necessary Python packages by running:

  ```bash
  pip install aiida-core aiida-vasp ase pymatgen numpy matplotlib seaborn pint
  ```

4. **Set Up AiiDA and AiiDA-VASP**

Ensure that your computer is properly configured for using AiiDA and AiiDA-VASP. Please refer to the official [AiiDA](https://aiida.readthedocs.io/) and [AiiDA-VASP documentation](https://aiida-vasp.readthedocs.io/) for configuration instructions.
This repository does **not** include instructions for configuring AiiDA or AiiDA-VASP.

5. **Update PYTHONPATH**

You need to add the path to `AiiDA_teros.py` in your PYTHONPATH. You can do this by running the following command in your terminal:

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/your/repo
```
Replace /path/to/your/repo with the actual path to the directory containing AiiDA_complete_thermo.py.

## Usage

The AiiDA-TEROS consists of two main scripts:

- AiiDA_teros.py: Defines the AiiDA-TEROS class, implementing the workflow.

- run_aiida.py: Submission script to run the AiiDA-TEROS with user-defined inputs.

1. **Preparing input files**

   - **Bulk Structure File**

Prepare a bulk structure file in VASP POSCAR format (e.g., bulk_structure.vasp). Place it in the same repertory as the main python scripts.

   - **Define INCAR parameters**

Customize the INCAR parameters for both bulk and slab relaxation in the submission script (run_aiida.py):

```python
# INCAR Parameters for Bulk Relaxations
INCAR_PARAMETERS_BULK = {'incar': {
    'ISMEAR': 0,
    'ENCUT': 550,
    'ISIF': 3,
    'IBRION': 2,
    'NSW': 100,
    'EDIFFG': -0.01,
    'PREC': 'Accurate',
    'EDIFF': 1e-5,
    }
}

# INCAR Parameters for Slab Relaxations
INCAR_PARAMETERS_SLAB = {'incar': {
    'ISMEAR': 0,
    'ENCUT': 550,
    'ISIF': 2,
    'IBRION': 2,
    'NSW': 1000,
    'EDIFFG': -0.05,
    'PREC': 'Accurate',
    'EDIFF': 1e-5,
    }
}
```

- **Set Workflow inputs**

    In run_aiida.py, configure the inputs for the workflow, such as:

- **Force Cutoff**: Convergence criterion for forces.

- **K-Points Precision**: The k-points mesh is automatically set using the `set_kpoints_mesh_from_density` method from AiiDA, ensuring an appropriate density of k-points for accurate Brillouin zone sampling based on the structure's size and geometry.

- **Potential Mapping and Family**: The potential family and mapping are set using an auxiliary function that defines the correct potentials for the calculation. For example, to use the 'PBE' family and map potentials for silver (Ag) and oxygen (O), you would configure it as follows in the code:

   ```python
   # Set potential family and mapping
   builder.potential_family = Str('PBE')  # Using the PBE potential family
   builder.potential_mapping = {'Ag': 'Ag', 'O': 'O', 'P': 'P'}  # Map potentials for Ag, O, and P elements
   ```

- **Slab Generation Parameters**: Set the Miller indices, slab thickness, and vacuum spacing to define the slab’s orientation, size, and separation between periodic images.

   ```python
   SLAB_PARAMETERS = {
       'miller_indices': [1, 1, 0],  # Example: [1, 1, 0]
       'min_slab_thickness': 10.0,   # Minimum slab thickness in Å
   }
   ```

- **Thermodynamic Parameters**:

  #### For Binary Oxides (e.g., Ag₂O)
  For binary oxides, you only need the enthalpy of formation (ΔH_f) in eV per formula unit. This can be calculated in advance or taken from experimental or theoretical sources. The thermodynamic model used the equations demonstrated in the following paper:
  [1] K. Reuter e M. Scheffler, “Composition, structure, and stability of RuO₂(110) as a function of oxygen pressure”, Phys. Rev. B, vol. 65, nº 3, p. 035406, dez. 2001, doi: 10.1103/PhysRevB.65.035406.

  #### For Ternary Oxides (e.g., Ag₂MoO₄)
  For ternary oxides, there are additional considerations:
  - Ensure the order of elements in the input file (e.g., POSCAR, CIF) is consistent with the standard notation (e.g., Ag₂MoO₄).
  - Besides the enthalpy of formation (ΔH_f), you also need the total energies (E_tot) calculated with DFT for the first element (in this example, Ag) and for oxygen. These should be obtained using the same functional as used in the TEROS workchain. The thermodynamic model used the equations demonstrated in the following paper:
  [1] E. Heifets, J. Ho, e B. Merinov, “Density functional simulation of the BaZrO₃ (011) surface structure”, Phys. Rev. B, vol. 75, nº 15, p. 155431, abr. 2007, doi: 10.1103/PhysRevB.75.155431.

  In future updates, we plan to automate the determination of the enthalpy of formation and the total energies of common elements within the code.

- **Computer Options**: Configure computer options such as resources and queue name:

   ```python
   COMPUTER_OPTIONS = {
       'resources': {
           'num_machines': 1,
           'num_cores_per_machine': 40
       },
       'queue_name': 'par40',
   }
   ```

## Running the WorkChain

Execute the submission script:

```bash
python run_aiida.py
```

This script performs the following actions:

- Starts and restarts the AiiDA daemon to ensure it's running.
- Loads the bulk structure and prepares it as an AiiDA [`StructureData`](https://aiida.readthedocs.io/projects/aiida-core/en/latest/topics/data_types.html#structuredata) node.
- Sets up all the necessary inputs for the AiiDA-TEROS.
- Submits the work chain to the AiiDA daemon for execution.
- Writes the process ID (PK) to `pks.txt` for future reference.

## Monitoring and Retrieving Results

- **Monitor the WorkChain**

  Use AiiDA's command-line tools to monitor the progress:
  ```bash
  verdi process list      # List running processes
  verdi process show PK   # Show details of the process with ID PK
  ```

- **Retrieve outputs**

  Upon completion, the workflow generates outputs including:

  - Relaxed bulk and slab structures.
  - Surface Gibbs free energies and phase diagrams.
  - Plots saved in the specified directory (`thermo_results/binary` or `thermo_results/ternary`).

- **Visualize results**

  The generated plots can be found in the appropriate subdirectory (`thermo_results/binary` or `thermo_results/ternary`). These include:

  - **surface_free_energies.pdf**

## Output Explanation

The **MasterThermoWorkChain** generates several key outputs organized into categories. These outputs include relaxed structures and associated data, essential for further thermodynamic analysis. Below is an overview:

### 1. Bulk Relaxation:
   Outputs from the relaxation of the bulk structure:
   - **structure**: The relaxed bulk structure.
   - **misc**: Miscellaneous information related to the relaxation.
   - **remote_folder**: Directory where the calculations were performed.
   - **retrieved**: Retrieved output files from the calculations.

### 2. Slab Relaxations:

Multiple symmetric slab terminations are generated from the bulk structure based on the specified Miller indices. These slabs represent different possible surface configurations of the material. Each slab is then relaxed to find its minimum energy configuration. The number of slabs (`ST_1`, `ST_2`, ...) depends on the material's symmetry and the possible terminations.

For each slab termination, the outputs include:

- **structure** (`StructureData`): The relaxed slab structure, containing the atomic positions and cell parameters after relaxation.
- **misc** (`Dict`): Additional information from the relaxation process, such as total energies, forces, stress tensors, and convergence details.
- **remote_folder** (`RemoteData`): The directory on the remote computer where the calculation was performed, useful for accessing raw data if needed.
- **retrieved** (`FolderData`): The files retrieved from the calculation, including output logs, charge densities, and other relevant data.

These outputs allow you to analyze the relaxed structures of different slab terminations and compare their properties.

### 3. Stable Structures:

After relaxing all slab terminations, the work chain identifies the most thermodynamically stable surface terminations based on calculated surface energies. These stable structures are critical for understanding surface phenomena and predicting material behavior under various conditions.

For each identified stable slab (`relax_slab_1`, `relax_slab_2`, ...), the outputs are:

- **structure** (`StructureData`): The relaxed stable slab structure, representing the most energetically favorable configuration.
- **misc** (`Dict`): Detailed information from the relaxation and analysis, including surface energy calculations and relevant thermodynamic quantities.
- **remote_folder** (`RemoteData`): The directory where the stable slab calculations were performed.
- **retrieved** (`FolderData`): Retrieved files from the calculation, which may include data specific to surface thermodynamics analysis.

### 4. Generated Plots for Surface Gibbs Free Energies vs. Chemical Potentials

#### For Binary Oxides (e.g., Ag₂O)
For binary oxides, a plot is generated showing the **Surface Gibbs free energy (γ) as a function of the chemical potential of oxygen (Δμ_O)**. This plot helps visualize how changes in the chemical potential of oxygen affect the surface stability.

The generated plot is saved in the directory `thermo_results/binary`.

#### For Ternary Oxides (e.g., Ag₂MoO₄)
For ternary oxides, two plots are generated:

- **Surface Gibbs free energy (γ) vs. chemical potential of oxygen (Δμ_O)**: Similar to the binary oxide case, this plot shows the relationship between the surface Gibbs free energy and the variation in the chemical potential of oxygen.
- **Surface Phase Diagram**: This plot shows the variation of the chemical potential of the first element (e.g., Ag, Δμ_Ag) and the variation of the chemical potential of oxygen (Δμ_O). The most stable surface terminations are identified in this diagram, providing insights into how both elements contribute to surface stability under different conditions.

The generated plots are saved in the directory `thermo_results/ternary`.

### Common Output Descriptions:
For all the categories above, the outputs are stored in consistent AiiDA node types:
- **structure** (`StructureData`): The relaxed structure (bulk or slab).
- **misc** (`Dict`): Additional data and metadata from the relaxation process.
- **remote_folder** (`RemoteData`): The directory on the remote computer where calculations were performed.
- **retrieved** (`FolderData`): The files retrieved from the remote calculation, including output data necessary for analysis.

These outputs provide all the essential information required for further thermodynamic analyses, such as calculating surface Gibbs free energies and constructing phase diagrams. The consistent structure of the outputs ensures that they can be easily accessed and processed, regardless of the specific material or the number of slabs generated.

## Acknowledgements

We gratefully acknowledge the financial support from the Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP) under project CEPID, multiusuário, grant number 2013/07296-2, 2016/23891-6, 2017/26105-4, and 2023/03447-8.
