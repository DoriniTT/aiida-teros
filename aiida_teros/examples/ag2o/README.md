# Ag2O Example - Binary Oxide

This example demonstrates how to use AiiDA-TEROS to study the surface thermodynamics of Ag2O, a binary oxide.

## Files in this Example

- `ag2o.vasp`: Bulk structure of Ag2O in VASP POSCAR format
- `ag.vasp`: Bulk structure of Ag metal in VASP POSCAR format
- `config.yaml`: Configuration file for the TEROS workflow
- `run_example.py`: Python script to run the workflow

## Tutorial: Running the Example

### Prerequisites

Before running this example, ensure you have:

1. AiiDA installed and configured with a running daemon
2. AiiDA-TEROS installed
3. VASP code configured in AiiDA

### Step 1: Modify the Configuration File

Edit the `config.yaml` file to set your specific computational parameters:

1. Set `code_label` to match your VASP code label in AiiDA:
   ```yaml
   code_label: "vasp@your_computer"
   ```

2. Adjust computation resources to match your system:
   ```yaml
   computer_options:
     resources:
       num_machines: 1
       num_cores_per_machine: 16  # Change as needed
     queue_name: "your_queue"     # Change as needed
   ```

3. Optionally modify other parameters like INCAR settings or thermodynamic parameters.

### Step 2: Run the Workflow

You have two options to run the workflow:

#### Option 1: Using the Command-Line Interface (Recommended)

```bash
cd /path/to/aiida-teros/aiida_teros/examples/ag2o
python run_example.py
```

#### Option 2: Using the Python API

```bash
cd /path/to/aiida-teros/aiida_teros/examples/ag2o
python run_example.py --api
```

### Step 3: Monitor the Workflow

After submitting the workflow, you can monitor its progress using AiiDA's CLI tools:

```bash
# List all running processes
verdi process list

# Show details of a specific process
verdi process show <PK>
```

The process PK (ID) is printed when you submit the workflow and is also saved in the `pks.txt` file.

### Step 4: Retrieve and View Results

Once the workflow completes, results will be available in:

1. AiiDA database: Structures, energies, and calculation details
2. `thermo_results/binary/` directory:
   - `surface_energy_plot.pdf`: Surface energies as a function of oxygen chemical potential
   - `termination_parameters_table.tex`: LaTeX table with surface energy parameters

## Understanding Ag2O Surface Thermodynamics

For binary oxides like Ag2O, the workflow:

1. Generates symmetric slabs from the bulk structure for the specified Miller indices `[1, 1, 0]`
2. Relaxes all terminations
3. Calculates surface Gibbs free energies as a function of oxygen chemical potential
4. Identifies the most stable surface terminations under different conditions
5. Generates plots and tables summarizing the results

The surface Gibbs free energy (γ) is calculated as:

γ = (Eslab - Nbulk·Ebulk/Nunit + ΔNO·μO) / 2A

Where:
- Eslab: Total energy of the slab
- Nbulk: Number of bulk units in the slab
- Ebulk: Energy per bulk unit
- Nunit: Number of atoms in the bulk unit
- ΔNO: Excess oxygen atoms
- μO: Oxygen chemical potential
- A: Surface area

## Additional Information

For more details, refer to the main AiiDA-TEROS documentation.