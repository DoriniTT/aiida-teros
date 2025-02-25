#!/usr/bin/env python
import os, pint, math, shutil, numpy as np, subprocess as sb
from aiida.engine import WorkChain, ToContext, append_
from aiida.orm import (
    Code, StructureData, KpointsData, Float, Dict, List, Str, Bool, SinglefileData, Int, load_node
)
from aiida.plugins import WorkflowFactory
from aiida.engine import if_, while_, return_
from aiida.engine import calcfunction
from aiida.common.extendeddicts import AttributeDict
from ase.io import read
from ase.build import molecule
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.outputs import Outcar
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from tabulate import tabulate
from aiida import load_profile
load_profile()
plt.rc('text', usetex=False)
plt.rc('font', family='sans-serif')
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
ureg = pint.UnitRegistry()
ureg.setup_matplotlib(True)

# Define workflows
VaspWorkflow = WorkflowFactory('vasp.vasp')

class AiiDATEROSWorkChain(WorkChain):
    """
    Master WorkChain for generating slab structures from a bulk structure
    and performing relaxation calculations on them.
    """

    @classmethod
    def define(cls, spec):
        """
        Define the input and output specifications for the AiiDA workflow.
        This method sets up the input parameters, namespaces, and computational resources required for the workflow.
        It also outlines the steps of the workflow and defines the exit codes for error handling.
        Inputs:
            - code (Code): VASP calculation code.
            - bulk_structure (StructureData): Bulk material structure.
            - incar_parameters_bulk (Namespace): INCAR parameters for relaxation of the bulk.
            - incar_parameters_slab (Namespace): INCAR parameters for relaxation of the slabs.
            - force_cutoff (Float): Force cutoff for relaxation.
            - steps (Int): Number of steps for relaxation.
            - kpoints_precision (Float): K-points mesh density precision.
            - potential_mapping (Dict): Mapping of potentials for each element.
            - potential_family (Str): Potential family for VASP calculations.
            - parser_settings (Dict): Settings for the parser.
            - computer_options (Dict): Computational resources and submission options.
            - miller_indices (List): Miller indices for slab generation.
            - min_slab_thickness (Float): Minimum slab thickness in Angstroms.
            - vacuum (Float): Vacuum spacing in Angstroms.
            - precision_phase_diagram (Int): Precision for the phase diagram.
            - HF_bulk (Float): Heat of formation of the bulk material.
            - path_to_graphs (Str): Path to store the graphs.
        Outline:
            1. run_relax_bulk: Relax bulk structure.
            2. inspect_relax_bulk: Inspect bulk relaxation results.
            3. generate_slabs: Generate slabs from relaxed bulk.
            4. run_relax_all_slabs: Relax all generated slabs.
            5. inspect_relax_all_slabs: Inspect relaxation results.
            6. result_binary: Create a dictionary with the results (if binary).
            7. plot_gammas_binary: Plot gamma vs delta_o with temperature (if binary).
            6. result_ternary: Create a dictionary with the results (if ternary).
            7. plot_gamma_delta_o_withtemp: Plot gamma vs delta_o with temperature (if ternary).
            8. plot_phase_diagram: Plot the combined figures (if ternary).
        Outputs:
            - bulk (Namespace): Exposed outputs from VaspWorkflow.
            - relaxations (Namespace): Dynamic namespace for relaxation outputs.
            - stable_structures (Namespace): Dynamic namespace for stable structures.
        Exit Codes:
            - 201: ERROR_GENERATE_SLABS_FAILED - Slab generation failed.
            - 202: ERROR_RELAX_SLABS_FAILED - Relaxation of slabs failed.
            - 203: ERROR_RELAX_BULK_FAILED - Relaxation of bulk structure failed.
        """
        super().define(spec)

        # Define inputs
        spec.input(
            'code', 
            valid_type=Code, 
            help='VASP calculation code.'
        )
        spec.input(
            'bulk_structure', 
            valid_type=StructureData, 
            help='Bulk material structure.'
        )

        spec.input(
            'bulk_metal', 
            valid_type=StructureData, 
            help='Bulk material of the metal.'
        )

        spec.input_namespace(
            'terminations', 
            valid_type=StructureData, 
            required=False,
            help='Specific terminations to calculate.',
        )
        
        # INCAR parameters
        spec.input_namespace(
            'incar_parameters_bulk', 
            help='INCAR parameters for relaxation of the bulk.', 
            dynamic=True
        )
        spec.input_namespace(
            'incar_parameters_bulk_metal', 
            help='INCAR parameters for relaxation of the bulk metal.', 
            dynamic=True
        )

        spec.input_namespace(
            'incar_parameters_o2',
            help='INCAR parameters for relaxation of the O2 molecule.',
            dynamic=True
        )

        spec.input_namespace(
            'incar_parameters_slab', 
            help='INCAR parameters for relaxation of the slabs.', 
            dynamic=True
        )
        
        # Relaxation settings
        spec.input(
            'force_cutoff', 
            valid_type=Float, 
            default=lambda: Float(0.01),
            help='Force cutoff for relaxation.'
        )
        spec.input(
            'steps', 
            valid_type=Int, 
            default=lambda: Int(1000), 
            help='Number of steps for relaxation.'
        )
        spec.input(
            'kpoints_precision', 
            valid_type=Float, 
            default=lambda: Float(0.3), 
            help='K-points mesh density precision.'
        )
        
        # Potential settings
        spec.input(
            'potential_mapping', 
            valid_type=Dict, 
            help='Mapping of potentials for each element.'
        )
        spec.input(
            'potential_family', 
            valid_type=Str, 
            help='Potential family for VASP calculations.'
        )
        
        # Parser settings
        spec.input(
            'parser_settings', 
            valid_type=Dict, 
            help='Settings for the parser.'
        )
        
        # Computational resources
        spec.input(
            'computer_options', 
            valid_type=Dict, 
            help='Computational resources and submission options.'
        )
        
        # Slab generation settings
        spec.input(
            'miller_indices', 
            valid_type=List, 
            help='Miller indices for slab generation.'
        )
        spec.input(
            'min_slab_thickness', 
            valid_type=Float, 
            help='Minimum slab thickness in Angstroms.'
        )
        spec.input(
            'vacuum', 
            valid_type=Float, 
            default=lambda: Float(15), 
            help='Vacuum spacing in Angstroms.'
        )
        
        # Phase diagram settings
        spec.input(
            'precision_phase_diagram', 
            valid_type=Int, 
            default=lambda: Int(500), 
            help='Precision for the phase diagram.'
        )
        # Formation enthalpy
        spec.input(
            'HF_bulk', 
            valid_type=Float, 
            help='Heat of formation of the bulk material.'
        )
        # Path for storing graphs
        spec.input(
            'path_to_graphs', 
            valid_type=Str, 
            default=lambda: Str(os.getcwd()), 
            help='Path to store the graphs.'
        )

        # Define the workflow outline
        spec.outline(
            cls.run_relax_bulk,            # Step 1: Relax bulk structure
            cls.inspect_relax_bulk,        # Step 2: Inspect bulk relaxation results
            cls.run_relax_bulk_metal,            # Step 1: Relax bulk structure
            cls.inspect_relax_bulk_metal,        # Step 2: Inspect bulk relaxation results
            cls.run_relax_o2,
            cls.inspect_relax_o2,
            cls.generate_slabs,            # Step 3: Generate slabs from relaxed bulk
            cls.run_relax_all_slabs,       # Step 4: Relax all generated slabs
            cls.inspect_relax_all_slabs,   # Step 5: Inspect relaxation results

            if_(cls.is_binary)(
            cls.result_binary,         # Step 6: Create a dictionary with the results
            cls.plot_gammas_binary,    # Step 7: Plot gamma vs delta_o with temperature
            ).else_(
            cls.result_ternary,        # Step 6: Create a dictionary with the results
            cls.plot_gammas_ternary,  # Step 7: Plot gamma vs delta_o with temperature
            cls.plot_phase_diagram,    # Step 8: Plot the combined figures
            )
        )

        # Define outputs
        spec.expose_outputs(
            VaspWorkflow, 
            namespace='bulk'
        )
        spec.expose_outputs(
            VaspWorkflow, 
            namespace='bulk_metal1',
        )

        spec.expose_outputs(
            VaspWorkflow,
            namespace='o2',
        )
        
        spec.output_namespace(
            'relaxations', 
            dynamic=True
        )
        spec.output_namespace(
            'stable_structures', 
            dynamic=True
        )

        # Define exit codes for error handling
        spec.exit_code(
            201, 
            'ERROR_GENERATE_SLABS_FAILED', 
            message='Slab generation failed.'
        )
        spec.exit_code(
            202, 
            'ERROR_RELAX_SLABS_FAILED', 
            message='Relaxation of slabs failed.'
        )
        spec.exit_code(
            203, 
            'ERROR_RELAX_BULK_FAILED', 
            message='Relaxation of bulk structure failed.'
        )

    def is_binary(self):

        #* Check if the bulk structure is binary or ternary
        num_atoms = len(set(self.inputs.bulk_structure.get_ase().get_chemical_symbols()))

        if num_atoms == 2:
            return True
        elif num_atoms == 3:
            return False
        else:
            raise ValueError("Bulk structure must contain two or three elements.")

    def run_relax_bulk(self):
        """
        Run a VASP relaxation of the bulk structure.
        """
        self.report('Running relaxation of the bulk structure.')

        try:
            # Get INCAR parameters for bulk relaxation
            incar_bulk = self.inputs.incar_parameters_bulk

            # Get the VASP builder with bulk INCAR parameters
            builder = self.get_vasp_builder(
                structure=self.inputs.bulk_structure,
                incar_parameters=incar_bulk,
                kpoint_density=self.inputs.kpoints_precision.value,
                label='relax_bulk_structure',
                description='Bulk relaxation of the material'
            )

            # Submit the bulk relaxation calculation
            future = self.submit(builder)
            #future = load_node(139563)
            self.report(f'Submitted VASP relaxation for bulk structure with PK {future.pk}')
            return ToContext(relax_bulk=future)

        except Exception as e:
            self.report(f'Failed to submit bulk relaxation: {e}')
            return self.exit_codes.ERROR_RELAX_BULK_FAILED

    def inspect_relax_bulk(self):
        """
        Inspect the result of the bulk relaxation.
        """
        if not self.ctx.relax_bulk.is_finished_ok:
            self.report(f'Relaxation of bulk failed with exit status {self.ctx.relax_bulk.exit_status}')
            return self.exit_codes.ERROR_RELAX_BULK_FAILED

        # Expose bulk relaxation outputs
        self.out_many(self.exposed_outputs(self.ctx.relax_bulk, VaspWorkflow, namespace='bulk'))
        self.report('Bulk relaxation completed successfully.')

    def run_relax_bulk_metal(self):
        """
        Run a VASP relaxation of the bulk structure.
        """
        self.report('Running relaxation of the bulk structure.')

        try:
            # Get INCAR parameters for bulk relaxation
            incar_bulk = self.inputs.incar_parameters_bulk_metal

            # Get the VASP builder with bulk INCAR parameters
            builder = self.get_vasp_builder(
                structure=self.inputs.bulk_metal,
                incar_parameters=incar_bulk,
                kpoint_density=self.inputs.kpoints_precision.value,
                label='relax_bulk_metal_structure',
                description='Bulk relaxation of the metal'
            )

            # Submit the bulk relaxation calculation
            future = self.submit(builder)
            #future = load_node(139563)
            self.report(f'Submitted VASP relaxation for bulk structure with PK {future.pk}')
            return ToContext(relax_bulk_metal=future)

        except Exception as e:
            self.report(f'Failed to submit bulk metal relaxation: {e}')
            return self.exit_codes.ERROR_RELAX_BULK_FAILED

    def inspect_relax_bulk_metal(self):
        """
        Inspect the result of the bulk relaxation.
        """
        if not self.ctx.relax_bulk_metal.is_finished_ok:
            self.report(f'Relaxation of bulk metal failed with exit status {self.ctx.relax_bulk_metal.exit_status}')
            return self.exit_codes.ERROR_RELAX_BULK_FAILED

        # Expose bulk relaxation outputs
        self.out_many(self.exposed_outputs(self.ctx.relax_bulk_metal, VaspWorkflow, namespace='bulk_metal1'))
        self.report('Metal bulk relaxation completed successfully.')

    def run_relax_o2(self):
        """
        Create an O2 molecule with ASE in a box and run a VASP relaxation.
        """
        self.report('Running relaxation of the O2 molecule.')

        try:
            # Create an O2 molecule using ASE
            o2_molecule = molecule('O2')
            o2_molecule.set_cell([13, 14, 15])
            o2_molecule.center()  # Center the molecule in the cell
            structure = StructureData(ase=o2_molecule)

            # Define default INCAR parameters for O2 relaxation
            incar_o2 = self.inputs.incar_parameters_o2

            # Get the VASP builder with O2 INCAR parameters
            builder = self.get_vasp_builder(
                structure=structure,
                incar_parameters=incar_o2,
                kpoint_density=self.inputs.kpoints_precision.value,
                label='relax_o2_molecule',
                description='Relaxation of the O2 molecule',
                slab=True #we treat the molecule as a slab
            )

            # Submit the O2 relaxation calculation
            future = self.submit(builder)
            self.report(f'Submitted VASP relaxation for O2 molecule with PK {future.pk}')
            return ToContext(relax_o2=future)

        except Exception as e:
            self.report(f'Failed to submit O2 relaxation: {e}')
            return self.exit_codes.ERROR_RELAX_O2_FAILED

    def inspect_relax_o2(self):
        """
        Inspect the result of the O2 relaxation.
        """
        if not self.ctx.relax_o2.is_finished_ok:
            self.report(f'Relaxation of O2 failed with exit status {self.ctx.relax_o2.exit_status}')
            return self.exit_codes.ERROR_RELAX_O2_FAILED

        # Expose O2 relaxation outputs
        self.out_many(self.exposed_outputs(self.ctx.relax_o2, VaspWorkflow, namespace='o2'))
        self.report('O2 relaxation completed successfully.')

    def generate_slabs(self):
        """
        Generate slab structures from the relaxed bulk structure using pymatgen's SlabGenerator.
        """

        try:

            if 'terminations' in self.inputs:  # Check if 'terminations' is provided as input
                self.ctx.slabs = list(self.inputs.terminations.values())
                self.report('The user requested specific terminations.')

            else:
                self.report('Generating slab structures from relaxed bulk structure.')
                # Retrieve the relaxed bulk structure
                relax_bulk_calc = self.ctx.relax_bulk
                if not relax_bulk_calc.is_finished_ok:
                    self.report('Bulk relaxation did not finish successfully.')
                    return self.exit_codes.ERROR_RELAX_BULK_FAILED

                relaxed_bulk = relax_bulk_calc.outputs.structure

                # Convert AiiDA StructureData to pymatgen Structure
                bulk_structure = relaxed_bulk.get_pymatgen()
                miller = tuple(self.inputs.miller_indices.get_list())  # e.g., [1, 1, 1]
                min_slab_thickness = self.inputs.min_slab_thickness.value  # in Angstroms
                vacuum = self.inputs.vacuum.value  # in Angstroms

                slab_generator = SlabGenerator(
                    bulk_structure,
                    miller,
                    min_slab_thickness,
                    vacuum,
                    lll_reduce=True,
                    center_slab=True
                )
                slabs = slab_generator.get_slabs(symmetrize=True)
                all_structures = []
                for slab in slabs:
                    slab = slab.get_orthogonal_c_slab()
                    ase_atoms = AseAtomsAdaptor().get_atoms(slab)
                    all_structures.append(StructureData(ase=ase_atoms))

                if not all_structures:
                    self.report('No slabs were generated.')
                    return self.exit_codes.ERROR_GENERATE_SLABS_FAILED

                self.ctx.slabs = all_structures
                self.report(f'Generated {len(all_structures)} slab structures.')

            # Assign slabs to dynamic output namespace
            structure_dict = {}
            for n, slab in enumerate(self.ctx.slabs, start=1):
                structure_dict[f'ST_{n}'] = slab
            self.out('relaxations.generated_terminations', structure_dict)

        except Exception as e:
            self.report(f'Failed to generate slabs: {e}')
            return self.exit_codes.ERROR_GENERATE_SLABS_FAILED

    def run_relax_all_slabs(self):
        """
        Run VASP relaxation calculations on all generated slabs.
        """
        self.report('Submitting relaxation calculations for all slabs.')

        self.ctx.relax_calculations = []

        for n, slab in enumerate(self.ctx.slabs, start=1):
            self.report(f'Running relaxation for slab structure_{n}.')

            try:
                # Get INCAR parameters for slab relaxation
                incar_slabs = self.inputs.incar_parameters_slab

                # Get the VASP builder with slab INCAR parameters
                builder = self.get_vasp_builder(
                    structure=slab,
                    incar_parameters=incar_slabs,
                    slab=True,
                    kpoint_density=self.inputs.kpoints_precision.value,
                    label=f'relax_slab_{n}',
                    description=f'Relaxation of slab structure_{n}'
                )

                # Submit the relaxation calculation
                future = self.submit(builder)
                self.report(f'Submitted relaxation for slab_{n} with PK {future.pk}')
                self.to_context(relax_calculations=append_(future))

            except Exception as e:
                self.report(f'Failed to submit relaxation for slab_{n}: {e}')
                return self.exit_codes.ERROR_RELAX_SLABS_FAILED

    def inspect_relax_all_slabs(self):
        """
        Inspect the results of all relaxation calculations.
        """
        self.report('Inspecting relaxation calculations for all slabs.')

        self.ctx.relaxed_structures = []
        self.ctx.relax_calcs = []

        for n, calc in enumerate(self.ctx.relax_calculations, start=1):
            calc = calc.called_descendants[-1]
            if not calc.is_finished_ok:
                self.report(f'Relaxation calculation for slab_{n} failed with exit status {calc.exit_status}')
                return self.exit_codes.ERROR_RELAX_SLABS_FAILED

            # Store the relaxed structure and calculation
            relaxed_structure = calc.outputs.structure
            self.ctx.relaxed_structures.append(relaxed_structure)
            self.ctx.relax_calcs.append(calc)
            
            # Aggregate outputs into a dictionary
            structure_dict = {}
            for output in calc.base.links.get_outgoing().all():
                structure_dict[output.link_label] = output.node
            # Assign the Dict node to the dynamic output namespace
            self.out(f'relaxations.ST_{n}', structure_dict)

            self.report(f'Relaxation for slab_{n} completed successfully.')

    def result_binary(self):
        """
        Calculate the surface total energy for all slab terminations and store the results.

        Parameters:
        workchain_ctx: AiiDA WorkChain context containing VaspCalculation nodes.
        bulk_structure (StructureData): AiiDA StructureData of the bulk.
        E_bulk (float): DFT total energy of the bulk.

        Returns:
        dict: Dictionary containing lower and upper surface Gibbs free energy values for each termination.
        """

        self.report('Calculating the surface total energy for a binary system.')

        self.report('Starting calculation of surface total energy for a binary system.')

        # Create the file thermo_results/binary if it does not exist
        os.makedirs(f'{self.inputs.path_to_graphs.value}/thermo_results/binary', exist_ok=True)

        delta_Hf = self.inputs.HF_bulk.value
        E_bulk = self.ctx.relax_bulk.outputs.misc.get_dict()['total_energies']['energy_extrapolated']
        bulk_structure = self.inputs.bulk_structure.get_ase()
        bulk_metal_structure = self.inputs.bulk_metal.get_ase()
        # Get the number of atoms in bulk_metal_structure as an integer
        num_atoms_bulk_metal_structure = len(bulk_metal_structure)
        E_bulk_metal = self.ctx.relax_bulk_metal.outputs.misc.get_dict()['total_energies']['energy_extrapolated'] / num_atoms_bulk_metal_structure

        element_counts = Counter(atom.symbol for atom in bulk_structure)
        gcd_value = np.gcd.reduce(list(element_counts.values()))
        
        for element, natoms in element_counts.items():
            if element == 'O':
                y = natoms / gcd_value
            else:
                x = natoms / gcd_value
        self.report('Determined the minimal composition of the bulk structure.')

        E_bulk_per_fu = E_bulk / gcd_value  # Energy per formula unit

        self.report('Calculated limits for the chemical potential of oxygen.')
        lower_limit = 1/y * (E_bulk_per_fu - x * E_bulk_metal)
        upper_limit = lower_limit + 1/y * delta_Hf
        #Report the all variables used in the calculation
        self.report('The total energy of the bulk structure is: {}'.format(E_bulk_per_fu))
        self.report('The total energy of the bulk metal structure is: {}'.format(E_bulk_metal))
        self.report('The heat of formation of the bulk structure is: {}'.format(delta_Hf))
        self.report('The minimal composition of the bulk structure is: {}'.format(gcd_value))
        self.report('The number of atoms of the first element is: {}'.format(x))
        self.report('The number of atoms of the second element is: {}'.format(y))
        self.report('The lower limit for the chemical potential of oxygen is: {}'.format(lower_limit))
        self.report('The upper limit for the chemical potential of oxygen is: {}'.format(upper_limit))

        mu_O_values = self.ctx.mu_O_values = [lower_limit, upper_limit]
        gammas = {}
        termination_data = []  # Data collection for LaTeX table
        for i, calc in enumerate(self.ctx.relax_calcs):
            self.report(f'Processing relaxation calculation {i+1}.')

            slab_structure = calc.outputs.structure
            slab_structure_ase = slab_structure.get_ase()
            E_slab = calc.outputs.misc.get_dict()['total_energies']['energy_extrapolated']

            element_symbols = list(set(atom.symbol for atom in slab_structure_ase if atom.symbol != 'O'))
            if len(element_symbols) != 1:
                raise ValueError("Structure must contain exactly one non-oxygen element.")
            element = element_symbols[0]

            N_element_slab = len([atom for atom in slab_structure_ase if atom.symbol == element])
            N_O_slab = len([atom for atom in slab_structure_ase if atom.symbol == 'O'])

            if N_element_slab == 0 or N_O_slab == 0:
                raise ValueError("Number of non-oxygen or oxygen atoms cannot be zero.")

            lengths = slab_structure_ase.get_cell_lengths_and_angles()
            if lengths[0] <= 0 or lengths[1] <= 0:
                raise ValueError("Cell dimensions must be positive values.")

            A = lengths[0] * lengths[1]  # Area is the product of the x and y coordinates of the cell

            if A == 0:
                raise ValueError("Surface area cannot be zero.")

            gamma_values = []
            for mu_O in mu_O_values:
                gamma = (2*A)**(-1) * (E_slab - (N_element_slab/x) * E_bulk_per_fu + ((y/x) * N_element_slab - N_O_slab) * mu_O)
                gamma = gamma * 1.602176634e-19 * 1e20 # eV/Å² to J/m²
                gamma_values.append(gamma)

            self.report(f'Calculated surface Gibbs free energy for termination {i+1}.')

            gammas[f'Termination {i+1}'] = {
            'gamma_lower': gamma_values[0],
            'gamma_upper': gamma_values[1]
            }

            # Collect data for the LaTeX table
            termination_data.append([
                f'Termination {i+1}',
                E_slab,
                N_element_slab,
                N_O_slab,
                A
            ])

        self.report(f'Calculated surface Gibbs free energy for termination {i+1}.')

        self.ctx.gammas_binary = gammas
        self.report('Completed calculation of surface total energy for all terminations.')

        # Define headers for the LaTeX table
        headers = [
            "Termination",
            "$E_{\\mathrm{slab}}$ (eV)",
            "$N_{\\mathrm{metal}}$",
            "$N_\\mathrm{O}$",
        ] 
        
        # Generate the LaTeX table using tabulate
        table_latex = tabulate(termination_data, headers=headers, tablefmt="latex_raw")
        
        # Save the LaTeX table to a file
        os.makedirs(f'{self.inputs.path_to_graphs.value}/thermo_results/binary', exist_ok=True)
        with open(f'{self.inputs.path_to_graphs.value}/thermo_results/binary/termination_parameters_table.tex', 'w') as f:
            f.write("\\documentclass{article}\n")
            f.write("\\usepackage{amsmath}\n")
            f.write("\\usepackage{geometry}\n")
            f.write("\\geometry{a4paper, margin=1in}\n")
            f.write("\\begin{document}\n")
            f.write("\\section*{Termination Parameters}\n")
            f.write(table_latex)
            f.write("\n")
            f.write("\\end{document}\n")

        self.report('Generated LaTeX table with termination parameters.')

    def plot_gammas_binary(self):
        """
        Plot surface energies as a function of oxygen chemical potential for different terminations.
        All data is defined within the function and saves the plot automatically.
        """
        # Define the data inside the function
        gammas = self.ctx.gammas_binary
        mu_O_values = self.ctx.mu_O_values

        # Plot settings
        figsize = (10, 6)
        line_width = 2.5
        marker_size = 8

        # Create figure and axis
        fig, ax = plt.subplots(figsize=figsize)

        # Create color cycle for different terminations
        colors = plt.cm.get_cmap('tab10', len(gammas))

        # Plot lines for each termination
        for n, ((term_name, term_data), color) in enumerate(zip(gammas.items(), colors.colors), start=1):
            # Create points for the line
            x_points = mu_O_values
            y_points = [term_data['gamma_lower'], term_data['gamma_upper']]
            # Plot the line with markers
            ax.plot(x_points, y_points, '-', label=f'ST-A{str(n)}', color=color, linewidth=line_width)
            ax.plot(x_points, y_points, 'o', color=color, markersize=marker_size, markeredgecolor='black', alpha=0.8)

        # Customize the plot
        ax.set_xlabel(r'$\mu_O$ (eV)', fontsize=18)
        ax.set_ylabel(r'$\gamma$ (J/m²)', fontsize=18)

        # Ticks settings for better readability
        ax.tick_params(axis='both', which='major', labelsize=16, length=6, width=1.5)
        ax.tick_params(axis='both', which='minor', length=4, width=1.2)

        # Add grid for easier data interpretation
        ax.grid(True, linestyle='--', alpha=0.6, linewidth=1)

        # Add legend outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=16, title='Terminations', title_fontsize=16)

        # Adjust layout to prevent cutting off the legend and labels
        plt.tight_layout()

        # Save the figure
        plt.savefig(f'{self.inputs.path_to_graphs.value}/thermo_results/binary/surface_energy_plot.pdf', 
                    bbox_inches='tight', dpi=200)

    def result_ternary(self):

        os.makedirs(f'{self.inputs.path_to_graphs.value}/thermo_results/ternary', exist_ok=True)

        self.ctx.dict_results = {}
        termination_data = []  # Data collection for LaTeX table
        for n, calc in enumerate(self.ctx.relax_calcs, start=1):

            #energy_ag = -2.8289*ureg.eV
            #energy_ag = self.inputs.total_energy_first_element.value * ureg.eV
            #energy_o2 = self.inputs.total_energy_o2.value * ureg.eV
            # Count the number of atoms in the bulk metal structure
            bulk_metal_structure = self.inputs.bulk_metal.get_ase()
            num_atoms_bulk_metal_structure = len(bulk_metal_structure)
            energy_ag = (self.ctx.relax_bulk_metal.outputs.misc.get_dict()['total_energies']['energy_extrapolated'] / num_atoms_bulk_metal_structure) * ureg.eV
            energy_o2 = (self.ctx.relax_o2.outputs.misc.get_dict()['total_energies']['energy_extrapolated'] / 2) * ureg.eV
            #energy_o2 = -9.82*ureg.eV
            logP_P0 = 1
            kB = 1.38064852e-23*ureg["J/K"] # Boltzmann's constant
            mu_ex = 0.4*ureg.eV
            delta_mu_ag = 0*ureg.eV

            def write_FolderData(folderdata, *args):
                for arg in args:
                    with folderdata.open(arg, mode='rb') as source:
                        with open(arg, mode='wb') as destination:
                            shutil.copyfileobj(source, destination)

            label = calc.label
            self.ctx.exemple_label = label
            self.ctx.dict_results[label] = {}

            write_FolderData(calc.outputs.retrieved, 'vasprun.xml')
            outcar_ase = read('vasprun.xml', format='vasp-xml')
            os.remove('vasprun.xml')

            write_FolderData(self.ctx.relax_bulk.outputs.retrieved, 'vasprun.xml')
            outcar_ase_bulk = read('vasprun.xml', format='vasp-xml')
            divide_minimal_comp = self.check_minimal_composition(StructureData(ase=outcar_ase_bulk))
            elements = outcar_ase_bulk.get_chemical_symbols()
            self.ctx.unique_symbols = unique_symbols = list(dict.fromkeys(elements))
            nAg_bulk = outcar_ase_bulk.get_chemical_symbols().count(unique_symbols[0])/divide_minimal_comp
            nMe_bulk = outcar_ase_bulk.get_chemical_symbols().count(unique_symbols[1])/divide_minimal_comp
            nO_bulk = outcar_ase_bulk.get_chemical_symbols().count(unique_symbols[2])/divide_minimal_comp
            stoichiometry = [nAg_bulk, nMe_bulk, nO_bulk]
            energy_bulk = outcar_ase_bulk.get_potential_energy()*ureg["eV"]/divide_minimal_comp
            composition_bulk = [nAg_bulk, nMe_bulk, nO_bulk]
            os.remove('vasprun.xml')

            self.ctx.dict_results[label]['energy_slab'] = outcar_ase.get_potential_energy()*ureg["eV"]
            self.ctx.dict_results[label]['a'] = outcar_ase.cell.cellpar()[0]*ureg["angstrom"]
            self.ctx.dict_results[label]['b'] = outcar_ase.cell.cellpar()[1]*ureg["angstrom"]
            self.ctx.dict_results[label]['elements'] = {}
            self.ctx.dict_results[label]['elements'][unique_symbols[0]] = outcar_ase.get_chemical_symbols().count(unique_symbols[0])
            self.ctx.dict_results[label]['elements'][unique_symbols[1]] = outcar_ase.get_chemical_symbols().count(unique_symbols[1])
            self.ctx.dict_results[label]['elements'][unique_symbols[2]] = outcar_ase.get_chemical_symbols().count(unique_symbols[2])

            self.ctx.dict_results[label]['Delta_Me_Ag'] = composition_bulk[1] * self.ctx.dict_results[label]['elements'][unique_symbols[0]] - composition_bulk[0] * self.ctx.dict_results[label]['elements'][unique_symbols[1]]
            self.ctx.dict_results[label]['Delta_Me_O'] = composition_bulk[1] * self.ctx.dict_results[label]['elements'][unique_symbols[2]] - composition_bulk[2] * self.ctx.dict_results[label]['elements'][unique_symbols[1]]

            self.ctx.dict_results[label]['psi'] = self.ctx.dict_results[label]['energy_slab'] - self.ctx.dict_results[label]['elements'][unique_symbols[1]] * energy_bulk / nMe_bulk - energy_ag * self.ctx.dict_results[label]['Delta_Me_Ag'] - 1/2 * self.ctx.dict_results[label]['Delta_Me_O'] * energy_o2

            range_of_T =  [t*ureg.K for t in range(100, 1300, 100)] # Temperature (K)
            Delta_Gs = [-0.15, -0.341, -0.548, -0.765, -0.991, -1.222, -1.460, -1.702, -1.949, -2.199, -2.453, -2.711]
            Delta_Gs = [g*ureg.eV for g in Delta_Gs]
            #delta_mu = [1/2 * (g + kB * T * logP_P0) + mu_ex for T, g in zip(range_of_T, Delta_Gs)]
            lim_delta_o = self.inputs.HF_bulk.value / stoichiometry[2]
            delta_mu = np.linspace(0, lim_delta_o, 20)
            delta_mu = [d*ureg.eV for d in delta_mu]
            gamma = []
            delta_mu_copy = delta_mu.copy()
            #delta_mu_copy =  [mu.magnitude for mu in delta_mu_copy] # Temperature (K)
            self.ctx.dict_results[label]['delta_mu'] = delta_mu_copy

            for m in delta_mu:
                g = (1 / ( 2 * self.ctx.dict_results[label]['a'] * self.ctx.dict_results[label]['b']) ) * (self.ctx.dict_results[label]['psi'] - self.ctx.dict_results[label]['Delta_Me_Ag'] * delta_mu_ag - self.ctx.dict_results[label]['Delta_Me_O'] * m)
                gamma.append(g.to('J/m^2').magnitude)
            self.ctx.dict_results[label]['gamma_delta_o'] = gamma

            # Collect data for the LaTeX table
            termination_data.append([
                label,
                self.ctx.dict_results[label]['Delta_Me_O'],
                self.ctx.dict_results[label]['Delta_Me_Ag'],
                self.ctx.dict_results[label]['psi']
            ])

            # Define headers for the LaTeX table
            headers = [
                "Termination Label",
                "$N_O - xN_B$",
                "$N_A - yN_B$",
                "$\Theta$ (eV)"
            ]

            # Generate the LaTeX table using tabulate
            table_latex = tabulate(termination_data, headers=headers, tablefmt="latex")

            # Save the LaTeX table to a file
            with open(f'{self.inputs.path_to_graphs.value}/thermo_results/ternary/ternary_parameters_table.tex', 'w') as f:
                f.write(r"\documentclass{article}\n")
                f.write(r"\usepackage{amsmath}\n")
                f.write(r"\usepackage{geometry}\n")
                f.write(r"\geometry{a4paper, margin=1in}\n")
                f.write(r"\begin{document}\n")
                f.write(r"\section*{Termination Parameters Table}\n")
                f.write(table_latex)
                f.write(r"\end{document}\n")

            self.report('Generated LaTeX table with ternary termination parameters and saved to file.')

            precision = self.inputs.precision_phase_diagram.value
            #* Calculating the range of Delta_Ag and Delta_O
            self.ctx.lim_delta_ag = lim_delta_ag = self.inputs.HF_bulk.value / stoichiometry[0]
            self.ctx.lim_delta_o = lim_delta_o = self.inputs.HF_bulk.value / stoichiometry[2]
            delta_ag = np.linspace(0, lim_delta_ag, precision)
            delta_o = np.linspace(0, lim_delta_o, precision)
            gamma_delta_ag_delta_o = np.zeros((precision, precision))
            for i, ag in enumerate(delta_ag):
                for j, o in enumerate(delta_o):
                    gamma_delta_ag_delta_o[i, j] = (1 / ( 2 * self.ctx.dict_results[label]['a'].magnitude * self.ctx.dict_results[label]['b'].magnitude) ) * (self.ctx.dict_results[label]['psi'].magnitude - self.ctx.dict_results[label]['Delta_Me_Ag'] * ag - self.ctx.dict_results[label]['Delta_Me_O'] * o)

            self.ctx.dict_results[label]['delta_ag'] = delta_ag
            self.ctx.dict_results[label]['delta_o'] = delta_o
            self.ctx.dict_results[label]['gamma_delta_ag_delta_o'] = gamma_delta_ag_delta_o

            for key, value in self.ctx.dict_results[label].items():
                if isinstance(value, ureg.Quantity):
                    self.ctx.dict_results[label][key] = value.magnitude

        self.ctx.dict_results = dict(sorted(self.ctx.dict_results.items(), key=lambda item: item[0]))

        #* Getting the most stable structures!
        delta_mu = np.array([dm.to('eV').magnitude for dm in self.ctx.dict_results[self.ctx.exemple_label]['delta_mu']])
        # Define target delta_mu values
        target_delta_mu = [0.0, -2.0]
        target_labels = ['delta_mu = 0 eV', 'delta_mu = -2 eV']
        most_stable = {}

        for target, label in zip(target_delta_mu, target_labels):
            # Find the index of the closest delta_mu value
            idx = np.argmin(np.abs(delta_mu - target))
            closest_delta_mu = delta_mu[idx]
            print(f"Closest delta_mu to {target} eV is {closest_delta_mu} eV at index {idx}")

            # Extract gamma_delta_o for all structures at this index
            gamma_at_target = {}
            for structure_name, info_structure in self.ctx.dict_results.items():
                gamma_value = info_structure['gamma_delta_o'][idx]
                gamma_at_target[structure_name] = gamma_value

            # Find the structure with the minimum gamma_delta_o
            most_stable_structure = min(gamma_at_target, key=gamma_at_target.get)
            most_stable[target] = most_stable_structure

            # Report the most stable structure
            self.report(f'Most stable structure at {label}: {most_stable_structure}')

        # Store the most stable structures in the context
        self.ctx.most_stable_delta_mu_0 = most_stable.get(0.0, None)
        self.ctx.most_stable_delta_mu_minus_2 = most_stable.get(-2.0, None)

        stable_structures = []
        for n, calc in enumerate(self.ctx.relax_calcs, start=1):

            if calc.label == self.ctx.most_stable_delta_mu_0:
                stable_structures.append(calc)
            elif calc.label == self.ctx.most_stable_delta_mu_minus_2:
                stable_structures.append(calc)

        for n, calc in enumerate(stable_structures):
            structure_dict = {}
            for output in calc.base.links.get_outgoing().all():
                structure_dict[output.link_label] = output.node
            # Assign the Dict node to the dynamic output namespace
            if n == 0:
                self.out(f'stable_structures.{calc.label}', structure_dict)
            elif n == 1:
                self.out(f'stable_structures.{calc.label}', structure_dict)

    def plot_gammas_ternary(self):
    
        import matplotlib.pyplot as plt
        import seaborn as sns
        import numpy as np
        import subprocess as sb

        colors = sns.color_palette("husl", min(max(len(self.ctx.dict_results), 8), 20))  # 8 distinct colors
        fig, ax1 = plt.subplots(figsize=(8, 6), dpi=200)

        # Convert delta_mu to a sorted numpy array for efficient processing
        delta_mu = np.array([dm.to('eV').magnitude for dm in self.ctx.dict_results[self.ctx.exemple_label]['delta_mu']])

        labels = []
        lines = []
        min_values = {}
        max_values = {}

        for n, (structure_name, info_structure) in enumerate(self.ctx.dict_results.items(), start=1):
            # Print the most stable structures
            print(f'{structure_name} possess a gamma of {str(info_structure["gamma_delta_o"][0])} at delta_mu_o = 0')

            # Plot the lines
            line, = ax1.plot(delta_mu, info_structure['gamma_delta_o'], label=f'ST-B{str(n)}', linewidth=2.5, color=colors[(n - 1) % len(colors)])
            labels.append(f'ST-B{str(n)}')
            lines.append(line)

        # Set axis limits based on min and max values
        all_gamma_values = [value for info_structure in self.ctx.dict_results.values() for value in info_structure['gamma_delta_o']]
        y_min = min(all_gamma_values) - 0.2
        y_max = max(all_gamma_values) + 0.2

        fontsize_ticks = 18
        ax1.set_xlabel(r'$\Delta \mu_{\mathrm{O}}$ (eV)', fontsize=fontsize_ticks + 2)
        ax1.set_ylabel(r'$\gamma$ (J/m$^2$)', fontsize=fontsize_ticks + 2)
        primary_xticks = np.linspace(-2, 0, 5)
        ax1.set_xticks(primary_xticks)
        ax1.set_xticklabels(primary_xticks, fontsize=fontsize_ticks - 2)
        yticks = np.round(np.linspace(y_min, y_max, 8), 0)
        ax1.set_yticks(yticks)
        ax1.set_yticklabels([int(tick) for tick in yticks], fontsize=fontsize_ticks - 2)
        ax1.set_xlim(-2.0, 0)  # Set the limits of the x-axis
        ax1.set_ylim(y_min, y_max)  # Set the limits of the y-axis

        # Add first secondary x-axis (e.g., Temperature)
        ax2 = ax1.twiny()
        ax2.spines['top'].set_position(('outward', 10))
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(np.linspace(-4, 0, 9))
        T = [1100, 1000, 900, 800, 700, 600, 500, 400, 300]
        ax2.set_xticklabels(T, fontsize=fontsize_ticks - 5)
        ax2.set_xlabel('Temperature (K) @ 1 atm', fontsize=fontsize_ticks - 5)

        # Add second secondary x-axis (e.g., Pressure)
        ax3 = ax1.twiny()
        ax3.spines['top'].set_position(('outward', 60))
        ax3.set_xlim(ax1.get_xlim())
        ax3.set_xticks(np.linspace(-4, 0, 12))
        ax3.set_xticklabels([None, None, r'$10^{-27}$', r'$10^{-24}$', r'$10^{-21}$', r'$10^{-18}$', 
                             r'$10^{-15}$', r'$10^{-12}$', r'$10^{-9}$', r'$10^{-6}$', 
                             r'$10^{-3}$', r'$1$'], fontsize=fontsize_ticks - 5)
        ax3.set_xlabel('Pressure (atm) @ 300 K', fontsize=fontsize_ticks - 5)

        # Create the legend outside of the plot area
        ax1.legend(lines, labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13, borderaxespad=0.)

        sb.run([f'mkdir -p {self.inputs.path_to_graphs.value}/thermo_results'], shell=True)
        plt.tight_layout()
        plt.savefig(f"{self.inputs.path_to_graphs.value}/thermo_results/surface_free_energies.pdf", bbox_inches="tight")

        # Store the minimum and maximum values found in the context
        self.ctx.min_values = min_values
        self.ctx.max_values = max_values

    def plot_phase_diagram(self):

        import matplotlib.pyplot as plt
        import numpy as np
        import seaborn as sns
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.patches as mpatches
        import pint  # Ensure you have pint for unit handling

        ureg = pint.UnitRegistry()  # Initialize unit registry

        print_most_stable_structure = True
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(13, 6), dpi=200,
            gridspec_kw={'wspace': 0, 'width_ratios': [2, 1]}
        )  # Adjust the width_ratios

        labels = sorted(self.ctx.dict_results.keys())
        delta_ag = self.ctx.dict_results[labels[0]]['delta_ag']
        delta_o = self.ctx.dict_results[labels[0]]['delta_o']
        gamma_min = np.full((len(delta_ag), len(delta_o)), np.inf)
        phase = np.zeros((len(delta_ag), len(delta_o)), dtype=int)

        most_stable_structure, sampling_rate = {}, 20
        for label in labels:
            gamma = np.array(self.ctx.dict_results[label]['gamma_delta_ag_delta_o'])
            mask = gamma < gamma_min
            gamma_min[mask] = gamma[mask]
            phase[mask] = int(label.replace('relax_slab_', ''))

        # Create a custom colormap from the seaborn husl palette
        husl_palette = sns.color_palette("husl", len(labels))  # Ensure enough colors for all phases
        cmap = LinearSegmentedColormap.from_list("custom_husl", husl_palette, N=len(labels))

        # Plot the phase diagram
        contour = ax1.contourf(delta_ag, delta_o, phase.T, levels=len(labels), cmap=cmap)

        ax1.spines['top'].set_visible(True)
        ax1.spines['right'].set_visible(False)  # Hide the right spine to glue it visually
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(True)
        ax1.xaxis.tick_top()
        ax1.xaxis.set_label_position('top')
        ax1.yaxis.set_ticks_position('left')
        ax1.set_ylabel(r'$\Delta \mu_{\mathrm{O}}$ (eV)', fontsize=16)
        ax1.set_xlabel(rf'$\Delta \mu_{{\mathrm{{{self.ctx.unique_symbols[0]}}}}}$ (eV)', fontsize=16, labelpad=10)
        ax1.set_ylim(-2, 0)
        ax1.set_xlim(self.ctx.lim_delta_ag, 0)
        ax1.yaxis.set_visible(False)  # Remove the y-axis

        # Uniform tick and label sizes
        ax1.tick_params(axis='both', which='major', labelsize=14)
        ax1.tick_params(axis='both', which='minor', labelsize=14)

        # Plot Delta Mu vs Pressure on the right (ax2)
        log_pressures = np.linspace(-30, 10, 10)
        T = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200]
        T = [i * ureg.K for i in T]
        G = [-0.15, -0.341, -0.548, -0.765, -0.991, -1.222, -1.460, -1.702, -1.949, -2.199, -2.453, -2.711]
        G = [i * ureg.eV for i in G]
        mu_ex = -0.03 * ureg.eV

        kB = 1.38064852e-23 * ureg["J/K"]  # Boltzmann's constant
        colors = sns.color_palette("husl", 12)  # 12 distinct colors for temperatures
        for color, t, g in zip(colors, T, G):
            delta_mu = [
                1 / 2 * (g + kB * t * np.log(10**p)) + mu_ex
                for p in log_pressures
            ]
            delta_mu = [i.magnitude for i in delta_mu]
            ax2.plot(log_pressures, delta_mu, label=f'T = {t.magnitude} K', linewidth=2.5, color=color)

        # Set axis labels
        ax2.set_xlabel(r'ln(p/p$^0$)', fontsize=16, labelpad=10)
        ax2.set_ylabel(r'$\Delta \mu_{\mathrm{O}}$ (eV)', fontsize=16)

        # Move y-axis to the right
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position('right')
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position('top')

        # Uniform tick and label sizes
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

        # Adjust layout and save the combined figure
        plt.tight_layout()
        fig.subplots_adjust(wspace=0)  # Remove space between subplots
        fig.savefig(f'{self.inputs.path_to_graphs.value}/thermo_results/surface_phase_diagram.pdf')


    # Helper Methods
    def get_vasp_builder(
        self, structure, incar_parameters, kpoint_density,
        label, description, clean_workdir=False, slab=False,
    ):
        """
        Helper method to configure the VASP builder.
        """
        builder = VaspWorkflow.get_builder()

        builder.structure = structure

        # Use only the custom INCAR parameters provided as input
        builder.parameters = Dict(dict=incar_parameters)

        # Set k-points
        kpoints = KpointsData()
        kpoints.set_cell_from_structure(structure)
        kpoints.set_kpoints_mesh_from_density(kpoint_density)
        # Adjust the k-point mesh for the slab: set the third direction to 1 (for slab calculations)
        kpoints_list = kpoints.get_kpoints_mesh()[0]
        if slab:
            kpoints_list[2] = 1
        kpoints.set_kpoints_mesh(kpoints_list)
        builder.kpoints = kpoints

        # Set code
        builder.code = self.inputs.code

        # Set computational options
        builder.options = self.inputs.computer_options

        # Set metadata
        builder.metadata.label = label
        builder.metadata.description = description

        # Set settings (parser settings)
        builder.settings = self.inputs.parser_settings

        builder.clean_workdir = Bool(clean_workdir)

        # Set potential family and mapping using the auxiliary function
        builder.potential_family = Str(self.inputs.potential_family.value)
        builder.potential_mapping = self.set_potentials()

        return builder

    def set_potentials(self):
        """
        Auxiliary method to set the potential mapping for VASP calculations.
        """
        self.report('Setting potential mapping.')
        try:
            potential_mapping = self.inputs.potential_mapping.get_dict()
            self.report('Potential mapping set successfully.')
            return Dict(dict=potential_mapping)
        except Exception as e:
            self.report(f'Failed to set potential mapping: {e}')
            raise

    def write_folder_data(self, folderdata, *files):
        for file in files:
            with folderdata.open(file, mode='rb') as source, open(file, mode='wb') as destination:
                shutil.copyfileobj(source, destination)

    def check_minimal_composition(self, bulk_structure):
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
