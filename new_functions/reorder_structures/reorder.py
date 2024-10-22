from aiida.orm import StructureData
from pymatgen.ext.matproj import MPRester
from ase.io import read, write
from ase import Atoms
from aiida import load_profile
load_profile()

# Define your Materials Project API key
MATERIALS_PROJECT_API_KEY = 'DE9BQot894Tah2QqkkCS9lFh7vpjkvWk'

def reorder_atoms(structure: StructureData) -> StructureData:
    """
    Reorder atoms in the input structure to match a common reference order.

    Parameters:
        structure (StructureData): The input structure in AiiDA's StructureData format.

    Returns:
        StructureData: A new StructureData object with atoms reordered.
    """
    # Get ASE structure from StructureData
    ase_structure: Atoms = structure.get_ase()

    # Get the chemical formula of the structure in reduced form
    formula = ase_structure.get_chemical_formula(mode='hill')

    # Use Materials Project API to find a common order for the atoms
    with MPRester(MATERIALS_PROJECT_API_KEY) as m:
        entries = m.get_entries(formula)
        if not entries:
            raise ValueError(f"No entries found for formula {formula} in Materials Project database.")
        
        # Use the first entry as reference (could be improved for more sophisticated selection)
        ref_structure = m.get_structure_by_material_id(entries[0].material_id)
        ref_species = [str(site.specie) for site in ref_structure]

    # Create a mapping from element to indices in the ASE structure
    ase_species = [ase_structure.get_chemical_symbols()[i] for i in range(len(ase_structure))]
    sorted_indices = sorted(range(len(ase_species)), key=lambda i: ref_species.index(ase_species[i]))

    # Reorder atoms in the ASE structure
    reordered_structure = ase_structure[sorted_indices]

    # Convert back to AiiDA StructureData
    reordered_structuredata = StructureData(ase=reordered_structure)

    write('reordered_POSCAR', reordered_structuredata.get_ase(), format='vasp')

    return reordered_structuredata

def main():
    """
    Main function to launch the reorder_atoms function from a POSCAR file.
    """
    # Read the POSCAR file to get the ASE structure
    poscar_path = 'POSCAR'  # Change this to your POSCAR file path
    ase_structure = read(poscar_path)

    # Convert ASE structure to AiiDA StructureData
    structure = StructureData(ase=ase_structure)

    # Reorder atoms in the structure
    reordered_structure = reorder_atoms(structure)

    # Print the reordered structure or save it to a file
    print(reordered_structure.get_ase())

if __name__ == "__main__":
    main()
