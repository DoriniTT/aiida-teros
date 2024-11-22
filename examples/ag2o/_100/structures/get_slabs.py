from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read, write
from ase.visualize import view

# Load the primitive structure from the .vasp file
primitive_structure = Poscar.from_file("ag2o.vasp").structure
# Define the Miller indices for the 110 surface
miller_indices = (1, 0, 0)
# Define the minimum and maximum slab thickness
min_slab_thickness = 15 # in Angstroms
vacuum = 15  # in Angstroms
# Generate the slabs
slab_generator = SlabGenerator(primitive_structure, miller_indices, min_slab_thickness, vacuum, lll_reduce=True, center_slab=True)
slabs = slab_generator.get_slabs(symmetrize=True)

all_structures = []
for n, slab in enumerate(slabs):
    slab = slab.get_orthogonal_c_slab()
    ase_atoms = AseAtomsAdaptor().get_atoms(slab) 
    all_structures.append(ase_atoms)
    write('structure_'+str(n+1)+'.vasp', ase_atoms)
#view(all_structures)
