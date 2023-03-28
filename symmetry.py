from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from vasprun import Vasprun

def pmg_structure(vasprun: Vasprun):
    '''Returns pymatgen structure object'''
    lattice = vasprun.structure.final_basis
    species = vasprun.atom_types
    coords = vasprun.structure.final_positions
    return Structure(lattice, species, coords)

def get_symmetry(vasprun: Vasprun) -> list[str,int]:
    '''Returns symmetry of the structure'''
    structure = pmg_structure(vasprun) 
    analyzer = SpacegroupAnalyzer(structure)
    return analyzer.get_space_group_symbol(), analyzer.get_space_group_number()

def get_high_symm_kpath(vasprun: Vasprun) -> dict:
    '''Returns high symmetry kpoints'''
    structure = pmg_structure(vasprun)
    kpoints = HighSymmKpath(structure).kpath
    return kpoints

