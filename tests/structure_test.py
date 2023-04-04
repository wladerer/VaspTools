from structure import Structure, pmg_structure
from pymatgen.core.structure import Structure as pmgStructure
import numpy as np
import pandas as pd

vasprun_xml_file = 'vasprun.xml'
structure = Structure(vasprun_xml_file)

def test_structure():

    #assert that atom count is a positive integer
    assert isinstance(structure.atom_count, int)
    assert structure.atom_count > 0

    #assert that unique atom count is a positive integer
    assert isinstance(structure.unique_atom_count, int)
    assert structure.unique_atom_count > 0

    #assert that atom types is a list of strings
    assert isinstance(structure.atom_types, list)
    assert all([isinstance(atom_type, str) for atom_type in structure.atom_types])

    #assert that formula is a string
    assert isinstance(structure.formula, str)

    #assert that initial basis is a 3x3 numpy array
    assert isinstance(structure.initial_basis, np.ndarray)
    assert structure.initial_basis.shape == (3, 3)

    #asser that the final basis is a 3x3 numpy array
    assert isinstance(structure.final_basis, np.ndarray)
    assert structure.final_basis.shape == (3, 3)

    #assert that initial positions is a 3xN numpy array
    assert isinstance(structure.initial_positions, np.ndarray)
    assert structure.initial_positions.shape == (structure.atom_count, 3)

    #assert that final positions is a 3xN numpy array
    assert isinstance(structure.final_positions, np.ndarray)
    assert structure.final_positions.shape == (structure.atom_count, 3)

    #assert that structure.info is a pandas dataframe
    assert isinstance(structure.info, pd.DataFrame)

    structure.plot()


def test_symmetry_tools():

    #assert that pmg_structure is a pymatgen structure
    assert isinstance(pmg_structure(structure), pmgStructure)


if __name__ == '__main__':
    test_structure()
    test_symmetry_tools()