from kspace import ElectronicStructure, DensityOfStates, make_labels
import numpy as np
import pandas as pd

vasprun_xml_file = '/home/wladerer/github/VaspTools/tests/vasprun.xml'
electronic_structure = ElectronicStructure(vasprun_xml_file)

def test_electronic_structure():

    #asser fermi energy is a float between -100 and 100
    assert isinstance(electronic_structure.fermi_energy, float)
    assert -100 < electronic_structure.fermi_energy < 100

    #assert kpoint_coordinates is a Nx3 numpy array
    assert isinstance(electronic_structure.kpoint_coordinates, np.ndarray)
    assert electronic_structure.kpoint_coordinates.shape[1] == 3

    #assert kpoint weights is a 1xN numpy array
    assert isinstance(electronic_structure.kpoint_weights, np.ndarray)
    assert electronic_structure.kpoint_weights.shape[1] == 1

    #assert kpoint grid type is a string
    assert isinstance(electronic_structure._kpoint_grid_type(), str)

    #assert kpath_linemode_divisions is an int
    assert isinstance(electronic_structure.kpath_linemode_divisions, int)

    #assert kpath is a Nx3 numpy array
    assert isinstance(electronic_structure.kpath, np.ndarray)

    #assert that data is a pandas dataframe that is not empty
    assert isinstance(electronic_structure.values, pd.DataFrame)
    assert not electronic_structure.values.empty
    

if __name__ == '__main__':
    test_electronic_structure()
