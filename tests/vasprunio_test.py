from vasprunio import read_vasprun, unpack_varray, unpack_rarray
import numpy as np

def test_read_vasprun():
    """Test read_vasprun function"""
    xml_file_path = '/home/wladerer/github/VaspTools/tests/vasprun.xml'
    root = read_vasprun(xml_file_path)
    assert root is not None

def test_unpack_varray():
    """Test unpack_varray function"""
    xml_file_path = '/home/wladerer/github/VaspTools/tests/vasprun.xml'
    root = read_vasprun(xml_file_path)
    #xpath is //modelling/structure/crystal/varray[@name='basis']
    varray = root.find('structure[@name="initialpos"]').find('crystal').find('varray[@name="basis"]')

    varray_array = unpack_varray(varray)

    #make sure that the varray is a numpy array
    assert isinstance(varray_array, np.ndarray)

    #make sure that the dimensions of the varray and the numpy array are the same
    v_elements = varray.findall('v') 
    v_strings = [v.text for v in v_elements]
    v_entries = [s.split() for s in v_strings]

    #make sure that the number of rows in the varray is the same as the number of rows in the numpy array
    assert len(v_entries) == varray_array.shape[0]

    #make sure that the number of columns in the varray is the same as the number of columns in the numpy array
    assert len(v_entries[0]) == varray_array.shape[1]

    print(varray_array)

def test_unpack_rarray():
    """Test unpack_rarray function"""
    xml_file_path = '/home/wladerer/github/VaspTools/tests/vasprun.xml'
    root = read_vasprun(xml_file_path)
    rarray = root.find('calculation/array/set/rarray')
    rarray_array = unpack_rarray(rarray)

    #make sure that the rarray is a numpy array
    assert isinstance(rarray_array, np.ndarray)

    #make sure that the dimensions of the rarray and the numpy array are the same
    r_elements = rarray.findall('r') 
    r_strings = [r.text for r in r_elements]
    r_entries = [s.split() for s in r_strings]

    #make sure that the number of rows in the rarray is the same as the number of rows in the numpy array
    assert len(r_entries) == rarray_array.shape[0]

    #make sure that the number of columns in the rarray is the same as the number of columns in the numpy array
    assert len(r_entries[0]) == rarray_array.shape[1]

    print(rarray_array)

if __name__ == '__main__':
    test_read_vasprun()
    test_unpack_varray()
