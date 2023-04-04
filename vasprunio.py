from pathlib import Path
import lxml.etree as ET
import numpy as np

def read_vasprun(xml_file_path: str | Path, debug = False) -> ET.Element:

    #make sure that the xml file exists
    if not Path(xml_file_path).exists():
        raise FileNotFoundError(f'Could not find file {xml_file_path}')
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    return root
        

def unpack_varray(varray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    v_strings = [r.text for r in varray.findall('v')]
    varray_array = np.array([np.fromstring(s, dtype=float, sep=' ') for s in v_strings], dtype=float)
    return varray_array


def unpack_rarray(rarray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    r_strings = [r.text for r in rarray.findall('r')]
    r_floats = np.array([np.fromstring(s, dtype=float, sep=' ') for s in r_strings], dtype=float)
    return r_floats

