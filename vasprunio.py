from pathlib import Path
import lxml.etree as ET
import numpy as np

def read_vasprun(xml_file_path: str | Path):

        try:
            tree = ET.parse(xml_file_path)
            root = tree.getroot()
            return root
        except ET.XMLSyntaxError:

            print('ERROR: Invalid XML file')
            print('ERROR: Check if the vasprun.xml file is complete')
            print('ERROR: This can either be due to a crash or an incomplete calculation')

        except FileNotFoundError:
                 
            print(f'ERROR: vasprun xml file {xml_file_path} not found. Maybe you\'re missing a forward slash')
        
        except Exception as e:
            print(e)


def unpack_varray(varray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    v_elements = varray.findall('v')
    v_strings = [r.text for r in v_elements]
    v_floats = np.array([np.fromstring(s, dtype=float, sep=' ')
                        for s in v_strings])
    varray_array = np.array(v_floats, dtype=float)

    return varray_array


def unpack_rarray(rarray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    r_elements = rarray.findall('r')
    r_strings = [r.text for r in r_elements]
    r_floats = np.array([np.fromstring(s, dtype=float, sep=' ')
                        for s in r_strings])
    rarray_array = np.array(r_floats, dtype=float)

    return rarray_array

