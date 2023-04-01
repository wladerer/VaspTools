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
    """Unpacks a varray element into a numpy array"""
    # Extract the text content of the <v> tags and split it into a list of strings
    v_strs = varray.xpath('./v/text()')

    # Convert the list of strings to a numpy array
    v_array = np.fromiter(v_strs, dtype=np.float)

    return v_array

def unpack_rarray(rarray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    r_elements = rarray.findall('r')
    r_strings = [r.text for r in r_elements]
    r_floats = np.array([np.fromstring(s, dtype=float, sep=' ')
                        for s in r_strings])
    rarray_array = np.array(r_floats, dtype=float)

    return rarray_array


