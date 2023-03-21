# import pathlib and io modules
from pathlib import Path
import io
import sys
import pandas as pd
import numpy as np


import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# import lxml utilities
import lxml.etree as ET
import lxml.objectify as objectify


def get_kpath_labels(path: Path) -> list:
    """Returns a list of kpath labels from a KPOINTS file"""
    with open(path, 'r') as f:
        lines = f.readlines()

    # get the line that contains the kpath labels
    kpath_labels_line = lines[4:]

    # split each line into a list of strings
    kpath_info = [label.split() for label in kpath_labels_line]

    kpath_info = [label for label in kpath_info if label]
    # labels are the last element of each list, the remainig elements are the kpoint coordinates
    kpath_labels = [label[-1] for label in kpath_info]
    kpath_positions = [label[:-1] for label in kpath_info]
    kpath_positions = [np.array([float(p) for p in pos])
                       for pos in kpath_positions]

    return kpath_labels, kpath_positions


def unpack_varray(varray: ET.Element) -> np.ndarray:
    """Unpacks a varray element into a numpy array"""
    # get all v elements
    v_elements = varray.findall('v')
    # split the text of each v element into a list of strings
    v_strings = [v.text.split() for v in v_elements]
    # convert the list of strings into a list of floats
    v_floats = [[float(s) for s in v] for v in v_strings]
    # convert the list of floats into a numpy array
    varray_array = np.array(v_floats)

    return varray_array


def unpack_rarray(rarray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    # get all v elements
    r_elements = rarray.findall('r')
    # split the text of each v element into a list of strings
    r_strings = [r.text.split() for r in r_elements]
    # convert the list of strings into a list of floats
    r_floats = [[float(s) for s in r] for r in r_strings]
    # convert the list of floats into a numpy array
    rarray_array = np.array(r_floats)

    return rarray_array


def merge_discontinuities(labels: list[str]) -> list[str]:
    '''Merges discontinuities in a list of labels'''
    # if label n == n+1, keep only label n
    # if label n != n+1, keep both labels
    labels = [label for i, label in enumerate(
        labels[:-1]) if label != labels[i+1]] + [labels[-1]]

    # if the label n starts with the same substring the same as the label n+1, then the label n+1 is a discontinuity
    # and should be merged with the label n using the character '|' (e.g 'S_0|S_2')
    # dont keep the label n+1
    merged_labels = []

    for i, label in enumerate(labels):
        # short circuit evaluation
        if (i != len(labels) - 1) and label[0] == labels[i+1][0]:
            merged_labels.append(label + '|' + labels[i+1])
            # remove the label n+1
            labels.pop(i+1)
        else:
            merged_labels.append(label)

    return merged_labels


class Vasprun:

    def __init__(self, path: Path):
        self.path = path
        self.root = self.read_vasprun()
        self.dos = DensityOfStates(self)

    def read_vasprun(self):
        
        try:
            tree = ET.parse(self.path)
            root = tree.getroot()
            return root
        except ET.XMLSyntaxError:
            print('Error: Invalid XML file')
            print('Error: Check if the vasprun.xml file is complete')
            sys.exit(1)

    @property
    def incar(self):
        incar = self.root.find('incar')
        return incar

    @property
    def initial_stucture_element(self):
        initial_pos_element = self.root.find('structure[@name="initialpos"]')
        initial_structure_element = initial_pos_element.find('crystal')
        return initial_structure_element

    @property
    def initial_basis(self) -> np.ndarray:
        initial_basis_varray = self.initial_stucture_element.find(
            'varray[@name="basis"]')
        initial_basis = unpack_varray(initial_basis_varray)

        return initial_basis

    @property
    def initial_reciprocal_basis(self) -> np.ndarray:
        initial_reciprocal_basis_varray = self.initial_stucture_element.find(
            'varray[@name="rec_basis"]')
        initial_reciprocal_basis = unpack_varray(
            initial_reciprocal_basis_varray)

        return initial_reciprocal_basis

    @property
    def initial_positions(self) -> np.ndarray:
        initial_positions_varray = self.root.find(
            'structure[@name="initialpos"]').find('varray[@name="positions"]')
        initial_positions = unpack_varray(initial_positions_varray)

        return initial_positions

    @property
    def final_stucture_element(self):
        final_pos_element = self.root.find('structure[@name="finalpos"]')
        final_structure_element = final_pos_element.find('crystal')

        return final_structure_element

    @property
    def final_basis(self) -> np.ndarray:
        final_basis_varray = self.final_stucture_element.find(
            'varray[@name="basis"]')
        final_basis = unpack_varray(final_basis_varray)

        return final_basis

    @property
    def final_reciprocal_basis(self) -> np.ndarray:
        final_reciprocal_basis_varray = self.final_stucture_element.find(
            'varray[@name="rec_basis"]')
        final_reciprocal_basis = unpack_varray(final_reciprocal_basis_varray)

        return final_reciprocal_basis

    @property
    def final_positions(self) -> np.ndarray:
        final_positions_varray = self.root.find(
            'structure[@name="finalpos"]').find('varray[@name="positions"]')
        final_positions = unpack_varray(final_positions_varray)

        return final_positions

    @property
    def incar(self) -> dict:
        incar_children = self.root.find('incar').getchildren()
        incar_dict = {child.attrib['name']                      : child.text for child in incar_children}
        # trim leading and trailing whitespace from the values
        incar_dict = {k: v.strip() for k, v in incar_dict.items()}

        return incar_dict

    @property
    def kpoints(self) -> np.ndarray:
        kpoints_varray = self.root.find('kpoints').find(
            'varray[@name="kpointlist"]')
        kpoints = unpack_varray(kpoints_varray)

        return kpoints

    @property
    def kpoint_weights(self) -> np.ndarray:
        kpoint_weights_varray = self.root.find(
            'kpoints').find('varray[@name="weights"]')
        kpoint_weights = unpack_varray(kpoint_weights_varray)

        return kpoint_weights

    @property
    def kpoint_grid_type(self) -> str:
        kpoint_grid_type = self.root.find(
            'kpoints').find('generation').attrib['param']

        return kpoint_grid_type

    @property
    def kpoint_linemode_divisions(self) -> int:
        kpoint_linemode_divisions = self.root.find(
            'kpoints').find('generation').find('i')
        kpoint_linemode_divisions = int(kpoint_linemode_divisions.text)

        return kpoint_linemode_divisions

    @property
    def kpath(self) -> np.ndarray:
        kpath_vlist = self.root.find('kpoints').find('generation').findall('v')
        kpath = np.array([np.array([float(s) for s in v.text.split()])
                         for v in kpath_vlist])

        return kpath

    @property
    def kpath_labels(self) -> list:

        try:

            return get_kpath_labels('KPOINTS')[0]

        except FileNotFoundError:

            return None

    @property
    def kpath_labels_merged(self) -> list:
        merged_labels = merge_discontinuities(self.kpath_labels)

        return merged_labels

    @property
    def atom_count(self) -> int:
        atom_count = int(self.root.find('atominfo').find('atoms').text)

        return atom_count

    @property
    def unique_atom_count(self) -> int:
        unique_atoms = int(self.root.find('atominfo').find('types').text)

        return unique_atoms

    @property
    def atom_types(self) -> list:
        atom_types = self.root.find('atominfo').find(
            'array[@name="atoms"]').find('set').findall('rc')
        atom_types = [atom_type.find('c').text for atom_type in atom_types]

        return atom_types

    @property
    def formula(self) -> str:
        formula_dict = {atom_type: self.atom_types.count(
            atom_type) for atom_type in self.atom_types}
        formula = ''.join(
            [f'{atom_type}{formula_dict[atom_type]}' for atom_type in formula_dict])
        # remove 1s
        formula = formula.replace('1', '')

        return formula

    @property
    def eigenvalues(self) -> pd.DataFrame:
        spins = self.root.find('calculation').find(
            'eigenvalues').find('array').find('set').findall('set')
        eigenvalues = pd.DataFrame()
        for spin in spins:
            spin_index = int(spin.attrib['comment'].split()[-1])

            kpoints = spin.findall('set')
            for kpoint in kpoints:
                kpoint_index = int(kpoint.attrib['comment'].split()[-1])
                rows = [val.text for val in kpoint.findall('r')]
                energies = [float(row.split()[0]) for row in rows]
                occupations = [float(row.split()[1]) for row in rows]
                # create a dataframe for this kpoint
                kpoint_df = pd.DataFrame(
                    {'kpoint': kpoint_index, 'spin': spin_index, 'energy': energies, 'occupation': occupations})
                # concatenate to the main dataframe
                eigenvalues = pd.concat([eigenvalues, kpoint_df])

        return eigenvalues

    @property
    def fermi_energy(self) -> float:
        fermi_energy = float(self.root.find('calculation').find(
            'dos').find('i[@name="efermi"]').text)

        return fermi_energy


class DensityOfStates:

    def __init__(self, vasprun: Vasprun):
        self.dos_element = vasprun.root.find('calculation').find('dos')
        self.total_dos_element = self.dos_element.find(
            'total').find('array').find('set')
        self.projected_dos_element = vasprun.root.find(
            'calculation').find('projected')

    @property
    def total(self) -> pd.DataFrame:
        total_dos = pd.DataFrame()
        spins = self.total_dos_element.findall('set')
        for spin in spins:
            spin_index = int(spin.attrib['comment'].split()[-1])
            vals = unpack_rarray(spin)
            headers = ['energy', 'total', 'integrated']
            df = pd.DataFrame(vals, columns=headers)
            df['spin'] = spin_index
            total_dos = pd.concat([total_dos, df])

        return total_dos

    @property
    def partial(self) -> pd.DataFrame:
        partial_dos = pd.DataFrame()
        partial_dos_element = self.dos_element.find('partial').find('array')
        headers = [
            field.text.strip() for field in partial_dos_element.findall('field')]
        ions = partial_dos_element.find('set').findall('set')
        for ion in ions:
            ion_index = int(ion.attrib['comment'].split()[-1])

            spins = ion.findall('set')
            for spin in spins:
                spin_index = int(spin.attrib['comment'].split()[-1])
                vals = unpack_rarray(spin)
                df = pd.DataFrame(vals, columns=headers)
                df['ion'] = ion_index
                df['spin'] = spin_index
                partial_dos = pd.concat([partial_dos, df])

        return partial_dos

    @property
    def projected(self) -> pd.DataFrame:
        projected_dos = pd.DataFrame()
        projected_dos_element = self.projected_dos_element.find('array')
        headers = [
            field.text.strip() for field in projected_dos_element.findall('field')]


        spins = projected_dos_element.find('set').findall('set')
        for spin in spins:
            # spin index is weird, it no longer has a space between the word and the number, so just take the last character
            spin_index = int(spin.attrib['comment'][-1])

            kpoints = spin.findall('set')
            for kpoint in kpoints:
                kpoint_index = int(kpoint.attrib['comment'].split()[-1])

                bands = kpoint.findall('set')
                for band in bands:
                    band_index = int(band.attrib['comment'].split()[-1])
                    vals = unpack_rarray(band)
                    df = pd.DataFrame(vals, columns=headers)
                    df['band'] = band_index
                    df['kpoint'] = kpoint_index
                    df['spin'] = spin_index
                    projected_dos = pd.concat([projected_dos, df])

        return projected_dos

