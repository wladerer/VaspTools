from vasprunio import read_vasprun, unpack_varray, unpack_rarray
from pathlib import Path
import numpy as np
import re
import os

def get_kpath_labels(path: Path) -> list:
    '''Returns a list of kpath labels from a KPOINTS file'''
    with open(path, 'r') as f:
        lines = f.readlines()

    kpath_info = [[label for label in line.split() if label != '!'] for line in lines[4:] if line.strip()]

    kpath_labels = [label[-1] for label in kpath_info]
    kpath_positions = [np.array([float(p) for p in label[:-1]]) for label in kpath_info]

    return kpath_labels, kpath_positions

def merge_discontinuities(labels: list[str]) -> list[str]:
    '''Merges discontinuities in a list of labels'''
    labels = [label for i, label in enumerate(
        labels[:-1]) if label != labels[i+1]] + [labels[-1]]

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


def handle_greek_characters(labels: list[str]) -> list[str]:
    '''Handles greek characters in labels'''
    initial_labels = [label.upper() for label in labels]
    #handles just gamma for now - can either be G or GAMMA
    labels = [label.replace('GAMMA', r'\Gamma') for label in initial_labels]
    labels = [label.replace('G', r'\Gamma') for label in labels if label != r'\Gamma']


    return labels

def make_labels(labels):
    '''Formats the labels for the kpath plot.'''
    #remove backslashes
    labels = [label.replace('\\', '') for label in labels]
    labels = handle_greek_characters(labels)
    labels = [r'$' + label + r'$' for label in labels]
    

    return labels


class ElectronicStructure:
    '''Class for reading in kpoint and band structure info from vasprun.xml file'''

    def __init__(self, vasprun_xml_file: str | Path) -> None:        
        self.vasprun_xml_file = vasprun_xml_file
        self.vasprun_root = read_vasprun(self.vasprun_xml_file)
        self.kpoints_root = self.vasprun_root.find('kpoints')
        
    @property
    def fermi_energy(self) -> float:
        fermi_energy = float(self.root.find('calculation').find(
            'dos').find('i[@name="efermi"]').text)

        return fermi_energy

    @property
    def kpoints(self) -> np.ndarray:
        kpoints_varray = self.kpoints_root.find('varray[@name="kpointlist"]')
        kpoints = unpack_varray(kpoints_varray)

        return kpoints

    @property
    def kpoint_weights(self) -> np.ndarray:
        kpoint_weights_varray = self.kpoints_root.find('varray[@name="weights"]')
        kpoint_weights = unpack_varray(kpoint_weights_varray)

        return kpoint_weights


    def _kpoint_grid_type(self) -> str:
        kpoint_grid_type = self.kpoints_root.find('generation').attrib['param']

        return kpoint_grid_type

    @property
    def kpath_linemode_divisions(self) -> int:
        kpath_linemode_divisions = int(self.kpoints_root.find('generation').find('i').text)

        return kpath_linemode_divisions

    @property
    def kpath(self) -> np.ndarray:
        kpath_vlist = self.kpoints_root.find('generation').findall('v')
        kpath = np.array([np.array([float(s) for s in v.text.split()])
                         for v in kpath_vlist])

        return kpath

    @property
    def kpath_labels(self) -> list:

        try:
            # append KPOINTS to self.path
            kpoints_path = os.path.join(os.path.dirname(self.path), 'KPOINTS')
            labels = get_kpath_labels(kpoints_path)[0]
            merged_labels = merge_discontinuities(labels)
            formatted_labels = make_labels(merged_labels)
            
            return formatted_labels

        except FileNotFoundError:

            return None

    def write_kpoints(self, path: Path):
        '''Writes an identical kpoint file for the current calculation'''
        kpoints = self.kpoints

        with open(path, 'w') as f:
            f.write('Resurrected KPOINTS file\n')
            f.write(f'{len(kpoints)}\n')
            f.write(f'Cartesian\n')
            for kpoint in kpoints:
                f.write(f'{kpoint[0]:12f} {kpoint[1]:12f} {kpoint[2]:12f}')
                f.write('\n')