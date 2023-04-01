from vasprunio import read_vasprun, unpack_varray, unpack_rarray
from pathlib import Path
import numpy as np
import pandas as pd
import pickle
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


class DensityOfStates:

    def __init__(self, vasprun_xml_file: Path | str):
        self.path = vasprun_xml_file
        self.vasprun_root = read_vasprun(vasprun_xml_file)

        try:
            
            self.dos_element = self.vasprun_root.find('calculation').find('dos')
            self.total_dos_element = self.dos_element.find('total').find('array').find('set')
            
        
        except AttributeError:
            print('No DOS data found in vasprun.xml')

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
        '''Returns a dataframe of the projected density of states.'''
        partial_dos = pd.DataFrame()
        total_dfs = []
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
                total_dfs.append(df)

        partial_dos = pd.concat([partial_dos, total_dfs])

        return partial_dos


class ElectronicStructure:
    '''Class for reading in kpoint and band structure info from vasprun.xml file'''

    def __init__(self, vasprun_xml_file: str | Path) -> None:        
        self.vasprun_root = read_vasprun(vasprun_xml_file)
        self.kpoints_root = self.vasprun_root.find('kpoints')
        self.projected_dos_element = self.vasprun_root.find('calculation').find('projected')
        
    @property
    def fermi_energy(self) -> float:
        fermi_energy = float(self.vasprun_root.find('calculation').find('dos').find('i[@name="efermi"]').text)

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

    
    @property
    def data(self) -> pd.DataFrame:
        '''Returns a dataframe of the electronic structure information from a vasprun.xml file.'''

    
        projected_dos_element = self.projected_dos_element.find('array')
        headers = [field.text.strip() for field in projected_dos_element.findall('field')]
        spins = [spin for spin in projected_dos_element.find('set').findall('set')]
        n_spins = len(spins)
        kpoints = [kpoint.findall('set') for kpoint in spins]
        n_kpoints = len(kpoints[0])
        band_lists = [band.findall('set') for kpoint in kpoints for band in kpoint]
        band_arrays = [unpack_rarray(band) for band_list in band_lists for band in band_list]
        n_bands = len(band_lists[0])
        n_ions = len(band_arrays[0])
        ion_indices = np.arange(1, len(band_arrays[0]) + 1)
        spin_indices = np.arange(1, n_spins + 1)
        kpoint_indices = np.arange(1, n_kpoints + 1)
        band_indices = np.arange(1, n_bands + 1)
        projected_dos = pd.DataFrame(np.concatenate(band_arrays), columns=headers)


        projected_dos['ion'] = np.tile(ion_indices, n_spins * n_kpoints * n_bands)
        projected_dos['spin'] = np.tile(np.repeat(spin_indices, n_kpoints * n_bands), n_ions)
        projected_dos['kpoint'] = np.tile(np.repeat(kpoint_indices, n_bands), n_ions * n_spins)
        projected_dos['band'] = np.tile(band_indices, n_ions * n_spins * n_kpoints)

        print(f'Ion order is: {np.tile(ion_indices, n_spins * n_kpoints * n_bands)}')
        print(f'Spin order is: {np.tile(np.repeat(spin_indices, n_kpoints * n_bands), n_ions)}')
        print(f'Kpoint order is: {np.tile(np.repeat(kpoint_indices, n_bands), n_ions * n_spins)}')
        print(f'Band order is: {np.tile(band_indices, n_ions * n_spins * n_kpoints)}')

        return projected_dos

    @property
    def cprojected(self) -> pd.DataFrame:
        ''' Checks to see if a vrproj.pkl file exists and if so, returns the projected density of states from that file. If not, it creates the file and returns the projected density of states.'''

        if os.path.isfile('vrproj.pkl') and os.path.getsize('vrproj.pkl') > 0:
            with open('vrproj.pkl', 'rb') as f:
                projected_dos = pickle.load(f)
        else:
            projected_dos = self.projected
            with open('vrproj.pkl', 'wb') as f:
                pickle.dump(projected_dos, f)

        return projected_dos


vasprun_xml_file = '/home/wladerer/github/VaspTools/tests/vasprun.xml'
es = ElectronicStructure(vasprun_xml_file).data.head(100)
print(es)