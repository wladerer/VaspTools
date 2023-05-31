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

    kpath_info = [[label for label in line.split() if label != '!']
                  for line in lines[4:] if line.strip()]

    kpath_labels = [label[-1] for label in kpath_info]
    kpath_positions = [np.array([float(p) for p in label[:-1]])
                       for label in kpath_info]

    return kpath_labels, kpath_positions


def merge_discontinuities(labels: list[str]) -> list[str]:
    '''Merges discontinuities in a list of labels'''

    new_labels = []
    last = ''
    for i, l in enumerate(labels):
        new_labels.append(l)
        if i % 2 or i == 0:
            last = l
            continue
        if l == last:
            new_labels.pop()
        else:
            z = new_labels.pop()
            new_labels[-1] = new_labels[-1] + '|' + z
        last = l
    
    return new_labels


def handle_greek_characters(labels: list[str]) -> list[str]:
    '''Handles greek characters in labels'''
    initial_labels = [label.upper() for label in labels]
    #replace GAMMA with \Gamma
    labels = [label.replace('GAMMA', r'\Gamma') for label in initial_labels]

    return labels


def make_labels(labels):
    '''Formats the labels for the kpath plot.'''
    # remove backslashes
    labels = [label.replace('\\', '') for label in labels]
    labels = handle_greek_characters(labels)
    labels = [r'$' + label + r'$' for label in labels]

    return labels


class ElectronicStructure:
    '''Class for reading in kpoint and band structure info from vasprun.xml file'''

    def __init__(self, vasprun_xml_file: str | Path) -> None:
        self.file = vasprun_xml_file
        self.vasprun_root = read_vasprun(vasprun_xml_file)
        self.kpoints_root = self.vasprun_root.find('kpoints')
        self.kpoint_energy_root = self.vasprun_root.xpath(
            "calculation/projected/eigenvalues/array/set/set/set")
        self.projected_dos_element = self.vasprun_root.find(
            'calculation').find('projected')

    def band(self, bands: int) -> pd.DataFrame:
        '''Returns a dataframe of the kpoints and energies for a given band'''
        return self.values[self.values['band'] == bands]

    @property
    def fermi_energy(self) -> float:
        fermi_energy = float(self.vasprun_root.find(
            'calculation').find('dos').find('i[@name="efermi"]').text)

        return fermi_energy

    @property
    def kpoint_coordinates(self) -> np.ndarray:
        kpoints_varray = self.kpoints_root.find('varray[@name="kpointlist"]')
        kpoint_coordinates = unpack_varray(kpoints_varray)

        return kpoint_coordinates

    @property
    def kpoint_weights(self) -> np.ndarray:
        kpoint_weights_varray = self.kpoints_root.find(
            'varray[@name="weights"]')
        kpoint_weights = unpack_varray(kpoint_weights_varray)

        return kpoint_weights

    def _kpoint_grid_type(self) -> str:
        kpoint_grid_type = self.kpoints_root.find('generation').attrib['param']

        return kpoint_grid_type

    @property
    def kpath_linemode_divisions(self) -> int:
        try:
            kpath_linemode_divisions = int(
                self.kpoints_root.find('generation').find('i').text)
        except:
            print('Warning: No linemode divisions found in KPOINTS file')
            return None

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
            kpoints_path = os.path.join(os.path.dirname(self.file), 'KPOINTS')
            labels = get_kpath_labels(kpoints_path)[0]
            merged_labels = merge_discontinuities(labels)
            formatted_labels = make_labels(merged_labels)

            return formatted_labels

        except FileNotFoundError:

            return None

    @property
    def kpoint_energies(self) -> pd.DataFrame:
        kpoint_arrays = [unpack_rarray(kpoint)
                         for kpoint in self.kpoint_energy_root]
        # add an index to the array starting at 1 and going to the length of the array

        dfs = []
        bands = np.arange(1, len(kpoint_arrays[0]) + 1)
        for index, kpoint_array in enumerate(kpoint_arrays):

            df = pd.DataFrame(kpoint_array, columns=['energy', 'occupation'])
            df['kpoint'] = index + 1
            df['band'] = bands
            dfs.append(df)

        kpoints = pd.concat(dfs)

        return kpoints

    def as_dict(self):
        '''Returns a dictionary of kpoint coordinates and generation scheme'''
        kx, ky, kz = self.kpoint_coordinates.T
        kpoint_dict = {'kx': kx, 'ky': ky, 'kz': kz,
                       'kpoint_grid_type': self._kpoint_grid_type()}

        if self.kpath_linemode_divisions != None:
            kpoint_dict['kpath_linemode_divisions'] = self.kpath_linemode_divisions

        return kpoint_dict

    def write_kpoints(self, path: Path):
        '''Writes an identical kpoint file for the current calculation'''
        kpoints = self.kpoint_coordinates

        with open(path, 'w') as f:
            f.write('Resurrected KPOINTS file\n')
            f.write(f'{len(kpoints)}\n')
            f.write(f'Cartesian\n')
            for kpoint in kpoints:
                f.write(f'{kpoint[0]:12f} {kpoint[1]:12f} {kpoint[2]:12f}')
                f.write('\n')

    @property
    def projected_dos(self) -> pd.DataFrame:
        '''Returns a dataframe of the electronic structure information from a vasprun.xml file.'''

        projected_dos_element = self.projected_dos_element.find('array')
        headers = [field.text.strip()
                   for field in projected_dos_element.findall('field')]
        spins = [spin for spin in projected_dos_element.find(
            'set').findall('set')]
        n_spins: int = len(spins)
        kpoints = [kpoint.findall('set') for kpoint in spins]
        n_kpoints: int = len(kpoints[0])
        band_lists = [band.findall('set')
                      for kpoint in kpoints for band in kpoint]
        band_arrays = [unpack_rarray(band)
                       for band_list in band_lists for band in band_list]
        n_bands: int = len(band_lists[0])
        n_ions: int = len(band_arrays[0])
        ion_indices = np.arange(1, len(band_arrays[0]) + 1)
        spin_indices = np.arange(1, n_spins + 1)
        projected_dos = pd.DataFrame(
            np.concatenate(band_arrays), columns=headers)

        band_index_arrays = [n * np.repeat(1, n_ions)
                             for n in range(1, n_bands + 1)]
        band_indices = np.tile(np.concatenate(
            band_index_arrays), n_spins * n_kpoints)

        kpoint_index_arrays = [
            n * np.repeat(1, n_ions * n_bands) for n in range(1, n_kpoints + 1)]
        kpoint_indices = np.tile(np.concatenate(kpoint_index_arrays), n_spins)

        spin_indices = [n * np.repeat(1, n_ions * n_bands * n_kpoints)
                        for n in range(1, n_spins + 1)]
        spin_indices = np.concatenate(spin_indices)

        projected_dos['ion'] = np.tile(
            np.arange(1, n_ions + 1), n_spins * n_kpoints * n_bands)
        projected_dos['spin'] = spin_indices
        projected_dos['kpoint'] = kpoint_indices
        projected_dos['band'] = band_indices

        return projected_dos

    @property
    def values(self):
        '''Mergeds projected_dos, kpoint_energies, and kpoint_positions into a single dataframe'''
        try:
            projected_dos = self.cprojected_dos
        except:
            projected_dos = self.projected_dos

            
        kpoint_energies = self.kpoint_energies
        kpoint_positions = self.kpoint_coordinates

        # convert kpoint_positions to a dataframe and add an index
        kpoint_positions = pd.DataFrame(
            kpoint_positions, columns=['x', 'y', 'z'])
        kpoint_positions['kpoint'] = kpoint_positions.index + 1

        # merge kpoint_energies and kpoint_positions
        kpoint_energies = kpoint_energies.merge(kpoint_positions, on='kpoint')

        # merge projected_dos and kpoint_energies
        values = projected_dos.merge(kpoint_energies, on=['kpoint', 'band'])

        return values

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

    def bands_crossing_fermi(self) -> pd.DataFrame:
        '''Returns a dataframe of the bands crossing the fermi level'''
        values = self.values
        fermi_energy = self.fermi_energy

        tolerance = 0.28  # set your tolerance value here
        # create a boolean mask
        mask = abs(values['energy'] - self.fermi_energy) <= tolerance
        df_new = values[mask]  # apply the mask to create a new dataframe
        unique_bands = df_new['band'].unique()  # get the unique bands
        return df_new, unique_bands


def ignore_spin(df: pd.DataFrame) -> pd.DataFrame:
    '''Returns a dataframe with the spin column removed'''
    df = df.drop('spin', axis=1)
    return df


def ignore_ion(df: pd.DataFrame) -> pd.DataFrame:
    '''Returns a dataframe with the ion column removed'''
    df = df.drop('ion', axis=1)
    return df


def speed_test():
    vasprun_xml_file = 'tests/vasprun.xml'
    # vasprun_xml_file2 = '/home/wladerer/github/VaspTools/tests/soc_kpoint_mesh_test_file.xml'

    electronic_structure = ElectronicStructure(vasprun_xml_file)
    # electronic_structure2 = ElectronicStructure(vasprun_xml_file2)

    electronic_structure.values
    # electronic_structure2.values

def generate_kpoint_mesh(kpoint_sampling: list[int]) -> np.ndarray:
    '''Generates a 3d mesh of kpoints in the BZ'''
    kx, ky, kz = kpoint_sampling
    x = np.linspace(0, 1, kx)
    y = np.linspace(0, 1, ky)
    z = np.linspace(0, 1, kz)

    xx, yy, zz = np.meshgrid(x, y, z)

    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

    return points

def generate_kpoint_plane(plane: list[int]) -> np.ndarray:
    '''Generates a plane of kpoints in the BZ'''
    x, y, z = plane
    p1 = np.array([x, 0, 0])
    p2 = np.array([0, y, 0])
    p3 = np.array([0, 0, z])

    # Calculate the normal vector of the plane
    normal = np.cross(p2 - p1, p3 - p1)

    # Define a grid of points in the x-y plane
    x, y = np.meshgrid(np.linspace(0, 1, 10), np.linspace(0, 1, 10))

    # Calculate the z-coordinate of each point on the plane
    z = (-normal[0] * x - normal[1] * y - np.dot(normal, p1)) / normal[2]

    # Stack the x, y, and z coordinates into a single array
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))


    # Check if [0, 0, 0] is in the array
    if not any(np.all(point == [0, 0, 0]) for point in points):
        # Add [0, 0, 0] to the array
        points = np.vstack((points, [0, 0, 0]))

    return points


def isosurface_bands(electronics: ElectronicStructure, energy: float) -> list[pd.DataFrame]:
    # get kpoints where the energy is roughly equal to an energy value

    isosurface = electronics.values[np.isclose(
        electronics.values['energy'], energy, rtol=0.05)]

    if isosurface.empty:
        print(f'No kpoints found wth an energy of {energy} eV')
        print('Retrying at a lower tolerance')
        isosurface = electronics.values[np.isclose(
            electronics.values['energy'], energy, rtol=0.08)]
        if isosurface.empty:
            raise ValueError(f'No kpoints found with an energy of {energy} eV')

    # get the bands that are on the surface
    bands = isosurface['band'].unique()

    isosurface_bands = []
    for band in bands:
        isosurface_bands.append(
            isosurface[isosurface['band'] == band])

    return isosurface_bands