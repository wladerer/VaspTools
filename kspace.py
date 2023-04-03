from vasprunio import read_vasprun, unpack_varray, unpack_rarray
from structure import Structure
from pathlib import Path
import numpy as np
import pandas as pd
import pickle
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
    

class BandStructure:

    def __init__(self, structure: Structure, ):
        self.structure = Structure
        self.kpoint_elements = structure.root.xpath("calculation/projected/eigenvalues/array/set/set/set")
        self.kpoint_values = self.get_kpoints()

    def get_kpoints(self) -> pd.DataFrame:
        '''Returns a dataframe of the kpoints and their energies, bands, positions, and occupations'''
        kpoint_arrays = [unpack_rarray(kpoint)
                         for kpoint in self.kpoint_elements]
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

    def merge_kpoint_pos_and_values(self) -> pd.DataFrame:

        # convert self.vasprun.kpoints to a dataframe with headers x y z
        kpoint_coords = pd.DataFrame(
            self.vasprun.kpoints, columns=['x', 'y', 'z'])
        kpoint_coords['kpoint'] = np.arange(1, len(kpoint_coords) + 1)

        # merge the two dataframes
        kpoints = pd.merge(self.kpoint_values, kpoint_coords, on='kpoint')

        return kpoints

    @property
    def kpoints(self) -> pd.DataFrame:
        # if vr.pkl exists and is not empty, return the kpoints from that file

        if os.path.isfile('vr.pkl') and os.path.getsize('vr.pkl') > 0:
            with open('vr.pkl', 'rb') as f:
                kpoints = pickle.load(f)
        else:
            kpoints = self.merge_kpoint_pos_and_values()

        return kpoints

    def band(self, bands: int) -> pd.DataFrame:
        return self.kpoints[self.kpoints['band'] == bands]

    @property
    def fermi_surface(self) -> list[pd.DataFrame]:
        # get kpoints where the energy is roughly equal to the fermi energy

        fermi_surface = self.kpoints[np.isclose(
            self.kpoints['energy'], self.vasprun.fermi_energy, rtol=0.05)]

        if fermi_surface.empty:
            print('No kpoints found on the fermi surface')
            print('Retrying at a lower tolerance')
            fermi_surface = self.kpoints[np.isclose(
                self.kpoints['energy'], self.vasprun.fermi_energy, rtol=0.08)]
            if fermi_surface.empty:
                raise ValueError('No kpoints found on the fermi surface')

        # get the bands that are on the fermi surface
        bands = fermi_surface['band'].unique()

        fermi_surface_bands = []
        for band in bands:
            fermi_surface_bands.append(
                fermi_surface[fermi_surface['band'] == band])

        return fermi_surface_bands

    def plot_3d(self, bands: list[int], title: str = None):
        import plotly.graph_objects as go

        # get the bands
        band_indices = bands
        bands = [self.band(band) for band in bands]

        # plot the bands together
        fig = go.Figure()
        fermi_energy = self.vasprun.fermi_energy
        for band in bands:

            mesh = go.Mesh3d(x=2*band['x'], y=2*band['y'], z=band['energy'] -
                             fermi_energy, showscale=False)
            fig.add_trace(mesh)

        # add a legend indicating the band number
        fig.update_scenes(xaxis_title='kx', yaxis_title='ky',
                          zaxis_title='E - Ef (eV)')


        if title == None:
            title = f'3D Band Structure of {self.vasprun.formula} (bands {", ".join(str(band) for band in band_indices)})'

        #update the color of the traces 
        for i in range(len(bands)):
            fig.data[i].update(colorscale='Viridis', showscale=False)

        fig.update_layout(title=title, title_x=0.5)

        fig.show()

    def plot_1d(self, axis: str = 'x', emin: float = None, emax: float = None, bands: list[int] = None, title: str = None):
        '''Plots the band structure along a given axis'''

        kpoints = self.kpoints
        bands = [self.kpoints[self.kpoints['band'] == band]
                 for band in kpoints['band'].unique()]

        if emin == None:
            emin = kpoints['energy'].min()
        if emax == None:
            emax = kpoints['energy'].max()

        # plot each band within the energy range
        import plotly.graph_objects as go

        bands = [band[band['energy'] < emax] for band in bands]
        bands = [band[band['energy'] > emin] for band in bands]

        fig = go.Figure()
        for band in bands:

            line_plot = go.Scatter(
                x=band['x'], y=band['energy'], mode='lines', line=dict(color='black', width=1))

            fig.add_trace(line_plot)

        fig.show()

    def plot(self, emin: float = None, emax: float = None, bands: list[int] = None, title: str = None):
        kpath = self.vasprun.kpath
        kpoints = self.kpoints
        bands = [self.kpoints[self.kpoints['band'] == band]
                 for band in kpoints['band'].unique()]

        if emin == None:
            emin = kpoints['energy'].min()
        if emax == None:
            emax = kpoints['energy'].max()

        # plot each band within the energy range
        import plotly.graph_objects as go

        n_subplots = len(kpath) - 1
        paths = [kpath_loc for kpath_loc in kpath]
        path_points = []
        for i, path in enumerate(paths):

            if i == len(paths) - 1:
                break

            start = path
            end = paths[i+1]
            path_point = (start, end)
            path_points.append(path_point)

    def contour_plot(self, band=None, emin: float = None, emax: float = None, bands: list[int] = None, title: str = None):
        '''Plots the band structure as a contour plot'''
        kpoints = self.kpoints

        if band == None:

            # get the band that is closest to the fermi energy
            band = kpoints[np.isclose(
                kpoints['energy'], self.vasprun.fermi_energy, rtol=0.05)]
            band = band[band['band'] == band['band'].max()]

        else:
            band = self.band(band)

        # plot kx vs ky with energy as the color in a contour plot with contour smoothing
        import plotly.graph_objects as go

        kx = band['x']
        ky = band['y']
        energy = band['energy'] - self.vasprun.fermi_energy

        fig = go.Figure(data=go.Contour(x=kx, y=ky, z=energy, colorscale='rdbu',
                        contours=dict(coloring='heatmap'), line_smoothing=1.2))
        # label kx and ky
        fig.update_xaxes(title_text=r'$k_x$')
        fig.update_yaxes(title_text='$k_y$')
        # label energy
        fig.update_layout(coloraxis_colorbar=dict(title=r'$E - E_f$ (eV)'))
        fig.show()

    def fermi_energy_heatmap(self, title: str = None):
        '''Plots the fermi surface'''
        kpoints = self.kpoints
        # get bands within the energy range fermi_energy - 0.2 and fermi_energy + 0.2
        fermi_surface = kpoints[np.isclose(
            kpoints['energy'], self.vasprun.fermi_energy, rtol=0.05)]
        # plot kx vs ky with energy as the color in a contour plot with contour smoothing


        kx = fermi_surface['x']
        ky = fermi_surface['y']
        kz = fermi_surface['z']
        energy = fermi_surface['energy'] - self.vasprun.fermi_energy

        
        fig = make_subplots(rows=1, cols=3)

        fig.add_trace(go.Heatmap(x=kx, y=ky, z=energy, colorscale='rdbu', showscale=True,
                      connectgaps=True, zsmooth='best', name=r'$k_x ^ k_y$'), row=1, col=1)
        fig.add_trace(go.Heatmap(x=kx, y=kz, z=energy, colorscale='rdbu', showscale=True,
                      connectgaps=True, zsmooth='best', name=r'$k_x ^ k_z$'), row=1, col=2)
        fig.add_trace(go.Heatmap(x=ky, y=kz, z=energy, colorscale='rdbu', showscale=True,
                      connectgaps=True, zsmooth='best', name=r'$k_y ^ k_z'), row=1, col=3)

        # label kx and ky
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=1)
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=2)
        fig.update_xaxes(title_text=r'$k_y$', row=1, col=3)
        fig.update_yaxes(title_text=r'$k_y$', row=1, col=1)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=2)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=3)

        fig.update_traces(
            hovertemplate='Energy: %{z:.2f}\n kx: %{x:.2f}\n ky: %{y:.2f}\n kz: %{z:.2f}')

        # add title to figure
        if title == None:
            title = f'Energy Heatmap for {self.vasprun.formula} near the Fermi Energy: {self.vasprun.fermi_energy:.2f} eV'
        fig.update_layout(title_text=title, title_x=0.5)
        fig.show()

    def fermi_heatmap(self, title: str = None):
        '''Plots the fermi surface'''
        kpoints = self.kpoints
        # get bands within the energy range fermi_energy - 0.2 and fermi_energy + 0.2
        fermi_surface = kpoints[np.isclose(
            kpoints['energy'], self.vasprun.fermi_energy, rtol=0.08)]
        # plot kx vs ky with energy as the color in a contour plot with contour smoothing
        import plotly.graph_objects as go

        kx = fermi_surface['x']
        ky = fermi_surface['y']
        kz = fermi_surface['z']
        occupation = fermi_surface['occupation']

        # make subplots for kx vs ky, kx vs kz, and ky vs kz
        fig = make_subplots(rows=1, cols=3)

        fig.add_trace(go.Heatmap(x=kx, y=ky, z=occupation, showscale=True,
                      connectgaps=True, zsmooth='best', name=r'$k_x ^ k_y$'), row=1, col=1)
        fig.add_trace(go.Heatmap(x=kx, y=kz, z=occupation, showscale=True,
                      connectgaps=True, zsmooth='best', name=r'$k_x ^ k_z$'), row=1, col=2)
        fig.add_trace(go.Heatmap(x=ky, y=kz, z=occupation, showscale=True,
                      connectgaps=True, zsmooth='best', name=r'$k_y ^ k_z'), row=1, col=3)

        # label kx and ky
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=1)
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=2)
        fig.update_xaxes(title_text=r'$k_y$', row=1, col=3)
        fig.update_yaxes(title_text=r'$k_y$', row=1, col=1)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=2)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=3)

        fig.update_traces(hovertemplate='Occupancy: %{z:.2f}')

        # add title to figure
        if title == None:
            title = f'Occupation Heatmap for {self.structure.formula} near the Fermi Energy: {self.vasprun.fermi_energy:.2f} eV'
        fig.update_layout(title_text=title, title_x=0.5)
        fig.show()


class ElectronicStructure:
    '''Class for reading in kpoint and band structure info from vasprun.xml file'''

    def __init__(self, vasprun_xml_file: str | Path) -> None:        
        self.vasprun_root = read_vasprun(vasprun_xml_file)
        self.kpoints_root = self.vasprun_root.find('kpoints')
        self.kpoint_energy_root = self.vasprun_root.xpath("calculation/projected/eigenvalues/array/set/set/set")
        self.projected_dos_element = self.vasprun_root.find('calculation').find('projected')
        
    @property
    def fermi_energy(self) -> float:
        fermi_energy = float(self.vasprun_root.find('calculation').find('dos').find('i[@name="efermi"]').text)

        return fermi_energy

    @property
    def kpoint_coordinates(self) -> np.ndarray:
        kpoints_varray = self.kpoints_root.find('varray[@name="kpointlist"]')
        kpoint_coordinates = unpack_varray(kpoints_varray)

        return kpoint_coordinates

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
        
    @property
    def kpoint_energies(self) -> pd.DataFrame:
        kpoint_arrays = [unpack_rarray(kpoint) for kpoint in self.kpoint_energy_root]
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
        kpoint_dict = {'kx': kx, 'ky': ky, 'kz': kz, 'kpoint_grid_type': self._kpoint_grid_type(), 'kpath_linemode_divisions': self.kpath_linemode_divisions}
    
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
        headers = [field.text.strip() for field in projected_dos_element.findall('field')]
        spins = [spin for spin in projected_dos_element.find('set').findall('set')]
        n_spins: int = len(spins)
        kpoints = [kpoint.findall('set') for kpoint in spins]
        n_kpoints: int = len(kpoints[0])
        band_lists = [band.findall('set') for kpoint in kpoints for band in kpoint]
        band_arrays = [unpack_rarray(band) for band_list in band_lists for band in band_list]
        n_bands: int = len(band_lists[0])
        n_ions: int = len(band_arrays[0])
        ion_indices = np.arange(1, len(band_arrays[0]) + 1)
        spin_indices = np.arange(1, n_spins + 1)
        projected_dos = pd.DataFrame(np.concatenate(band_arrays), columns=headers)


        band_index_arrays = [ n * np.repeat(1, n_ions) for n in range(1, n_bands + 1)]
        band_indices = np.tile(np.concatenate(band_index_arrays), n_spins * n_kpoints)

        kpoint_index_arrays = [ n * np.repeat(1, n_ions * n_bands) for n in range(1, n_kpoints + 1)]
        kpoint_indices = np.tile(np.concatenate(kpoint_index_arrays), n_spins)

        spin_indices = [n * np.repeat(1, n_ions * n_bands * n_kpoints) for n in range(1, n_spins + 1)] 
        spin_indices = np.concatenate(spin_indices)

        projected_dos['ion'] = np.tile(np.arange(1, n_ions + 1), n_spins * n_kpoints * n_bands)
        projected_dos['spin'] = spin_indices
        projected_dos['kpoint'] = kpoint_indices 
        projected_dos['band'] = band_indices
 
        return projected_dos
    
    @property
    def values(self):
        '''Mergeds projected_dos, kpoint_energies, and kpoint_positions into a single dataframe'''
        projected_dos = self.projected_dos
        kpoint_energies = self.kpoint_energies
        kpoint_positions = self.kpoint_coordinates

        #convert kpoint_positions to a dataframe and add an index
        kpoint_positions = pd.DataFrame(kpoint_positions, columns=['x', 'y', 'z'])
        kpoint_positions['kpoint'] = kpoint_positions.index + 1

        #merge kpoint_energies and kpoint_positions
        kpoint_energies = kpoint_energies.merge(kpoint_positions, on='kpoint')

        #merge projected_dos and kpoint_energies
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
        mask = abs(values['energy'] - self.fermi_energy) <= tolerance  # create a boolean mask
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

