from pathlib import Path
import pandas as pd
import numpy as np
import os

from vasprunio import read_vasprun, unpack_varray

from plotly.subplots import make_subplots
import plotly.express as px

import pickle









class BandStructure:

    def __init__(self, vasprun: Vasprun):
        self.vasprun = vasprun
        self.kpoint_elements = vasprun.root.xpath(
            "calculation/projected/eigenvalues/array/set/set/set")
        self.kpoint_values = self.get_kpoints()

    def get_kpoints(self) -> pd.DataFrame:
        '''Returns a dataframe of the kpoints and their energies, bands, positions, and occupations'''
        kpoint_arrays = [unpack_rarray(kpoint)
                         for kpoint in self.kpoint_elements]
        # add an index to the array starting at 1 and going to the length of the array

        kpoints = pd.DataFrame()

        bands = np.arange(1, len(kpoint_arrays[0]) + 1)
        for index, kpoint_array in enumerate(kpoint_arrays):

            df = pd.DataFrame(kpoint_array, columns=['energy', 'occupation'])
            df['kpoint'] = index + 1
            df['band'] = bands
            kpoints = pd.concat([kpoints, df])

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
        import plotly.graph_objects as go

        kx = fermi_surface['x']
        ky = fermi_surface['y']
        kz = fermi_surface['z']
        energy = fermi_surface['energy'] - self.vasprun.fermi_energy

        # make subplots for kx vs ky, kx vs kz, and ky vs kz
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
            title = f'Occupation Heatmap for {self.vasprun.formula} near the Fermi Energy: {self.vasprun.fermi_energy:.2f} eV'
        fig.update_layout(title_text=title, title_x=0.5)
        fig.show()


class Geometry:

    def __init__(self, vasprun):
        self.vasprun = vasprun
        self.structure = vasprun.structure
        self.initial_positions = self.structure.initial_positions
        self.final_positions = self.structure.final_positions
        self.atoms = self.structure.atom_types
        self.positions: pd.DataFrame = self.positions
        self.final_basis = self.structure.final_basis
        self.formula = self.structure.formula

    @property
    def positions(self) -> pd.DataFrame:
        from periodictable import elements
        '''Returns a dataframe with the initial and final positions of each atom'''
        initial_positions = pd.DataFrame(
            self.initial_positions, columns=['x_i', 'y_i', 'z_i'])
        final_positions = pd.DataFrame(
            self.final_positions, columns=['x', 'y', 'z'])
        atoms = pd.DataFrame(self.atoms, columns=['atom'])
        # radii = pd.DataFrame([elements.symbol(atom).covalent_radius for atom in atoms['atom']], columns=['radii'])
        structure = pd.concat(
            [atoms, initial_positions, final_positions], axis=1)
        return structure

    def plot(self, title=None):
        import plotly.graph_objects as go

        fig = go.Figure()

        colors = [(51, 92, 103), (153, 168, 140), (255, 243, 176),
                  (224, 159, 62), (158, 42, 43), (84, 11, 14)]
        colors = [f'rgb{color}' for color in colors]
        color_dict = {atom: color for atom, color in zip(
            self.structure['atom'].unique(), colors)}

        # scale the xyz coordinates to the size of the unit cell
        scaled_x = self.structure['x']
        scaled_y = self.structure['y']
        scaled_z = self.structure['z']
        scaled_pos = pd.DataFrame([scaled_x, scaled_y, scaled_z]).T
        # plot a 3d scatter plot of the atoms, color coded by atom type

        # create a trace for each atom type
        for atom in self.structure['atom'].unique():

            atom_positions = scaled_pos[self.structure['atom'] == atom]
            atom_colors = [color_dict[atom]
                           for i in range(len(atom_positions))]

            trace = go.Scatter3d(x=atom_positions['x'], y=atom_positions['y'], z=atom_positions['z'], mode='markers', marker=dict(
                size=25, color=atom_colors), name=atom)
            fig.add_trace(trace)

        # remove background grid
        fig.update_layout(scene=dict(xaxis=dict(showbackground=False, showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(
            showbackground=False, showgrid=False, zeroline=False, showticklabels=False), zaxis=dict(showbackground=False, showgrid=False, zeroline=False, showticklabels=False)))

        fig.show()


class BsDos:

    def __init__(self, vasprun):
        self.vasprun = vasprun

    @property
    def dataframe(self):
        # check if the vr.pkl file exists, return the dataframe if it does. Else, use the .merge() method to create the dataframe
        if os.path.exists('vr.pkl') and os.path.getsize('vr.pkl') > 0:
            with open('vr.pkl', 'rb') as f:
                df = pickle.load(f)
                return df
        else:
            estruct = EStruct(self.vasprun)
            df = estruct.merge()
            return df

    def plot_3d(self, bands: list[int], title: str = None):
        import plotly.graph_objects as go

        # get the dataframe
        df = self.dataframe

        # sum all orbitals to get total density of states
        orbital_titles = ['s', 'p', 'd', 'f', 'px', 'py',
                          'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
        orbital_titles = [
            orbital for orbital in orbital_titles if orbital in df.columns]
        df['total_dos'] = df[orbital_titles].sum(axis=1)
        # normalize the total density of states
        df['total_dos'] = df['total_dos']/df['total_dos'].max()

        # get unique bands from the dataframe
        unique_bands = df['band'].unique()
        bands_indices = bands
        bands = [df[df['band'] == band] for band in bands]

        # plot the bands together
        fig = go.Figure()
        fermi_energy = self.vasprun.fermi_energy
        for band in bands:

            mesh = go.Mesh3d(x=2*band['x'], y=2*band['y'], z=band['energy'] -
                             fermi_energy, intensity=band['total_dos'], showscale=True, cmin=df['total_dos'].min(), cmax=df['total_dos'].max(), colorbar=dict(title='DOS', titleside='right'), colorscale='inferno')
            fig.add_trace(mesh)

        # add a legend indicating the band number
        fig.update_scenes(xaxis_title='kx', yaxis_title='ky',
                          zaxis_title='E - Ef (eV)')

        if title == None:
            title = f'3D Band Structure of {self.vasprun.formula} (bands {", ".join(str(band) for band in bands_indices)})'

        fig.update_layout(title=title, title_x=0.5)

        fig.show()


class EStruct:
    '''Leverages both BandStructure and DensityOfStates to provide a more complete picture of the electronic structure of a material'''

    def __init__(self, vasprun):
        self.vasprun = vasprun
        self.dos = DensityOfStates(vasprun)
        self.bs = BandStructure(vasprun)

    def merge(self) -> pd.DataFrame:
        '''Merges the band structure and density of states dataframes'''

        kpoints = self.bs.kpoints
        dos = self.dos.cprojected
        orbital_titles = ['s', 'p', 'd', 'f', 'px', 'py',
                          'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
        orbital_titles = [
            orbital for orbital in orbital_titles if orbital in dos.columns]
        dos = dos.drop('ion', axis=1)
        reduced_dos = dos.groupby(['kpoint', 'band'])[orbital_titles].sum()
        merged = pd.merge(kpoints, reduced_dos, on=['kpoint', 'band'])

        return merged

    def cache(self):
        '''Caches the merged dataframe'''
        import pickle

        # warn user that this will overwrite the cached file
        if os.path.exists('vr.pkl'):
            print('Warning: This will overwrite the cached file')

        # merge the dataframes
        merged = self.merge()

        # save the merged dataframe as a pickle file
        with open('vr.pkl', 'wb') as f:
            pickle.dump(merged, f)

    def load(self):
        '''Loads the cached dataframe'''
        import pickle
        with open('vr.pkl', 'rb') as f:
            merged = pickle.load(f)
        return merged

    def plot_1dband(self):

        # get the merged dataframe
        try:

            merged = self.load()

        except:

            merged = self.merge()

        # plot the band structure
        fig = px.line(merged, x='kpoint', y='energy', color='band', color_continuous_scale='viridis',
                      line_group='band', hover_name='band', hover_data=['kpoint', 'energy'])
        fig.update_layout(title=f'Band Structure of {self.vasprun.formula}',
                          title_x=0.5, xaxis_title='kx', yaxis_title='E - Ef (eV)')
        fig.show()


def sort_dos_by_ion(projected: pd.DataFrame, ion: int) -> pd.DataFrame:
    '''Processes the projected dos data into a dataframe.'''
    # grab ion from projected dos
    ion_dos = projected[projected['ion'] == ion]
    # get orbital titles (anything that isnt ion kpoint spin or band)
    orbital_titles = [orbital for orbital in ion_dos.columns if orbital not in [
        'ion', 'kpoint', 'spin', 'band']]
    group_by_titles = [orbital for orbital in ion_dos.columns if orbital in [
        'kpoint', 'spin', 'band']]  # makes dataframes more extensible
    reduced_dos = ion_dos.groupby(group_by_titles)[orbital_titles].sum()
    # keep spin column
    reduced_dos = reduced_dos.reset_index()

    return reduced_dos


def sort_dos_by_spin(projected: pd.DataFrame, spin: int) -> pd.DataFrame:
    '''Processes the projected dos data into a dataframe.'''
    # grab ion from projected dos
    spin_dos = projected[projected['spin'] == spin]
    # get orbital titles (anything that isnt ion kpoint spin or band)
    orbital_titles = [orbital for orbital in spin_dos.columns if orbital not in [
        'ion', 'kpoint', 'spin', 'band']]
    group_by_titles = [orbital for orbital in spin_dos.columns if orbital in [
        'ion', 'kpoint', 'spin', 'band']]  # makes dataframes more extensible
    reduced_dos = spin_dos.groupby(group_by_titles)[orbital_titles].sum()
    # keep the ion column
    reduced_dos = reduced_dos.reset_index()

    return reduced_dos


