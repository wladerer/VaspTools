# import pathlib and io modules
from pathlib import Path
import io
import sys
import pandas as pd
import numpy as np
import os


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


def funpack_varray(varray: ET.Element) -> np.ndarray:
    """Unpacks a varray element into a numpy array"""
    # Extract the text content of the <v> tags and split it into a list of strings
    v_strs = varray.xpath('./v/text()')

    # Convert the list of strings to a numpy array
    v_array = np.fromiter(v_strs, dtype=np.float)

    return v_array


def funpack_rarray(rarray: ET.Element) -> np.ndarray:
    """Unpacks an rarray element into a numpy array"""
    r_elements = rarray.findall('r')
    r_strings = [r.text.split() for r in r_elements]
    r_floats = [[float(s) for s in r] for r in r_strings]

    # get the original shape of the rarray
    shape_str = rarray.attrib.get('shape')
    shape = tuple(map(int, shape_str.split()))

    # convert the list of floats into a numpy array
    rarray_array = np.array(r_floats).reshape(shape)

    return rarray_array


def unpack_varray(varray: ET.Element) -> np.ndarray:
    """Unpacks a varray element into a numpy array"""
    v_elements = varray.findall('v')
    v_strings = [v.text.split() for v in v_elements]
    v_floats = [[float(s) for s in v] for v in v_strings]
    varray_array = np.array(v_floats)

    return varray_array


def unpack_rarray(rarray: ET.Element) -> np.ndarray:
    """Unpacks a rarray element into a numpy array"""
    r_elements = rarray.findall('r')
    r_strings = [r.text.split() for r in r_elements]
    r_floats = [[float(s) for s in r] for r in r_strings]
    rarray_array = np.array(r_floats)

    return rarray_array


def merge_discontinuities(labels: list[str]) -> list[str]:
    '''Merges discontinuities in a list of labels'''
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

    def read_vasprun(self):

        try:
            tree = ET.parse(self.path)
            root = tree.getroot()
            return root
        except ET.XMLSyntaxError:
            print('Error: Invalid XML file')
            print('Error: Check if the vasprun.xml file is complete')

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
        incar_dict = {child.attrib['name']: child.text for child in incar_children}
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
            # append KPOINTS to self.path
            kpoints_path = os.path.join(os.path.dirname(self.path), 'KPOINTS')
            return get_kpath_labels(kpoints_path)[0]

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

    @property
    def dos(self):
        return DensityOfStates(self)
    
    @property
    def bands(self):
        return BandStructure(self)


class DensityOfStates:

    def __init__(self, vasprun: Vasprun):

        try:
            self.dos_element = vasprun.root.find('calculation').find('dos')

            self.total_dos_element = self.dos_element.find(
                'total').find('array').find('set')
            self.projected_dos_element = vasprun.root.find(
                'calculation').find('projected')
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

    @property
    def projected(self) -> pd.DataFrame:
        ''' This is still not functional'''
        projected_dos = pd.DataFrame()
        projected_dos_element = self.projected_dos_element.find('array')
        headers = [
            field.text.strip() for field in projected_dos_element.findall('field')]

        total_dfs = []
        spin_elements = projected_dos_element.find('set').findall('set')
        spins = [ spin.findall('set') for spin in spin_elements]
        

        


        # return pd.concat(total_dfs, ignore_index=True)

    




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
        return self.merge_kpoint_pos_and_values()

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
                             fermi_energy, intensity=band['occupation'], showscale=False)
            fig.add_trace(mesh)

        # add a legend indicating the band number
        fig.update_scenes(xaxis_title='kx', yaxis_title='ky',
                          zaxis_title='E - Ef (eV)')

        if title == None:
            title = f'3D Band Structure of {self.vasprun.formula} (bands {", ".join(str(band) for band in band_indices)})'

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

    def contour_plot(self, band = None, emin: float = None, emax: float = None, bands: list[int] = None, title: str = None):
        '''Plots the band structure as a contour plot'''
        kpoints = self.kpoints
        
        if band == None:

            #get the band that is closest to the fermi energy
            band = kpoints[np.isclose(kpoints['energy'], self.vasprun.fermi_energy, rtol=0.05)]
            band = band[band['band'] == band['band'].max()]

        else:
            band = self.band(band)

        #plot kx vs ky with energy as the color in a contour plot with contour smoothing
        import plotly.graph_objects as go

        kx = band['x']
        ky = band['y']
        energy = band['energy'] - self.vasprun.fermi_energy

        fig = go.Figure(data=go.Contour( x=kx, y=ky, z=energy, colorscale='rdbu', contours=dict(coloring='heatmap'), line_smoothing=1.2))
        #label kx and ky
        fig.update_xaxes(title_text=r'$k_x$')
        fig.update_yaxes(title_text='$k_y$')
        #label energy
        fig.update_layout(coloraxis_colorbar=dict(title=r'$E - E_f$ (eV)'))
        fig.show()

    def fermi_energy_heatmap(self, title: str = None):
        '''Plots the fermi surface'''
        kpoints = self.kpoints
        #get bands within the energy range fermi_energy - 0.2 and fermi_energy + 0.2
        fermi_surface = kpoints[np.isclose(kpoints['energy'], self.vasprun.fermi_energy, rtol=0.05)]
        #plot kx vs ky with energy as the color in a contour plot with contour smoothing
        import plotly.graph_objects as go

        kx = fermi_surface['x']
        ky = fermi_surface['y']
        kz = fermi_surface['z']
        energy = fermi_surface['energy'] - self.vasprun.fermi_energy

        #make subplots for kx vs ky, kx vs kz, and ky vs kz
        fig = make_subplots(rows=1, cols=3)

        fig.add_trace(go.Heatmap(x=kx, y=ky, z=energy, colorscale='rdbu', showscale=True, connectgaps=True, zsmooth='best', name=r'$k_x ^ k_y$'), row=1, col=1)
        fig.add_trace(go.Heatmap(x=kx, y=kz, z=energy, colorscale='rdbu', showscale=True, connectgaps=True, zsmooth='best', name=r'$k_x ^ k_z$'), row=1, col=2)
        fig.add_trace(go.Heatmap(x=ky, y=kz, z=energy, colorscale='rdbu', showscale=True, connectgaps=True, zsmooth='best', name=r'$k_y ^ k_z'), row=1, col=3)
        
        #label kx and ky
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=1)
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=2)
        fig.update_xaxes(title_text=r'$k_y$', row=1, col=3)
        fig.update_yaxes(title_text=r'$k_y$', row=1, col=1)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=2)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=3)
        
        fig.update_traces(hovertemplate='Energy: %{z:.2f}\n kx: %{x:.2f}\n ky: %{y:.2f}\n kz: %{z:.2f}')

        #add title to figure
        if title == None:
            title = f'Energy Heatmap for {self.vasprun.formula} near the Fermi Energy: {self.vasprun.fermi_energy:.2f} eV'
        fig.update_layout(title_text=title, title_x=0.5)
        fig.show()

    def fermi_heatmap(self, title: str = None):
        '''Plots the fermi surface'''
        kpoints = self.kpoints
        #get bands within the energy range fermi_energy - 0.2 and fermi_energy + 0.2
        fermi_surface = kpoints[np.isclose(kpoints['energy'], self.vasprun.fermi_energy, rtol=0.08)]
        #plot kx vs ky with energy as the color in a contour plot with contour smoothing
        import plotly.graph_objects as go

        kx = fermi_surface['x']
        ky = fermi_surface['y']
        kz = fermi_surface['z']
        occupation = fermi_surface['occupation'] 

        #make subplots for kx vs ky, kx vs kz, and ky vs kz
        fig = make_subplots(rows=1, cols=3)

        fig.add_trace(go.Heatmap(x=kx, y=ky, z=occupation, showscale=True, connectgaps=True, zsmooth='best', name=r'$k_x ^ k_y$'), row=1, col=1)
        fig.add_trace(go.Heatmap(x=kx, y=kz, z=occupation, showscale=True, connectgaps=True, zsmooth='best', name=r'$k_x ^ k_z$'), row=1, col=2)
        fig.add_trace(go.Heatmap(x=ky, y=kz, z=occupation, showscale=True, connectgaps=True, zsmooth='best', name=r'$k_y ^ k_z'), row=1, col=3)
        
        #label kx and ky
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=1)
        fig.update_xaxes(title_text=r'$k_x$', row=1, col=2)
        fig.update_xaxes(title_text=r'$k_y$', row=1, col=3)
        fig.update_yaxes(title_text=r'$k_y$', row=1, col=1)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=2)
        fig.update_yaxes(title_text=r'$k_z$', row=1, col=3)
        
        fig.update_traces(hovertemplate='Occupancy: %{z:.2f}')

        #add title to figure
        if title == None:
            title = f'Occupation Heatmap for {self.vasprun.formula} near the Fermi Energy: {self.vasprun.fermi_energy:.2f} eV'
        fig.update_layout(title_text=title, title_x=0.5)
        fig.show()





class Geometry:

    def __init__(self, vasprun):
        self.vasprun = vasprun
        self.initial_positions = self.vasprun.initial_positions
        self.final_positions = self.vasprun.final_positions
        self.atoms = self.vasprun.atom_types
        self.structure: pd.DataFrame = self.get_structure()
        self.final_basis = self.vasprun.final_basis
        self.formula = self.vasprun.formula

    def get_structure(self) -> pd.DataFrame:
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
