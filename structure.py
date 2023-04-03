from vasprunio import read_vasprun, unpack_varray
from pathlib import Path
import numpy as np
import pandas as pd

from pymatgen.core import Structure as pmgStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

import plotly.graph_objects as go


class Structure:

    def __init__(self, xml_path: Path | str):
        '''Creates a Structure object from a xml file.'''
        self.root = read_vasprun(xml_path)

    def __str__(self) -> str:
        return f'{self.formula} \t '

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
        #remove whitespace
        atom_types = [atom_type.strip() for atom_type in atom_types]

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
    def _final_stucture_element(self):
        final_pos_element = self.root.find('structure[@name="finalpos"]')
        final_structure_element = final_pos_element.find('crystal')

        return final_structure_element

    @property
    def final_basis(self) -> np.ndarray:
        final_basis_varray = self._final_stucture_element.find(
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
    def formula(self) -> str:
        formula_dict = {atom_type: self.atom_types.count(
            atom_type) for atom_type in self.atom_types}
        formula = ''.join(
            [f'{atom_type}{formula_dict[atom_type]}' for atom_type in formula_dict])
        # remove 1s
        formula = formula.replace('1', '')

        return formula

    @property
    def selective_dynamics(self) -> list[list[str]]:
        '''Returns a list of lists of T/F strings indicating whether the atom is fixed in that direction'''
        selective_dynamics_element = self.root.find(
            'structure[@name="initialpos"]').find('varray[@name="selective"]')
        selective_dynamics = [v.text.split()
                              for v in selective_dynamics_element.findall('v')]

        return selective_dynamics

    def write_poscar(self, filename: str = 'POSCAR', final: bool = False, scale: float = 1.0):
        '''Writes a POSCAR file with the final or initial structure'''
        basis = self.final_basis if final else self.initial_basis
        positions = self.final_positions if final else self.initial_positions

        with open(filename, 'w') as f:
            f.write(f'{self.formula}\n{scale}\n')
            for basis_vector in basis:
                f.write(
                    f'{basis_vector[0]:8f} {basis_vector[1]:8f} {basis_vector[2]:8f}\n')

            atom_dict = {atom_type: self.atom_types.count(
                atom_type) for atom_type in self.atom_types}
            atom_line = ' '.join(atom_dict.keys())
            atom_counts = ' '.join(str(count) for count in atom_dict.values())

            f.write(f'{atom_line}\n{atom_counts}\nDirect\n')
            selective_dynamics = self.selective_dynamics
            for atom, position in enumerate(positions):
                selective = ' '.join(selective_dynamics[atom])
                f.write(
                    f'{position[0]:.8f} {position[1]:.8f} {position[2]:.8f} {selective} {self.atom_types[atom]}\n')

    @property
    def info(self) -> pd.DataFrame:
        '''Returns a pandas DataFrame with the structural information from the xml file'''
        df = pd.DataFrame(columns=['atom', 'type', 'ix', 'iy', 'iz', 'fx', 'fy', 'fz', 'seldyn_x', 'seldyn_y', 'seldyn_z'])
        for atom, atom_type in enumerate(self.atom_types):
            df.loc[atom] = [atom, atom_type, *self.initial_positions[atom], *self.final_positions[atom], *self.selective_dynamics[atom]]

        # calculate the displacements
        df['dx'] = df['fx'] - df['ix']
        df['dy'] = df['fy'] - df['iy']
        df['dz'] = df['fz'] - df['iz']

        return df 
    
    def plot(self, is_final: bool = True):
        '''Creates a 3D plot of the structure using plotly'''
        positions = self.final_positions if is_final else self.initial_positions
        fig = go.Figure(data=[go.Scatter3d(x=positions[:,0], y=positions[:,1], z=positions[:,2], mode='markers', marker=dict(size=10))])
        
        from plotly_helpers import remove_background_and_axes, add_3d_xyz_vectors
        remove_background_and_axes(fig)
        add_3d_xyz_vectors(fig)

        fig.show()

    def as_dict(self, is_final: bool = True) -> dict:
        '''Returns a dictionary of the structure'''
        positions = self.final_positions if is_final else self.initial_positions
        return {'basis': self.final_basis.tolist(), 'positions': positions.tolist(), 'atom_types': self.atom_types}


        

def pmg_structure(structure: Structure) -> pmgStructure:
    '''Returns pymatgen structure object'''
    lattice = structure.final_basis
    species = structure.atom_types
    coords = structure.final_positions
    return pmgStructure(lattice, species, coords)

def get_symmetry(structure: Structure) -> list[str,int]:
    '''Returns symmetry of the structure'''
    structure: pmgStructure = pmg_structure(structure) 
    analyzer = SpacegroupAnalyzer(structure)
    return analyzer.get_space_group_symbol(), analyzer.get_space_group_number()

def get_high_symm_kpath(structure: Structure) -> dict:
    '''Returns high symmetry kpoints'''
    structure = pmg_structure(structure)
    kpoints = HighSymmKpath(structure).kpath
    return kpoints

def create_paralleliped(structure: Structure, scale: int = 1) -> list:
    '''Creates a paralleliped using the basis vectors'''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #get basis vectors
    b1, b2, b3 = structure.final_reciprocal_basis
    #create vertices
    v = scale * np.array([[0,0,0], b1, b2, b3, b1+b2, b1+b3, b2+b3, b1+b2+b3])
    
    #create a scatter3d
    ax.scatter3D(v[:,0], v[:,1], v[:,2], color='black')

    #sides of the polygons
    sides = [
    [v[0], v[1], v[4], v[2]],
    [v[0], v[1], v[5], v[3]],
    [v[0], v[2], v[6], v[3]],
    [v[1], v[4], v[7], v[5]],
    [v[2], v[4], v[7], v[6]],
    [v[3], v[5], v[7], v[6]]
    ]

    #create a Poly3DCollection
    poly3d = Poly3DCollection(sides, linewidths=1, edgecolors='black', alpha=0.1)
    poly3d.set_facecolor('red')
    
    return ax, fig, poly3d

def format_axes(ax,poly3d):
    ax.add_collection3d(poly3d)

    #set the axes
    ax.set_xlabel(r'$\mathbf{k_x}$')
    ax.set_ylabel(r'$\mathbf{k_y}$')
    ax.set_zlabel(r'$\mathbf{k_z}$')


def add_kpoint_labels(structure: Structure, ax):
    '''Creates a scatter object of the high symmetry kpoints to be overlain on the unit cell'''
    kpoints = get_high_symm_kpath(structure)['kpoints']

    #kpoints = {label: [x,y,z]}
    
    coords = np.array(list(kpoints.values()))
    labels = list(kpoints.keys())
    
    #update labels to be surrounded by $$
    labels = [f'${label}$' for label in labels]

    ax.scatter3D(coords[:,0], coords[:,1], coords[:,2], color='black')
    for label, x, y, z in zip(labels, coords[:,0], coords[:,1], coords[:,2]):
        ax.text(x, y, z, label, color='purple', fontsize=10)

    return ax


def plot_unit_cell(structure: Structure, scale: int = 1):


    ax, fig, poly3d = create_paralleliped(structure, scale=scale)
    format_axes(ax, poly3d)
    
    #add kpoint labels
    add_kpoint_labels(structure, ax)
    #plot the unit cell
    plt.show()


