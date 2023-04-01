from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from vasprun import Vasprun

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt


def pmg_structure(vasprun: Vasprun):
    '''Returns pymatgen structure object'''
    lattice = vasprun.structure.final_basis
    species = vasprun.atom_types
    coords = vasprun.structure.final_positions
    return Structure(lattice, species, coords)

def get_symmetry(vasprun: Vasprun) -> list[str,int]:
    '''Returns symmetry of the structure'''
    structure = pmg_structure(vasprun) 
    analyzer = SpacegroupAnalyzer(structure)
    return analyzer.get_space_group_symbol(), analyzer.get_space_group_number()

def get_high_symm_kpath(vasprun: Vasprun) -> dict:
    '''Returns high symmetry kpoints'''
    structure = pmg_structure(vasprun)
    kpoints = HighSymmKpath(structure).kpath
    return kpoints

def create_paralleliped(vasprun: Vasprun, scale: int = 1) -> list:
    '''Creates a paralleliped using the basis vectors'''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #get basis vectors
    b1, b2, b3 = vasprun.structure.final_reciprocal_basis
    #create vertices
    v = scale * np.array([[0,0,0], b1, b2, b3, b1+b2, b1+b3, b2+b3, b1+b2+b3])
    
    #create a scatter3d
    ax.scatter3D(v[:,0], v[:,1], v[:,2], color='black')

    #sides of the polygons
    sides = [[v[0],v[1],v[4],v[2]], [v[0],v[1],v[5],v[3]], [v[0],v[2],v[6],v[3]], [v[1],v[4],v[7],v[5]], [v[2],v[4],v[7],v[6]], [v[3],v[5],v[7],v[6]]]

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


def add_kpoint_labels(vasprun: Vasprun, ax):
    '''Creates a scatter object of the high symmetry kpoints to be overlain on the unit cell'''
    kpoints = get_high_symm_kpath(vasprun)['kpoints']

    #kpoints = {label: [x,y,z]}
    
    coords = np.array(list(kpoints.values()))
    labels = list(kpoints.keys())
    
    #update labels to be surrounded by $$
    labels = [f'${label}$' for label in labels]

    ax.scatter3D(coords[:,0], coords[:,1], coords[:,2], color='black')
    for label, x, y, z in zip(labels, coords[:,0], coords[:,1], coords[:,2]):
        ax.text(x, y, z, label, color='purple', fontsize=10)

    return ax


def plot_unit_cell(vasprun: Vasprun, scale: int = 1):


    ax, fig, poly3d = create_paralleliped(vasprun, scale=scale)
    format_axes(ax, poly3d)
    
    #add kpoint labels
    add_kpoint_labels(vasprun, ax)
    #plot the unit cell
    plt.show()


