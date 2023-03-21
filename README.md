# VaspTools

VaspTools is a relatively file agnostic utility for VASP output files. You are only required to have the `vasprun.xml` file at a minium. 
This package can be used to pull properties from a completed VASP calculation. Most properties are returned as either numpy arrays or pandas dataframes for ease of use. 

## Current Callable Attributes

| Title | Property | Type | Description |
| --- | --- | --- | --- |
| Initial Basis Vectors | `initial_basis` | `np.array` | The initial basis vectors of the unit cell. |
| Initial Reciprocal Vectors | `initial_reciprocal_basis` | `np.array` | The initial reciprocal vectors of the unit cell. |
| Initial Positions | `initial_positions` | `np.array` | The initial positions of the atoms in the unit cell (x, y, z). |
| Final Basis Vectors | `final_basis` | `np.array` | The final basis vectors of the unit cell. |
| Final Reciprocal Vectors | `final_reciprocal_basis` | `np.array` | The final reciprocal vectors of the unit cell. |
| Final Positions | `final_positions` | `np.array` | The final positions of the atoms in the unit cell (x, y, z). |
| Incar Tags | 'incar_dict' | `dict` | User supplied INCAR tags and their values. |
| Kpoints | `kpoints` | `np.array` | The kpoints used in the calculation. |
| Kpoint Weights | `kpoint_weights` | `np.array` | The kpoint weights used in the calculation. |
| Kpoint Grid Type | `kpoint_grid_type` | `str` | The type of kpoint grid used in the calculation. |
| Kpoint linemode Divisions | `kpoint_linemode_divisions` | `np.array` | The kpoint linemode divisions used in the calculation. |
| Kpath | `kpath` | `np.array` | The kpoint values used to build the kpath |
| Kpath Labels | `kpath_labels` | `list[str]` | The labels for the kpoints used to build the kpath (only available if KPOINTS file is present)|
| Merged Kpath Labels | `merged_kpath_labels` | `list[str]` | Kpath labels used for plotting (handles discontinuities). |
| Atom Count | `atom_count` | `int` | The number of atoms in the unit cell. |
| Atom Types | `atom_types` | `list[str]` | The types of atoms in the unit cell. |
| Number of Unique Atoms | `unique_atom_count` | `int` | The number of unique atoms in the unit cell. |
| Reduced Formula | `formula` | `str` | The reduced formula of the unit cell (does not handle subscripts). |
| Eigenvalues | `eigenvalues` | `pd.DataFrame` | Spin resolved eigenvalues |
| Fermi Energy | `fermi_energy` | `float` | Fermi Energy. |
| Density of States | `dos` | `DensityOfStates` | DensityOfStates is a custom class that handles total, partial, and projected DOS |


Eigenvalues, Kpoints, and Poscar/Contcar will eventually be moved to separate classes for ease of use.
___



## DensityOfStates Class

The DensityOfStates class is used to handle the density of states data. It is returned as the `dos` attribute of the `Vasprun` class. 

Warning: The projected DOS is currently very slow to parse 

### DensityOfStates Attributes
| Title | Property | Type | Description | Columns |
| --- | --- | --- | --- | --- |
| Total DOS | `total` | `pd.DataFrame` | Total DOS | `energy`, `total`, `integrated` | 
| Partial DOS | `partial` | `pd.DataFrame` | Partial DOS | `ion`, `spin`,  `energy`, and all orbital types | 
| Projected DOS | `projected` | `pd.DataFrame` | Projected DOS | `ion`, `spin`, `band`, `energy`, and all orbital types |
___

## Plotting Utilities

Since most attributes are returned as numpy arrays or pandas dataframes, plotting is relatively easy. Regardless, there are a few plotting utilities included in the package for convenience.

### Plotting Utility Functions
| Title | Function | Description | Notable Arguments |
| --- | --- | --- | --- |
| Total DOS Plot | `plot_total_dos` | Plots the total DOS. | `is_spin` for spin polarized calculations |
| Partial DOS Plot | `plot_partial_dos` | Plots the partial DOS. | `ion` specify the ion as an integer, `is_spin` for spin polarized, `is_lm` for specifying l or lm resolved plots  |

___

## Example Usage


### Plot total dos of a spin polarized calculation
```python
from vasprun import Vasprun
from plotly_plotter import plot_total_dos, plot_partial_dos

vasprun = Vasprun("vasprun.xml")
plot_total_dos(vasprun.dos.total, is_spin=True)
 ```


