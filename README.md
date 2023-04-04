# VaspTools

VaspTools is a relatively file agnostic utility for VASP output files. You are only required to have the `vasprun.xml` file at a minium. 
This package can be used to pull properties from a completed VASP calculation. Most properties are returned as either numpy arrays or pandas dataframes for ease of use. There are other implementations of the same concept, but I have found many of them have completely abstracted away the actual data that a physicist or chemist might want. My goal is to offer a simple interface to return data in an unperturbed and straightforward manner. 

VaspTools has a straightforward syntax that would hopefully align with the way computational chemists/physicists think. Additionally, there are some quality of life features that do in fact abstract some processes away. Quick plotting tools for density of states and band structures are provided in the `plotter` module. The formatting and energy scale conventions are primarily oriented towards chemists (energy scale relative to the fermi energy, y-axis is energy). 
___
## Installation
Currently, the package is in its infancy and I dare not put it on PyPI until it is more mature. Nevertheless, the installation process is quite straightforward.

1. Clone the git repository (you may download the zip file instead)

```bash
git clone https://github.com/wladerer/VaspTools.git
``` 

2. Change into the VaspTools directory

```bash
cd VaspTools 
```

3. Install the package in editable mode

```bash
pip install -e .
```

Hopefully there will be no errors or package conflicts. 
___

## Quick Start

All objects are initialized using information from a `vasprun.xml` file. Currently, there exist a few modules that contain classes that encapsulate the input and output of a VASP calculation.

- incar 
- jobs 
- plotter 
- structure
- kspace

Each module has one or more classes that contain the appropriate and expected data. Note that kspace is ambiguously named and I apologize for that, but `electronicstructure.py` is way too long of a file name. 

### Handling IO

There is a module called vasprunio which handles parsing the xml files passed to it. The vasprunio module then hands down an xml tree root element which is used by all other classes to pull information from the xml file. 

___
### Basic Structural Information

Structural data is given with the assumption that the user is asking for converged structural data. Retrieving information goes as follows

```python

from structure import Structure

vasprun_xml_file_path = 'vasprun.xml' #can be a path object if you're nervous about it
structure = Structure(vasprun_xml_file_path)

#calculate the volume of the inital and final structures
init = structure.initial_basis
final = structure.final_basis


init_vol = np.linalg.det(init)
final_vol = np.linalg.det(final)

#get the volume difference
vol_change = final_vol - init_vol
print(vol_change)
```

Real space coordinates, reciprocal basis vectors, and much more are available as `numpy.ndarray` objects for ease of use. 
___

## Electronic Structure 

This package differs from others in that I have decided to combine band structure and density of states information. There is too much overlap between the two, and in fact, it is significantly more informative when they are considered simultaneously. Combining the two comes at a performance cost, but is a relatively minor inconvenience. I'm currently working on Rust based Python extension for parsing the xml file using xml-rs, but that is ways away. 

To view the electronic structure information, all you need to do is initialize an `ElectronicStructure` object and call the 'values' property.

```python

from kspace import ElectronicStructure

es = ElectronicStructure('vasprun.xml')
print(es.values)

```

An example output would appear as the following

```raw
             s      py      pz      px     dxy     dyz     dz2     dxz   x2-y2  ion  spin  kpoint  band   energy  occupation         x         y    z
0       0.0000  0.0000  0.0000  0.0000  0.0084  0.0021  0.0004  0.0063  0.0027    1     1       1     1 -26.5593         1.0  0.000000  0.000000  0.0
1       0.0000  0.0000  0.0000  0.0000  0.1036  0.0174  0.0183  0.0526  0.0336    2     1       1     1 -26.5593         1.0  0.000000  0.000000  0.0
2       0.0000  0.0000  0.0000  0.0000  0.0084  0.0021  0.0004  0.0063  0.0027    3     1       1     1 -26.5593         1.0  0.000000  0.000000  0.0
3       0.0000  0.0000  0.0000  0.0000  0.1040  0.0173  0.0182  0.0527  0.0335    4     1       1     1 -26.5593         1.0  0.000000  0.000000  0.0
4       0.0000  0.0000  0.0000  0.0000  0.0084  0.0021  0.0004  0.0063  0.0027    5     1       1     1 -26.5593         1.0  0.000000  0.000000  0.0
...        ...     ...     ...     ...     ...     ...     ...     ...     ...  ...   ...     ...   ...      ...         ...       ...       ...  ...
921595  0.0001  0.0002  0.0014  0.0016  0.0000  0.0000  0.0000  0.0004  0.0000   16     1     360   128   2.5619         0.0  0.333333  0.333333  0.0
921596  0.0002  0.0006  0.0027  0.0006  0.0001  0.0003  0.0000  0.0003  0.0000   17     1     360   128   2.5619         0.0  0.333333  0.333333  0.0
921597  0.0000  0.0005  0.0004  0.0013  0.0005  0.0003  0.0000  0.0002  0.0013   18     1     360   128   2.5619         0.0  0.333333  0.333333  0.0
921598  0.0002  0.0010  0.0018  0.0007  0.0000  0.0003  0.0000  0.0002  0.0000   19     1     360   128   2.5619         0.0  0.333333  0.333333  0.0
921599  0.0000  0.0004  0.0004  0.0012  0.0000  0.0000  0.0000  0.0002  0.0000   20     1     360   128   2.5619         0.0  0.333333  0.333333  0.0

[921600 rows x 18 columns]
```

It is assumed the user knows that the `x`, `y`, and `z` values are not real space coordinates. The columns expand and contract to accommodate all lm-decomposed orbitals present in the system.

### Cached output

Due to the nature of these dataframes getting incredibly large, there are built in features to create a `.pkl` file. This should drastically improve runtime. You are welcome to contact me for advice getting this to work -- I am sympathetic to the owners of slower computers. 
___
## RAM Mangament
This entire package has been profiled using the [scalene](https://github.com/emeryberger/scalene) package created by Emery Berger, Sam Stern, and Juan Altmayer Pizzorno. I have used it to manage the performance of each module. 

Unfortunately, I have not completely addressed the exorbitant RAM consumption that comes with parsing large xml files and pandas dataframes. This should not be an issue for modern computers, but it may impact the speed in which you switch through the 200 Firefox tabs you have open. 