from vasprun import Vasprun
from pymatgen.io.vasp import Poscar
import pandas as pd

def vacancies(vasprun: Vasprun, atom: str, start: int = 1, end: int = 1, step: float|int = 0.1):
    """Creates a list of dataframes with the appropriate stoichiometry for each"""
    structure = vasprun.structure
    atom_types: list[str] = vasprun.atom_types
    positions: np.ndarray = structure.final_positions
    
    #create a dataframe with positions and atom types
    df = pd.DataFrame(positions, columns=["x", "y", "z"])
    df["atom"] = atom_types

    #get atomic ratios
    ratios = df["atom"].value_counts(normalize=True)

    #get only the atoms of the specified type
    df = df[df["atom"] == atom]

    #remove a row at random 
    df = df.sample(frac=1).reset_index(drop=True)
    df = df.drop(df.index[0])

    #get the atomic ratios again
    ratios = df["atom"].value_counts(normalize=True)
    print(ratios)
    print(df)

vasprun = Vasprun("tests/vasprun.xml")
vacancies(vasprun, "Se", 1, 1, 0.1)


