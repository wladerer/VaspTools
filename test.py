from vasprun import Vasprun
import pandas as pd
vasprun = Vasprun('vasprun.xml')
kpoints: pd.DataFrame = vasprun.bands.kpoints

# print values of bands 389, 390, and 391
print(kpoints.iloc[389:392])
