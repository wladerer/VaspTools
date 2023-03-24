from vasprun import Vasprun
import numpy as np
vasprun = Vasprun("test.xml")
total_dos = vasprun.dos.projected
