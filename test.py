from vasprun import Vasprun
vasprun = Vasprun("test.xml")
vasprun.bands.plot_3d([48, 50 ,52])