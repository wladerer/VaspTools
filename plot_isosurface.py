from plotter import plot_fermi_heatmap, plot_band_contour, plot_band_volume
from kspace import ElectronicStructure
import numpy as np

file = '/home/wladerer/Projects/vasprun.xml'
es = ElectronicStructure(file)
bands = np.arange(64, 80)
plot_band_volume(es)
