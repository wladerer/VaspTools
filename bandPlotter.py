from vasprun import Vasprun
import pandas as pd
import plotly.express as px
px.defaults.template = 'simple_white'

def plot_bandstructure(vasprun: Vasprun, title=None, dark_mode=False, labels=None):
    '''Plot the bandstructure.'''
    if labels == None:

        try: 
            labels = vasprun.kpath_labels_merged
        except FileNotFoundError:
            print('Error: No KPOINTS file found. Please provide a KPOINTS file or provide a list of labels.')
    kpath = vasprun.kpath
    eigenvalues: pd.DataFrame = vasprun.eigenvalues
    #ignore the spin
    eigenvalues = eigenvalues.drop(columns=['spin'])
    



    

vasp = Vasprun('/home/wladerer/Projects/mpPoscar_mpKpath/SOD/vasprun.xml')
plot_bandstructure(vasp)