from vasprun import Vasprun
import pandas as pd
import plotly.express as px
import numpy as np

def filter_bands(df: pd.DataFrame, emin: float = None,  emax: float = None):
    
    if emin == None and emax == None:
        return df

    if emin == None:
        emin = np.min(df['energy'])
    if emax == None:
        emax = np.max(df['energy'])

    #use the bands dataframe and keep only bands that have an average energy between emin and emax
    bands = []
    for band in df['band'].unique():
        band_df = df[df['band'] == band]
        if np.mean(band_df['energy']) > emin and np.mean(band_df['energy']) < emax:
            bands.append(band_df)

    return pd.concat(bands)
    

def plot_bandstructure(vasprun: Vasprun, emin: float = None,  emax = None, labels=None):
    '''Plot the bandstructure.'''
    kpath = vasprun.kpath
    kpoints = vasprun.bands.kpoints
     
    kpoints['energy'] = kpoints['energy'] - vasprun.fermi_energy
    
    kpoints = filter_bands(kpoints, emin, emax)

    bands = []
    for band in kpoints['band'].unique():
        band_df = kpoints[kpoints['band'] == band]
        if np.mean(band_df['energy']) > emin and np.mean(band_df['energy']) < emax:
            bands.append(band_df)
    
    bands = pd.concat(bands)

    fig = px.line(bands, x='kpoint', y='energy', color='band')

    #add x axis as r'$k$' and y axis as r'$E-E_f$ (eV)'
    fig.update_xaxes(title_text=r'$\textbf{k}$')
    fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$')

    #remove x axis ticks
    fig.update_xaxes(showticklabels=False)

    #if labels are provided, add them to the plot
    if labels != None:
        #get the add a label at each multiple of the number of vasprun.kpath_linemode_divisions
        labels = [labels[i] for i in range(0, len(labels), vasprun.kpath_linemode_divisions)]
        #add the labels to the plot
        fig.update_xaxes(ticktext=labels, tickvals=vasprun.kpath_linemode_divisions*np.arange(len(labels)))

    #add a legend
    fig.update_layout(showlegend=True)

    fig.show()



vasprun = Vasprun('tests/linemode.xml')
plot_bandstructure(vasprun, emin=-2, emax=2)



def validate_kpoints_from_procar(vasprun: Vasprun):
    '''Brute force validation of kpoints from PROCAR and vasprun.xml.'''
    #Requires you to grep "k-point" > proc.txt PROCAR and have the same vasprun.xml as the PROCAR file.

    #read file into a list of lines
    lines = open('tests/proc.txt').readlines()
    
    #remove the first line
    lines = lines[1:]
    
    #add a space in front of each - sign
    lines = [line.replace('-', ' -') for line in lines]
    lines = [line.split() for line in lines]
    #turn line[4], line[5], line[6] into a numpy array
    array = [ np.array([float(line[4]), float(line[5]), float(line[6])]) for line in lines]
    #turn array into a numpy array
    array = np.array(array)
    
    vasprun = Vasprun('tests/linemode.xml')
    
    #check if vasprun.kpoints and array are the same
    
    print(vasprun.kpoints == array)
