from vasprun import Vasprun, merge_discontinuities
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

from dosPlotter import plot_total_dos


def filter_bands(df: pd.DataFrame, emin: float = None,  emax: float = None):

    if emin == None and emax == None:  # saves time if we don't need to filter
        return df

    if emin == None:
        emin = np.min(df['energy'])
    if emax == None:
        emax = np.max(df['energy'])

    # use the bands dataframe and keep only bands that have an average energy between emin and emax
    bands = []
    for band in df['band'].unique():
        band_df = df[df['band'] == band]
        if np.mean(band_df['energy']) > emin and np.mean(band_df['energy']) < emax:
            bands.append(band_df)

    return pd.concat(bands)


def create_band_traces(vasprun: Vasprun, emin: float = None, emax: float = None):
    '''Create the band traces for the plotly figure.'''
    kpoints = vasprun.bands.kpoints
    kpoints['energy'] = kpoints['energy'] - vasprun.fermi_energy
    kpoints = filter_bands(kpoints, emin, emax)


    bands = []
    for band in kpoints['band'].unique():
        band_df = kpoints[kpoints['band'] == band]
        if np.mean(band_df['energy']) > emin and np.mean(band_df['energy']) < emax:
            bands.append(band_df)

    bands = pd.concat(bands)
    

    return bands


def plot_bandstructure(vasprun: Vasprun, emin: float = None,  emax=None, labels=None, show=True) -> px.line:
    '''Plot the bandstructure.'''
    bands = create_band_traces(vasprun, emin, emax)

    trace = go.Scatter(x=bands['kpoint'], y=bands['energy'], line_color='blue', line_width=2, line_shape='spline', line_smoothing=0.5)
    fig = px.line(bands, x='kpoint', y='energy',
                  color_discrete_sequence=['blue'], color='band')

    # add x axis as r'$k$' and y axis as r'$E-E_f$ (eV)'
    fig.update_xaxes(title_text=r'$\textbf{k}$')
    fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$')

    # add a horizontal line at the fermi energy
    fig.add_hline(y=0, line_width=0.5, line_color="grey")

    # make background white
    fig.update_layout(plot_bgcolor='white')

    # make paper white
    fig.update_layout(paper_bgcolor='white')

    # make figure axes black lines
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black')
    # add vertical lines at every n*kpoint_linemode_divisions kpoint
    divisions = vasprun.kpoint_linemode_divisions
    for i in range(1, len(vasprun.kpoints) + 1, divisions):
        fig.add_vline(x=i - 0.5, line_width=0.5, line_color="black")

    if labels != None:
        # create array of kpoint vals
        kpoint_tick_val = np.arange(
            0, len(vasprun.kpoints) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))

    fig.update_layout(font=dict(size=18))

    if show:
        fig.show()

    return trace

def compare_bands(vasprun1: Vasprun, vasprun2: Vasprun, emin: float = -2.0, emax: float = 2.0, show: bool = True, labels: list[str] = None, title: str = None):
    '''Plots two bandstructures on the same plot in different colors.'''
    bands1 = create_band_traces(vasprun1, emin, emax)
    bands2 = create_band_traces(vasprun2, emin, emax)

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    #use plotly express to create the bandstructures
    lines1 = px.line(bands1, x='kpoint', y='energy', color_discrete_sequence=['blue'], color='band')
    #make red and dashed
    lines2 = px.line(bands2, x='kpoint', y='energy', color_discrete_sequence=['red'], color='band')

    #upate lines2 to be thin and dashed
    for i in range(len(lines2.data)):
        lines2.data[i].line.width = 1
        lines2.data[i].line.dash = 'dash'

    # add the bandstructures to the figure
    fig.add_traces(lines1.data + lines2.data)


    # add x axis as r'$k$' and y axis as r'$E-E_f$ (eV)'
    fig.update_xaxes(title_text=r'$\textbf{k}$')
    fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$')

    # add a horizontal line at the fermi energy
    fig.add_hline(y=0, line_width=0.5, line_color="grey")

    # make background white
    fig.update_layout(plot_bgcolor='white')

    # make paper white
    fig.update_layout(paper_bgcolor='white')

    # make figure axes black lines
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black')
    # add vertical lines at every n*kpoint_linemode_divisions kpoint
    divisions = vasprun1.kpoint_linemode_divisions
    for i in range(1, len(vasprun1.kpoints) + 1, divisions):
        fig.add_vline(x=i - 0.5, line_width=0.5, line_color="black")

    if labels != None:
        # create array of kpoint vals
        kpoint_tick_val = np.arange(0, len(vasprun1.kpoints) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict( tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))

    fig.update_layout(font=dict(size=18))

    if title != None:
        #use latex for title

        fig.update_layout(title_text=r'$\text{'+title+'}$', title_font_size=18, title_x=0.5)

    if show:
        fig.show()


def validate_kpoints_from_procar(vasprun: Vasprun):
    '''Brute force validation of kpoints from PROCAR and vasprun.xml.'''
    # Requires you to grep "k-point" > proc.txt PROCAR and have the same vasprun.xml as the PROCAR file.

    # read file into a list of lines
    lines = open('tests/proc.txt').readlines()

    # remove the first line
    lines = lines[1:]

    # add a space in front of each - sign
    lines = [line.replace('-', ' -') for line in lines]
    lines = [line.split() for line in lines]
    # turn line[4], line[5], line[6] into a numpy array
    array = [np.array([float(line[4]), float(line[5]), float(line[6])])
             for line in lines]
    # turn array into a numpy array
    array = np.array(array)

    vasprun = Vasprun('tests/linemode.xml')

    # check if vasprun.kpoints and array are the same

    print(vasprun.kpoints == array)


def plot_bsdos(vasprun: Vasprun, emin: float , emax: float, show=True):
    '''Plot the bandstructure and DOS.'''
    band_trace = plot_bandstructure(vasprun, emin=emin, emax=emax, show=False)
    
    dos_trace = plot_total_dos(vasprun, emin=emin - 2, emax=emax + 2, show=False)

    bsdos_fig = make_subplots(rows=1, cols=2, shared_yaxes=True, shared_xaxes=False)

    bsdos_fig.add_trace(band_trace, row=1, col=1)
    bsdos_fig.add_trace(dos_trace, row=1, col=2)

    
    bsdos_fig.update_xaxes(title_text=r'$\textbf{k}$', row=1, col=1)
    bsdos_fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$', row=1, col=1)

    bsdos_fig.update_xaxes(title_text=r'$\textbf{DOS}$', row=1, col=2)
    bsdos_fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$', row=1, col=2)

    bsdos_fig.update_layout(plot_bgcolor='white')
    bsdos_fig.update_layout(paper_bgcolor='white')

    bsdos_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', row=1, col=1)
    bsdos_fig.update_yaxes(showline=True, linewidth=1, linecolor='black', row=1, col=1)

    bsdos_fig.show()
    
    



vasprun1 = Vasprun('/home/wladerer/Projects/mpPoscar_mpKpath/SOC/vasprun.xml')
vasprun2 = Vasprun('/home/wladerer/Projects/mpPoscar_mpKpath/SOD/vasprun.xml')
vasprun3 = Vasprun('/home/wladerer/Projects/tmdPoscar_tmdKpath/SOC/vasprun.xml')
vasprun4 = Vasprun('/home/wladerer/Projects/tmdPoscar_tmdKpath/SOD/vasprun.xml')
vasprun5 = Vasprun('/home/wladerer/Projects/mpPoscar_tmdKpath/SOC/vasprun.xml')
vasprun6 = Vasprun('/home/wladerer/Projects/mpPoscar_tmdKpath/SOD/vasprun.xml')
vasprun7 = Vasprun('/home/wladerer/Projects/tmdPoscar_mpKpath/SOC/vasprun.xml')
vasprun8 = Vasprun('/home/wladerer/Projects/tmdPoscar_mpKpath/SOD/vasprun.xml')
# plot_bandstructure(vasprun1, emin=-2.0, emax=2.0, show=True)

def make_labels(labels):
    labels = [r'$' + label + r'$' for label in labels]
    #replace G with r'\Gamma'
    labels = [label.replace('G', r'\Gamma') for label in labels]
    labels = [label.replace('Gamma', r'\Gamma') for label in labels]
    return labels

mpKpath_labels = make_labels(vasprun1.kpath_labels_merged)
tmdKpath_labels = make_labels(vasprun3.kpath_labels_merged)



compare_bands(vasprun1=vasprun1, vasprun2=vasprun2, emin=-2.0, emax=2.0, show=True, labels=mpKpath_labels, title='mpPoscar-mpKpath SOC vs SOD')
compare_bands(vasprun1=vasprun3, vasprun2=vasprun4, emin=-2.0, emax=2.0, show=True, labels=tmdKpath_labels, title='tmdPoscar-tmdKpath SOC vs SOD')
compare_bands(vasprun1=vasprun5, vasprun2=vasprun6, emin=-2.0, emax=2.0, show=True, labels=tmdKpath_labels, title='mpPoscar-tmdKpath SOC vs SOD')
compare_bands(vasprun1=vasprun7, vasprun2=vasprun8, emin=-2.0, emax=2.0, show=True, labels=mpKpath_labels, title='tmdPoscar-mpKpath SOC vs SOD')

    