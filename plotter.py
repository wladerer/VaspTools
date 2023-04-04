import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

from kspace import ElectronicStructure
from structure import Structure

def remove_background_and_axes(fig: go.Figure) -> go.Figure:
    fig.update_layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        scene=dict(
            xaxis=dict(
                visible=False,
                showticklabels=False,
                showgrid=False,
                showspikes=False,
                showbackground=False,
                zeroline=False,
                ticks='',
            ),
            yaxis=dict(
                visible=False,
                showticklabels=False,
                showgrid=False,
                showspikes=False,
                showbackground=False,
                zeroline=False,
                ticks='',
            ),
            zaxis=dict(
                visible=False,
                showticklabels=False,
                showgrid=False,
                showspikes=False,
                showbackground=False,
                zeroline=False,
                ticks='',
            ),
        ),
    )

    return fig

def add_3d_xyz_vectors(fig: go.Figure) -> go.Figure:
    '''Adds vectors to show the x,y,z axes'''
    fig.add_trace(go.Scatter3d(x=[0, 0.1], y=[0, 0], z=[0, 0], mode='lines', line=dict(color='red', width=5)))
    fig.add_trace(go.Scatter3d(x=[0, 0], y=[0, 0.1], z=[0, 0], mode='lines', line=dict(color='green', width=5)))
    fig.add_trace(go.Scatter3d(x=[0, 0], y=[0, 0], z=[0, 0.1], mode='lines', line=dict(color='blue', width=5)))

    return fig

def filter_bands(df: pd.DataFrame, emin: float = None,  emax: float = None):

    if emin is None and emax is None:  # saves time if we don't need to filter
        return df

    if emin is None:
        emin = np.min(df['energy'])
    if emax is None:
        emax = np.max(df['energy'])

    # use the bands dataframe and keep only bands that have an average energy between emin and emax
    bands = df.groupby('band').filter(
        lambda x: emin < x['energy'].mean() < emax)

    return bands


def create_band_traces(electronics: ElectronicStructure, emin: float = None, emax: float = None):
    '''Create the band traces for the plotly figure.'''
    kpoints = electronics.kpoints
    kpoints['energy'] -= electronics.fermi_energy
    kpoints = filter_bands(kpoints, emin, emax)

    bands = kpoints.groupby('band').filter(
        lambda x: emin < x['energy'].mean() < emax)

    return bands


def add_fermi_line(fig):
    '''Adds a horizontal line at the fermi energy.'''
    fig.add_hline(y=0, line_width=0.5, line_color="grey")


def add_kpoint_divisions(fig, electronics: ElectronicStructure):
    '''Adds vertical lines at the kpoint divisions to seperate paths.'''
    divisions = electronics.kpoint_linemode_divisions
    for i in range(1, len(electronics.kpoints) + 1, divisions):
        fig.add_vline(x=i - 0.5, line_width=0.5, line_color="black")

    return divisions


def customize_bandstructure_layout(fig):
    '''Creates a custom layout for the bandstructure plot.'''
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=18),
        xaxis=dict(showline=True, linewidth=1, linecolor='black',
                   title_text=r'$\textbf{k}$'),
        yaxis=dict(showline=True, linewidth=1, linecolor='black',
                   title_text=r'$E-E_f \quad \mbox{[eV]}$'),
        showlegend=False
    )


def plot_bandstructure(electronics: ElectronicStructure, emin: float = None,  emax=None, labels=None, show=True, legend: bool = True) -> px.line:
    '''Plot the bandstructure.'''
    bands = create_band_traces(electronics, emin, emax)

    trace = go.Scatter(x=bands['kpoint'], y=bands['energy'], line_color='blue', line_width=2, line_shape='spline', line_smoothing=0.5)
    fig = px.line(bands, x='kpoint', y='energy', color_discrete_sequence=['blue'], color='band')

    add_fermi_line(fig)
    divisions = add_kpoint_divisions(fig, electronics)

    if labels != None:
        kpoint_tick_val = np.arange(
            0, len(electronics.kpoints) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))

    customize_bandstructure_layout(fig)

    if show:
        fig.show()

    return trace


def plotly_compare_bands(electronics1: ElectronicStructure, electronics2: ElectronicStructure, emin: float = -2.0, emax: float = 2.0, show: bool = True, save: bool = False, labels: list[str] = None, title: str = None, legend: bool = True):
    '''Plots two bandstructures on the same plot in different colors.'''
    bands1 = create_band_traces(electronics1, emin, emax)
    bands2 = create_band_traces(electronics2, emin, emax)

    fig = go.Figure()

    # use plotly graph objects to create the bandstructures
    for band in bands1.band.unique():
        fig.add_trace(go.Scatter(x=bands1[bands1.band==band].kpoint, y=bands1[bands1.band==band].energy, mode='lines', line=dict(color='blue'), name=f'Band {band}'))
    for band in bands2.band.unique():
        fig.add_trace(go.Scatter(x=bands2[bands2.band==band].kpoint, y=bands2[bands2.band==band].energy, mode='lines', line=dict(color='red', dash='dash'), name=f'Band {band}'))

    # add x axis as r'$k$' and y axis as r'$E-E_f$ (eV)'
    fig.update_xaxes(title_text=r'$\textbf{k}$')
    fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$')

    # add a horizontal line at the fermi energy
    fig.add_hline(y=0, line_width=0.5, line_color="grey")

    # make background white
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', font=dict(size=18))

    # make figure axes black lines
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black')

    # add vertical lines at every n*kpoint_linemode_divisions kpoint
    divisions = electronics1.kpoint_linemode_divisions
    vlines = [i - 0.5 for i in range(1, len(electronics1.kpoints) + 1, divisions)]
    fig.add_vline(x=vlines, line_width=0.5, line_color="black")

    if labels is not None:
        # create array of kpoint vals
        kpoint_tick_val = np.arange(
            0, len(electronics1.kpoints) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))


    if not legend:
        fig.update_layout(showlegend=False)

    if title is not None:
        # use latex for title

        fig.update_layout(
            title_text=r'$\text{'+title+'}$', title_font_size=18, title_x=0.5)

    if save:
        fig.write_image(title + '.png')

    if show:
        fig.show()



def compare_bands(electronics1: ElectronicStructure, electronics2: ElectronicStructure, emin: float = -2.0, emax: float = 2.0, show: bool = True, save: bool = False, labels: list[str] = None, title: str = None, legend: bool = True):
    '''Plots two bandstructures on the same plot in different colors.'''
    bands1 = create_band_traces(electronics1, emin, emax)
    bands2 = create_band_traces(electronics2, emin, emax)

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # use plotly express to create the bandstructures
    lines1 = px.line(bands1, x='kpoint', y='energy',
                     color_discrete_sequence=['blue'], color='band')
    # make red and dashed
    lines2 = px.line(bands2, x='kpoint', y='energy',
                     color_discrete_sequence=['red'], color='band')

    # upate lines2 to be thin and dashed
    for line in lines2.data:
        line.line.width = 1
        line.line.dash = 'dash'

    # add the bandstructures to the figure
    fig.add_traces(lines1.data)
    fig.add_traces(lines2.data)

    # add x axis as r'$k$' and y axis as r'$E-E_f$ (eV)'
    fig.update_xaxes(title_text=r'$\textbf{k}$')
    fig.update_yaxes(title_text=r'$E-E_f \quad \mbox{[eV]}$')

    # add a horizontal line at the fermi energy
    fig.add_hline(y=0, line_width=0.5, line_color="grey")

    # make background white
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    # make figure axes black lines
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black')
    # add vertical lines at every n*kpoint_linemode_divisions kpoint
    divisions = electronics1.kpoint_linemode_divisions
    vlines = [i - 0.5 for i in range(1, len(electronics1.kpoints) + 1, divisions)]
    fig.add_vline(x=vlines, line_width=0.5, line_color="black")

    if labels is not None:
        # create array of kpoint vals
        kpoint_tick_val = np.arange(
            0, len(electronics1.kpoints) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))

    fig.update_layout(font=dict(size=18))

    if not legend:
        fig.update_layout(showlegend=False)

    if title is not None:
        # use latex for title

        fig.update_layout(
            title_text=r'$\text{'+title+'}$', title_font_size=18, title_x=0.5)

    if save:
        fig.write_image(title + '.png')

    if show:
        fig.show()

