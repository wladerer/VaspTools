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



def create_band_traces(electronics: ElectronicStructure, emin: float = None, emax: float = None) -> pd.DataFrame:
    '''Create the band traces for the plotly figure.'''
    kpoints = electronics.values
    kpoints['energy'] -= electronics.fermi_energy
    kpoints = filter_bands(kpoints, emin, emax)

    bands = kpoints.groupby('band').filter(
        lambda x: emin < x['energy'].mean() < emax)

    return bands

def create_band_trace(electronics: ElectronicStructure, band: int):
    '''Create the band traces for the plotly figure.'''
    kpoints = electronics.values
    kpoints['energy'] -= electronics.fermi_energy

    band = kpoints[kpoints['band'] == band]

    return band

def add_fermi_line(fig):
    '''Adds a horizontal line at the fermi energy.'''
    fig.add_hline(y=0, line_width=0.5, line_color="grey")


def add_kpoint_divisions(fig, electronics: ElectronicStructure):
    '''Adds vertical lines at the kpoint divisions to seperate paths.'''
    divisions = electronics.kpath_linemode_divisions
    for i in range(1, len(electronics.kpoint_coordinates) + 2, divisions):
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
            0, len(electronics.kpoint_coordinates) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))

    customize_bandstructure_layout(fig)

    if show:
        fig.show()

    return trace

def plot_specific_bands(electronics: ElectronicStructure, bands: list[int], labels=None, show=True, legend: bool = True):
    '''Plot 1D bandstructure for the given bands.'''
    fig = go.Figure()
    bands = pd.concat([create_band_trace(electronics, band) for band in bands])

    #create the band traces , color them according to band index
    fig.add_trace(go.Scatter(x=bands['kpoint'], y=bands['energy'], line_width=2, line_shape='spline', line_smoothing=0.5, colorscale='Viridis', color=bands['energy']))



    add_fermi_line(fig)
    divisions = add_kpoint_divisions(fig, electronics)

    if labels != None:
        kpoint_tick_val = np.arange(
            0, len(electronics.kpoint_coordinates) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))
        
    customize_bandstructure_layout(fig)

    if show:
        fig.show()



def plot_2d(electronics: ElectronicStructure, bands: list[int], title: str = None):
    '''Plot the 2D band structure for the given bands.'''

    # get the bands
    band_indices = bands
    bands = [electronics.band(band) for band in bands]

    # plot the bands together
    fig = go.Figure()
    fermi_energy = electronics.fermi_energy
    for band in bands:

        mesh = go.Mesh3d(x=2*band['x'], y=2*band['y'], z=band['energy'] -
                         fermi_energy, showscale=False)
        fig.add_trace(mesh)

    # add a legend indicating the band number
    fig.update_scenes(xaxis_title='kx', yaxis_title='ky',
                      zaxis_title='E - Ef (eV)')


    if title == None:
        title = f'3D Band Structure (bands {", ".join(str(band) for band in band_indices)})'

    #update the color of the traces 
    for i in range(len(bands)):
        fig.data[i].update(colorscale='Viridis', showscale=False)

    fig.update_layout(title=title, title_x=0.5)

    fig.show()


def plotly_compare_bands(electronics1: ElectronicStructure, electronics2: ElectronicStructure, emin: float = -2.0, emax: float = 2.0, show: bool = True, save: bool = False, labels: list[str] = None, title: str = None, legend: bool = True):
    '''Plots two 1D bandstructures on the same plot in different colors.'''
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

def plot_kpath_against_energy(es: ElectronicStructure, band: int):
    '''Plot the kpath against the energy of a given band'''

    path = es.kpoint_coordinates
    energies: pd.DataFrame = es.kpoint_energies

    energies = energies[energies['band'] == band].sort_values(by=['kpoint'])['energy'].to_numpy() - es.fermi_energy

    fig = go.Figure(data=[go.Scatter3d(x=path[:, 0], y=path[:, 1], z=path[:, 2], mode='markers', marker=dict(color=energies, colorscale='inferno', showscale=True, size=5))])
    fig.update_traces(connectgaps=False)
    #add title "Band 94"
    fig.update_layout(title=f'Band {band}', title_x=0.5)

    #add axis labels x= kx, y=ky, z=kz
    fig.update_layout(scene = dict( xaxis_title='kx', yaxis_title='ky', zaxis_title='kz'))
    fig.show()

def compare_kpaths(kpath1, kpath2):
    '''Plot two kpaths against each other'''
    fig = go.Figure(data=[go.Scatter3d(x=kpath1[:, 0], y=kpath1[:, 1], z=kpath1[:, 2], mode='markers', marker=dict(color='red', size=5))])
    fig.add_trace(go.Scatter3d(x=kpath2[:, 0], y=kpath2[:, 1], z=kpath2[:, 2], mode='markers', marker=dict(color='blue', size=5)))
    fig.update_traces(connectgaps=False)
    fig.update_layout(title='Kpaths', title_x=0.5)
    #add legend for colors
    fig.update_layout(legend=dict( yanchor="top", y=0.99, xanchor="left", x=0.01))

    #add axis labels x= kx, y=ky, z=kz
    fig.update_layout(scene = dict( xaxis_title='kx', yaxis_title='ky', zaxis_title='kz'))
    fig.show()


def plot_energy_isosurface(es: ElectronicStructure, isoval: float = 0.0, tolerance: float = None):
    '''Plot the energy isosurface at a given energy value'''

    #get the kpoint coordinates
    path = es.kpoint_coordinates
    #get the kpoint energies
    energies: pd.DataFrame = es.kpoint_energies

    if tolerance is None:
        #get the average difference between kpoint energies
        tolerance = np.mean(np.diff(energies['energy'].to_numpy()))
    
    #plot x, y, and z coordinates of kpoints with energy rougly equal to isoval as a mesh grid
    fig = go.Figure(data=[go.Mesh3d(x=path[:, 0], y=path[:, 1], z=path[:, 2]) ] )

    fig.update_layout(title=f'Energy isosurface at {isoval} eV', title_x=0.5)
    #make aspect ratio equal
    fig.update_layout(scene_aspectmode='cube')
    fig.show()
    

def plot_isoplane(es: ElectronicStructure, isoval: float = 0.0, tolerance: float = None):
    '''Plot the energy isoplane at a given energy value'''
    from kspace import isosurface_bands


    if tolerance is None:
        #get the average difference between kpoint energies
        tolerance = np.mean(np.diff(es.kpoint_energies))

    bands: list[pd.DataFrame] = isosurface_bands(es, isoval)

    #plot x and y coordinates as a 2d scatter plot (line mode) with color corresponding to energy, for each band
    fig = go.Figure()
    for band in bands:
        fig.add_trace(go.Scatter(x=band['x'], y=band['y'], mode='lines'))

    fig.update_layout(title=f'Energy isoplane at {isoval} eV', title_x=0.5)
    #make aspect ratio equal
    fig.update_layout(scene_aspectmode='cube')
    fig.show()


def plot_band_contour(electronics: ElectronicStructure, band=None, title: str = None):
    '''Plots the band structure as a contour plot'''
    kpoints = electronics.values

    if band == None:

        # get the band that is closest to the fermi energy
        band = kpoints[np.isclose(
            kpoints['energy'], electronics.fermi_energy, rtol=0.05)]
        band = band[band['band'] == band['band'].max()]

    else:
        band = electronics.band(band)

    # plot kx vs ky with energy as the color in a contour plot with contour smoothing
    import plotly.graph_objects as go

    kx = band['x']
    ky = band['y']
    energy = band['energy'] - electronics.fermi_energy

    if len(kx) == 0 or len(ky) == 0 or len(energy) == 0:
        print("Cannot plot an empty plot")
    else:
        fig = go.Figure(data=go.Contour(x=kx, y=ky, z=energy, colorscale='rdbu',
                        contours=dict(coloring='heatmap'), line_smoothing=1.2))
        # label kx and ky
        fig.update_xaxes(title_text=r'$k_x$')
        fig.update_yaxes(title_text='$k_y$')
        # label energy
        fig.update_layout(coloraxis_colorbar=dict(title=r'$E - E_f$ (eV)'))
        fig.show()



