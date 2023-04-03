import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np


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


def create_band_traces(vasprun: Vasprun, emin: float = None, emax: float = None):
    '''Create the band traces for the plotly figure.'''
    kpoints = vasprun.bands.kpoints
    kpoints['energy'] -= vasprun.fermi_energy
    kpoints = filter_bands(kpoints, emin, emax)

    bands = kpoints.groupby('band').filter(
        lambda x: emin < x['energy'].mean() < emax)

    return bands


def add_fermi_line(fig):
    '''Adds a horizontal line at the fermi energy.'''
    fig.add_hline(y=0, line_width=0.5, line_color="grey")


def add_kpoint_divisions(fig, vasprun):
    '''Adds vertical lines at the kpoint divisions to seperate paths.'''
    divisions = vasprun.kpoint_linemode_divisions
    for i in range(1, len(vasprun.kpoints) + 1, divisions):
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


def plot_bandstructure(vasprun: Vasprun, emin: float = None,  emax=None, labels=None, show=True, legend: bool = True) -> px.line:
    '''Plot the bandstructure.'''
    bands = create_band_traces(vasprun, emin, emax)

    trace = go.Scatter(x=bands['kpoint'], y=bands['energy'], line_color='blue',
                       line_width=2, line_shape='spline', line_smoothing=0.5)
    fig = px.line(bands, x='kpoint', y='energy',
                  color_discrete_sequence=['blue'], color='band')

    add_fermi_line(fig)
    divisions = add_kpoint_divisions(fig, vasprun)

    if labels != None:
        kpoint_tick_val = np.arange(
            0, len(vasprun.kpoints) + 1, divisions) + 0.5
        fig.update_layout(xaxis=dict(
            tickmode='array', tickvals=kpoint_tick_val, ticktext=labels, tickfont=dict(size=18)))

    customize_bandstructure_layout(fig)

    if show:
        fig.show()

    return trace


def compare_bands(vasprun1: Vasprun, vasprun2: Vasprun, emin: float = -2.0, emax: float = 2.0, show: bool = True, save: bool = False, labels: list[str] = None, title: str = None, legend: bool = True):
    '''Plots two bandstructures on the same plot in different colors.'''
    bands1 = create_band_traces(vasprun1, emin, emax)
    bands2 = create_band_traces(vasprun2, emin, emax)

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
    fig.update_layout(plot_bgcolor='white')

    # make paper white
    fig.update_layout(paper_bgcolor='white')

    # make figure axes black lines
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=1, linecolor='black')
    # add vertical lines at every n*kpoint_linemode_divisions kpoint
    divisions = vasprun1.kpoint_linemode_divisions
    vlines = [i - 0.5 for i in range(1, len(vasprun1.kpoints) + 1, divisions)]
    fig.add_vline(x=vlines, line_width=0.5, line_color="black")

    if labels is not None:
        # create array of kpoint vals
        kpoint_tick_val = np.arange(
            0, len(vasprun1.kpoints) + 1, divisions) + 0.5
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


def plot_bsdos(vasprun: Vasprun, emin: float, emax: float, show: bool = True, save: bool = False, title: str = None):
    '''Plot the bandstructure and DOS.'''
    band_trace: px.line = plot_bandstructure(vasprun, emin=emin, emax=emax, show=False)
    dos_trace = plot_total_dos(
        vasprun, emin=emin - 2, emax=emax + 2, show=False)

    bsdos_fig = make_subplots(
        rows=1, cols=2, shared_yaxes=True, shared_xaxes=False)

    bsdos_fig.add_traces([band_trace, dos_trace], rows=[1, 1], cols=[1, 2])

    bsdos_fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    bsdos_fig.update_xaxes(
        title_text=r'$\textbf{k}$', showline=True, linewidth=1, linecolor='black', row=1, col=1)
    bsdos_fig.update_yaxes(
        title_text=r'$E-E_f \quad \mbox{[eV]}$', showline=True, linewidth=1, linecolor='black', row=1, col=1)

    bsdos_fig.update_xaxes(
        title_text=r'$\textbf{DOS}$', showline=True, linewidth=1, linecolor='black', row=1, col=2)
    bsdos_fig.update_yaxes(
        title_text=r'$E-E_f \quad \mbox{[eV]}$', row=1, col=2)

    if show:
        bsdos_fig.show()

    if save:
        bsdos_fig.write_image('bsdos.png')


y_axis = r'$E - E_f$ [eV]'
x_axis = r'$DOS$ [states/eV]'


def plot_total_dos(vasprun: Vasprun, emin=None, emax=None, title=None, show=True) -> go.Scatter:
    '''Plot the total density of states.'''
    data = vasprun.dos.total
    # add all spin channels together
    data = data.groupby('energy').sum().reset_index()

    # drop spin column
    data = data.drop(columns=['spin'])

    # subtract the fermi energy
    data['energy'] = data['energy'] - vasprun.fermi_energy

    # normalize the DOS
    data['total'] = data['total'] / data['total'].max()

    # constrain the energy range
    if emin != None:
        data = data[data['energy'] > emin]
    if emax != None:
        data = data[data['energy'] < emax]

    # sort data by energy
    data = data.sort_values(by='energy')

    # plot energy as the y axis and total as the x axis, add shaded area that is the total DOS
    trace = go.Scatter(x=data['total'], y=data['energy'], fill='tozerox',
                       line_color='blue', line_width=2, line_shape='spline', line_smoothing=0.5)
    fig = go.Figure(trace)

    # make line blue and thicker
    fig.update_traces(line_color='blue', line_width=2, overwrite=False)

    # add a horizontal line at the fermi energy
    fig.add_hline(y=0, line_width=0.5, line_color="grey")

    # set the axis labels
    fig.update_xaxes(title_text=x_axis)
    fig.update_yaxes(title_text=y_axis)

    # set the title
    if title == None:
        title = f'Total DOS for {vasprun.formula}'
    fig.update_layout(title=title, title_x=0.5, font_family='Times New Roman')

    if show:
        fig.show()

    return trace 

def plot_partial_dos(vasprun: Vasprun, ion: int, is_lm=True, title=None, dark_mode=False):
    '''Plot the partial density of states for a given ion.'''
    import plotly.express as px
    partial_dos = vasprun.dos.partial
    partial_dos = partial_dos[partial_dos['ion'] == ion]
    partial_dos['energy'] = partial_dos['energy'] - vasprun.fermi_energy
    # plot all orbitals of the ion
    x = [label for label in partial_dos.columns if label not in [
        'energy', 'ion', 'spin']]

    if is_lm == False:
        # combine the orbitals into l
        partial_dos['s'] = partial_dos['s']
        partial_dos['p'] = partial_dos['py'] + \
            partial_dos['pz'] + partial_dos['px']
        partial_dos['d'] = partial_dos['dxy'] + partial_dos['dyz'] + \
            partial_dos['dz2'] + partial_dos['dxz'] + partial_dos['x2-y2']
        x = ['s', 'p', 'd']

    fig = px.line(partial_dos, x=x, y='energy')

    # set the axis labels
    fig.update_xaxes(title_text=x_axis)
    fig.update_yaxes(title_text=y_axis)

    # set the title
    if title == None:
        title = f'Partial DOS of {vasprun.atom_types[ion]}{ion}'

    if dark_mode == True:
        fig.update_layout(template='plotly_dark')

    fig.update_layout(title=title, title_x=0.5, font_family='Times New Roman')

    # set the legend
    fig.update_layout(legend_title_text='Orbital')

    # allow autoscaling to fit the window
    fig.update_layout(autosize=True)
    fig.show()