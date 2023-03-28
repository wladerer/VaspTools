from vasprun import Vasprun
import pandas as pd
import plotly.graph_objects as go


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
