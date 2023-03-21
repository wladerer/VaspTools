from vasprun import Vasprun
import pandas as pd
import plotly.express as px
#use simple_white theme
px.defaults.template = 'simple_white'


y_axis = r'$E - E_f$ [eV]'
x_axis = r'$DOS$ [states/eV]'



def plot_total_dos(vasprun: Vasprun, is_spin=False, title=None, dark_mode=False):
    '''Plot the total density of states.'''
    data = vasprun.dos.total
    #subtract the fermi energy
    data['energy'] = data['energy'] - vasprun.fermi_energy
    if is_spin == False:
        fig = px.line(data, y='energy', x ='total')
    else:
        #plot spin ==1 and spin == 2
        fig = px.line(data['spin' == 1], y='energy', x='total')
        fig = px.line(data['spin' == 2], y='energy', x='total')

    #set the axis labels
    fig.update_xaxes(title_text=x_axis)
    fig.update_yaxes(title_text=y_axis)

    if dark_mode == True:
        fig.update_layout(template='plotly_dark')

    #set the title
    if title == None:
        title = f'Total DOS for {vasprun.formula}'
    fig.update_layout(title=title, title_x=0.5, font_family='Times New Roman')
    

    fig.show()

def plot_partial_dos(vasprun: Vasprun, ion: int, is_lm=True, title=None, dark_mode=False):
    '''Plot the partial density of states for a given ion.'''

    partial_dos = vasprun.dos.partial
    partial_dos = partial_dos[partial_dos['ion'] == ion]
    partial_dos['energy'] = partial_dos['energy'] - vasprun.fermi_energy
    #plot all orbitals of the ion
    x = [ label for label in partial_dos.columns if label not in ['energy', 'ion', 'spin']]

    if is_lm == False:
        #combine the orbitals into l 
        partial_dos['s'] = partial_dos['s']
        partial_dos['p'] = partial_dos['py'] + partial_dos['pz'] + partial_dos['px']
        partial_dos['d'] = partial_dos['dxy'] + partial_dos['dyz'] + partial_dos['dz2'] + partial_dos['dxz'] + partial_dos['x2-y2']
        x = ['s', 'p', 'd']


    fig = px.line(partial_dos, x=x, y='energy')

    #set the axis labels
    fig.update_xaxes(title_text=x_axis)
    fig.update_yaxes(title_text=y_axis)

    #set the title
    if title == None:
        title = f'Partial DOS of {vasprun.atom_types[ion]}{ion}'

    if dark_mode == True:
        fig.update_layout(template='plotly_dark')

    fig.update_layout(title=title, title_x=0.5, font_family='Times New Roman')

    #set the legend
    fig.update_layout(legend_title_text='Orbital')

    #allow autoscaling to fit the window
    fig.update_layout(autosize=True)
    fig.show()



