import plotly.graph_objects as go

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
