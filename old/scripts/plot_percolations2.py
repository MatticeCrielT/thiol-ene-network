import sys, os, re, json, io, itertools
import numpy as np
import pandas as pd
import signac
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.pyplot import cm
from matplotlib.colors import rgb2hex
import gsd, gsd.hoomd 
from collections import defaultdict
from collections import OrderedDict
from collections import Counter

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import seaborn as sns
import plotly.io as pio   
pio.kaleido.scope.mathjax = None


# make a default figure box 
def make_plotly_fig(xaxis="", xaxis_range=None, yaxis="", yaxis_range=None):
    fig = go.Figure()
    fig.update_layout(
        font_family="Avenir Medium",
        font_color="black",
    )
    fig.update_layout(
        plot_bgcolor='white'
    )
    fig.update_layout(showlegend=False,
                      width=500, 
                      height =500,
                      xaxis_title=xaxis,
                      yaxis_title=yaxis)
    
    fig.update_xaxes(linecolor='black', mirror=True,
                     zeroline=False,
                     title_font_size=25,
                     tickfont_size = 25,
                     ticks='inside'
                    )
    
    fig.update_yaxes(mirror=True,ticks='inside',linecolor='black', zeroline=False,
                     title_font_size=25,
                     tickfont_size = 25,
                    )
    if xaxis_range != None:
        fig.update_xaxes(range=xaxis_range)
    if yaxis_range != None:
        fig.update_yaxes(range=yaxis_range)
    return fig

# our signac project
project = signac.get_project()
# look through each POLYMERIZED job in the project
csr_p =  0
crosslinkers = [0,12.5,25,37.5,50,62.5,75,87.5,100]
color = iter(cm.rainbow(np.linspace(0, 1, 11)))
fig = make_plotly_fig(xaxis="Polymerization Cycle", yaxis = "Percolation Dimension")
fig_2 = go.Figure()

fig_2.update_layout(
        font_family="Avenir Medium",
        font_color="black",
        plot_bgcolor='white',
        legend_font_size=15,
        xaxis_title="Polymerization Cycle")

fig_2.update_xaxes(range=[0,1500],showgrid=False, tickfont_size = 25, title_font_size=25)
fig_2.update_yaxes(showgrid=False, 
                 zeroline=True, zerolinecolor='black', zerolinewidth=3,
                 showticklabels=False)
fig_2.update_layout(height=200, width = 500, plot_bgcolor='white', showlegend = False)

for cl in crosslinkers:
    c = next(color)
    c = rgb2hex(c)
    percolation_dimension_set = []
    percolation_dimension_frames = []
    first = []
    # collect all replicas
    for job in project: 
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)
            ### only analyze the plot group:
            if (statepoint["chain_side_reaction_probability"]!=csr_p):
                continue
            if (statepoint["crosslinker_percent"]!=cl):
                continue
        # if the job has been deformed
        if not job.isfile('percolation_data.txt'):
            continue
        ''' initial look '''
        df = pd.read_table(job.fn('percolation_data.txt'), sep=" ")
        df.dropna(inplace=True, axis=1)
        df.columns = ["frame", "percolation_dimension", "largest_cluster", "molecule_num"]
        #percolation_dimension_set.append(list(df["largest_cluster"]))
        #percolation_dimension_frames.append(list(df["frame"]))
        first.append(df[df["percolation_dimension"]>2]["frame"].min())
        fig.add_trace(go.Scatter(
            x=np.array(range(0,len(df["frame"]))), y = df["largest_cluster"], mode = "markers", marker = dict(color=c, size=8, opacity=.6),
        ))
        break
    '''
    minimum_frames = min([len(pds) for pds in percolation_dimension_set])
    percolation_dimension_set = [pds[:minimum_frames] for pds in percolation_dimension_set]
    percolation_dimensions = np.array(percolation_dimension_set).mean(axis=0).tolist()
    percolation_dimensions_std = np.array(percolation_dimension_set).std(axis=0).tolist()

    
    
    print(sum(first)/len(first)))))'''
    fig_2.add_trace(go.Scatter(x=[sum(first)/len(first)], y=[0], marker = dict(color=c, size=12)))
fig.write_image("./plots/ThiolEnePercolation.pdf")
fig_2.write_image("./plots/TEPercolation_nline.pdf")
fig.show()
fig_2.show()

    
        
        


