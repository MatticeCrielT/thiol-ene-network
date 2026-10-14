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
fig = make_plotly_fig(xaxis="Crosslinker", xaxis_range=[-2,102], yaxis = "Number Dangling Ends")


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

        # gel fractions
        df = pd.read_table(job.fn('cluster_sizes.txt'), sep=" ")
        df.dropna(inplace=True, axis=1)
        df.columns = ["cluster_size", "count"]
        #percolation_dimension_set.append(list(df["largest_cluster"]))
        #percolation_dimension_frames.append(list(df["frame"]))
        gf = float(df["cluster_size"].max() / sum(df["cluster_size"]*df["count"]))
        if cl==87.5 and gf < 0.1:
            continue

        df = pd.read_table(job.fn('dangle_sizes.txt'), sep=" ")
        df.columns = ["cluster_size", "count"]
        first.append(sum(df["count"]))

        
        fig2 = make_plotly_fig(xaxis="Dangle Size", yaxis = "Number")
        fig2.add_trace(go.Bar(x=[str(i) for i in df["cluster_size"]], y=df["count"],marker= dict(color = c, line=dict(color=c))))
    
        df = pd.read_table(job.fn('strand_sizes.txt'), sep=" ")
        df.dropna(inplace=True, axis=1)
        df.columns = ["strand_sizes", "count"]

    fig2.write_image("./plots/dangling_ends/TE_" + str(cl) + ".pdf")
    #fig2.show()
    
    print(first)
    first_std = np.array(first).std(axis=0).tolist()
    
    fig.add_trace(go.Scatter(
        x=[cl],
        y=[sum(first)/len(first)],
        marker = dict(color=c, size=8),
        error_y=dict(
            type='data',
            symmetric=True,
            array=[first_std])
        ))
    

fig.write_image("./plots/dangling_TE.pdf")
fig.show()
    


