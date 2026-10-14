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
import networkx as nx 
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
                     title_font_size=20,
                     tickfont_size = 18,
                     ticks='inside'
                    )
    
    fig.update_yaxes(mirror=True,ticks='inside',linecolor='black', zeroline=False,
                     title_font_size=20,
                     tickfont_size = 20,
                    )
    if xaxis_range != None:
        fig.update_xaxes(range=xaxis_range)
    if yaxis_range != None:
        fig.update_yaxes(range=yaxis_range)
    return fig

plt.rcParams["font.family"] = "Avenir"
color = cm.rainbow(np.linspace(0, 1, 12))
colors = [rgb2hex(c) for c in color]
'''
WHAT TO EDIT**:

------------------------------------------------------

'''
# our signac project
project = signac.get_project()
# look through each POLYMERIZED job in the project

# group_statepoint = {"replica_index": 8, "density": 0.9, "temperature": 0.9, "crosslinker_percent": 50, 
#                     "N_monomers": 500, "monomer_size": 0, "extender_size": 0, "radical_percent": 1.0, 
#                     "chain_side_reaction_probability": 0, "chain_transfer_probability": 0.5, 
#                     "thiol_reaction_probability": 0.5, "polymerize_period": 100, "r_cut_reaction": 1.1, 
#                     "angle_constant": 5.0}
# # which group to plot
# group_plot_name = "0_mon_500_cg"
# group_plot_dir = "./plots/" + group_plot_name + "/"
# if not os.path.isdir(group_plot_dir):
#     os.mkdir(group_plot_dir)
# cg_only = True


'''
END**
------------------------------------------------------
'''

polymerization_methods = ['custom_action_GPU','custom_action_GPU_bulk']
TPS = [[],[]]
avg_polymerization_time = [[],[]]
avg_integration_time = [[],[]]
integration_tps = [[],[]]

extent_of_reaction = [[],[]]

for i,polymerization_method in enumerate(polymerization_methods):
    for job in project:
        print(job.id)
        # if the job has been polymerized to any extent
        if job.isfile('polymerize.gsd'):
            with open(job.fn('signac_statepoint.json')) as f:
                statepoint = json.load(f)
            with open(job.fn('signac_job_document.json')) as f:
                job_document = json.load(f)
            print(polymerization_method)
            print(statepoint["polymerization_method"])
            if statepoint["polymerization_method"] != polymerization_method:
                continue
            
            try:
                TPS[i].append(job_document["TPS"])
                avg_polymerization_time[i].append(job_document["avg_polymerization_time"])
                avg_integration_time[i].append(job_document["avg_integration_time"])
                integration_tps[i].append(job_document["integration_tps"])
                extent_of_reaction[i].append(job_document["reacted_monomers"])
                # print(avg_integration_time)
            except:
                print("file exist but 0 size", job.id)

fig, ax = plt.subplots(1,5)

width = 0.25
r = np.arange(5)
for i in range(len(polymerization_methods)):

    avg_TPS = np.mean(TPS[i])
    avg_avg_polymerization_time = np.mean(avg_polymerization_time[i])
    avg_avg_integration_time = np.mean(avg_integration_time[i])
    avg_integration_tps = np.mean(integration_tps[i])
    avg_extent_of_reaction = np.mean(extent_of_reaction[i])

    data = [avg_TPS,avg_integration_tps,avg_avg_polymerization_time,avg_avg_integration_time,avg_extent_of_reaction]
    error = [np.std(TPS[i]),np.std(integration_tps[i]),np.std(avg_polymerization_time[i]),np.std(avg_integration_time[i]),np.std(extent_of_reaction[i])]

    shift = i*width
    # plt.bar(r+shift,data,width=width,label=polymerization_methods[i])
    # plt.errorbar(r+shift, data, yerr=error, fmt="o", color="r")
    for j in range(5):
        ax[j].bar(0+shift,data[j],width=width,label=polymerization_methods[i])
        ax[j].errorbar(0+shift,data[j], yerr=error[j], fmt="o", color="r")    

ax[0].set_title("TPS",wrap=True)
ax[1].set_title("Integration TPS",wrap=True)
ax[2].set_title("Polymerization time",wrap=True)
ax[3].set_title("Integration time",wrap=True)
ax[4].set_title("Extent of Reaction",wrap=True)
for axis in ax:
    axis.set_xticks([])

# plt.ylabel("Number of people voted") 
# plt.title("Number of people voted in each year") 
  
# plt.grid(linestyle='--') 
# plt.xticks(r + width/2,['TPS','Integration TPS','Polymerization time','Integration time','Extent of Reaction']) 
plt.legend()
plt.tight_layout()
  
plt.show()
    
