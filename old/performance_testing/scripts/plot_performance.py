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
                      height=500,
                      xaxis_title=xaxis,
                      yaxis_title=yaxis)
    
    fig.update_xaxes(linecolor='black', linewidth=2,mirror=True,
                     zeroline=False,
                     title_font_size=25,tickwidth=2,tickcolor='black',
                     tickfont_size = 25,
                     ticks='inside'
                    )
    
    fig.update_yaxes(linewidth=2,mirror=True,ticks='inside',linecolor='black', zeroline=False,
                     title_font_size=25,tickwidth=2,tickcolor='black',
                     tickfont_size = 25,
                    )
    if xaxis_range != None:
        fig.update_xaxes(range=xaxis_range)
    if yaxis_range != None:
        fig.update_yaxes(range=yaxis_range)
    return fig

# plt.rcParams["font.family"] = "Avenir"
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

n_replicas = 2
polymerization_methods = ['cpu_local_snapshot','custom_action_GPU','custom_action_CPU']
polymerization_periods = [1,10,100,1000,10000,100000]
TPS = np.zeros(shape=(3,6))
avg_polymerization_time = np.zeros(shape=(3,6))
avg_integration_time = np.zeros(shape=(3,6))
integration_tps = np.zeros(shape=(3,6))

extent_of_reaction = np.zeros(shape=(3,6))

for job in project:
    print(job.id)
    # if the job has been polymerized to any extent
    if job.isfile('polymerize.gsd'):
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)
        with open(job.fn('signac_job_document.json')) as f:
            job_document = json.load(f)
            print(job_document)
        
        i = polymerization_methods.index(statepoint["polymerization_method"])
        j = polymerization_periods.index(statepoint["polymerize_period"])

        # try:
        # TPS[i][j] += job_document["TPS"]/n_replicas
        avg_polymerization_time[i][j] += float(job_document["avg_polymerization_time"])/n_replicas
        avg_integration_time[i][j] += float(job_document["avg_integration_time"])/n_replicas
        integration_tps[i][j] += float(job_document["integration_tps"])/n_replicas
        extent_of_reaction[i][j] += float(job_document["reacted_monomers"])/n_replicas
        if float(job_document["avg_integration_time"]) > 1000000:
            print(job.id)
            exit()
            # TPS[i].append(job_document["TPS"])
            # avg_polymerization_time[i].append(job_document["avg_polymerization_time"])
            # avg_integration_time[i].append(job_document["avg_integration_time"])
            # integration_tps[i].append(job_document["integration_tps"])
            # extent_of_reaction[i].append(job_document["reacted_monomers"])
            # print(avg_integration_time)
        # except:
        #     print("file exist but 0 size", job.id)

        #Get the TPS from the final timestep:
        data = np.genfromtxt(job.fn("polymerize.txt"),dtype=float, skip_header=1)
            # params = input_arr[i].split('_')

        # Get the column headers from the first line of the txt file
        with open(job.fn("polymerize.txt")) as f:
            line = f.readline()
            line = line.replace('\n','')
        headers = line.split()
        print(headers)
        #Some formatting to make the headers look nice when used for plotting,
        # and also removing any preamble that hoomd 4 adds to its log files
        headers[0] = headers[0].replace('#','')
        for k,header in enumerate(headers):
            # headers[i] = headers[i].replace('_',' ')
            headers[k] = headers[k].replace('md.compute.ThermodynamicQuantities.','')
            headers[k] = headers[k].replace('Simulation.','')
        print("headers:", headers)

        final_timestep = data[-1][0]
        print(final_timestep)
        TPS[i][j] += final_timestep/(60*60*2)/n_replicas

fig, ax = plt.subplots(3,2,sharex=True)
colors = plt.get_cmap('viridis')(np.linspace(0,1,num=3))
colors = [rgb2hex(c) for c in colors]

width = 0.25
r = np.arange(5)
for i in range(len(polymerization_methods)):
    '''
    ax[0][0].plot(np.log10(polymerization_periods),TPS[i],color=colors[i],label=polymerization_methods[i],linestyle='--',marker='o')
    ax[1][0].plot(np.log10(polymerization_periods),integration_tps[i],color=colors[i],label=polymerization_methods[i],linestyle='--',marker='o')
    ax[2][0].plot(np.log10(polymerization_periods),avg_polymerization_time[i],color=colors[i],label=polymerization_methods[i],linestyle='--',marker='o')
    ax[0][1].plot(np.log10(polymerization_periods),avg_integration_time[i],color=colors[i],label=polymerization_methods[i],linestyle='--',marker='o')
    ax[1][1].plot(np.log10(polymerization_periods),extent_of_reaction[i],color=colors[i],label=polymerization_methods[i],linestyle='--',marker='o')
    '''
    figure = make_plotly_fig(xaxis="log10 Polymerization Period", yaxis="Timesteps per second")
    figure.add_trace(go.Scatter(x=np.log10(polymerization_periods), y=TPS[i], 
                                 mode='lines+markers', name=polymerization_methods[i],
                                 line=dict(color=colors[i], width=2, dash='dash')))
    figure.show()
    exit()

    # # plt.bar(r+shift,data,width=width,label=polymerization_methods[i])
    # # plt.errorbar(r+shift, data, yerr=error, fmt="o", color="r")
    # for j in range(5):
    #     ax[j].bar(0+shift,data[j],width=width,label=polymerization_methods[i])
    #     ax[j].errorbar(0+shift,data[j], yerr=error[j], fmt="o", color="r")    

ax[0][0].set_title("TPS",wrap=True)
ax[0][0].set_ylabel("Timesteps per second")
ax[1][0].set_title("Integration TPS",wrap=True)
ax[1][0].set_ylabel("Timesteps per second during simulation")
ax[2][0].set_title("Polymerization time",wrap=True)
ax[2][0].set_ylabel("Average time (s)")
ax[0][1].set_title("Integration time",wrap=True)
ax[0][1].set_ylabel("Average time (s)")
ax[1][1].set_title("Extent of Reaction",wrap=True)
ax[1][1].set_ylabel("")
# for axis in ax:
    # axis.set_xticks(polymerization_periods)
    # axis.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# ax[0].legend(draggable=True)

ax[2][0].set_xlabel(r"$log_{10}$ Polymerization period")
ax[2][1].set_xlabel(r"$log_{10}$ Polymerization period")
# plt.title("Number of people voted in each year") 
  
# plt.grid(linestyle='--') 
# plt.xticks(r + width/2,['TPS','Integration TPS','Polymerization time','Integration time','Extent of Reaction']) 
#plt.legend()
ax[0][1].legend(draggable=True)
#fig.legend()
plt.tight_layout()

plt.show()
    
