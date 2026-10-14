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

"""
frame.bonds.group gives a np.array that is of form [[a,b],[c,d],[e,f],[a,g],...] 
where a is bonded to b, c is bonded to d, etc.

atom_bonds finds all atoms to which a given atom 'a' is bonded to
args: frame.bonds.group

returns: dictionary of atoms and the atoms they are bonded to
{atom_no: [bonded1, bonded2, ...], ...}
"""
def atom_bonds(bonds):
    # create a default dictionary with an element returning a list if no neighbors
    neighbors = defaultdict(list) 
    # add every instance of an atom in a bond to the dictionary
    for each in bonds:
        for item in each:
            neighbors[item].extend(each)
            # remove the atom itself from the list of atoms bonded to it
            neighbors[item].remove(item)
    return neighbors

### a-b-c
###  \ /
###   V
def cyclopropyl_groups(a_bonds):
    cyclopropyl_count = 0
    # get keys for atoms with more than one bond
    keys = [i for i in a_bonds.keys() if len(a_bonds[i])>1]
    cyclo_keys = []
    cyclo_keys_flat = []
    # a
    for a in keys:
        if a in cyclo_keys_flat:
            continue
        # b
        for b in a_bonds[a]:
            if b in cyclo_keys_flat:
                continue
            # c
            for c in a_bonds[b]:
                if c in cyclo_keys_flat:
                    continue
                if a in a_bonds[c]:
                    cyclo_keys = cyclo_keys + [[a,b,c]]
                    cyclopropyl_count += 1
                    cyclo_keys_flat = [x for xs in cyclo_keys for x in xs]
    return(cyclopropyl_count, cyclo_keys)

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
color = cm.rainbow(np.linspace(0, 1, 11))
colors = [rgb2hex(c) for c in color]
'''
WHAT TO EDIT**:

------------------------------------------------------

'''
# our signac project
project = signac.get_project()
# look through each POLYMERIZED job in the project

group_statepoint = {"replica_index": 8, "density": 0.9, "temperature": 0.9, "crosslinker_percent": 50, 
                    "N_monomers": 500, "monomer_size": 0, "extender_size": 0, "radical_percent": 1.0, 
                    "chain_side_reaction_probability": 0, "chain_transfer_probability": 0.5, 
                    "thiol_reaction_probability": 0.5, "polymerize_period": 100, "r_cut_reaction": 1.1, 
                    "angle_constant": 5.0}
# which group to plot
group_plot_name = "0_mon_500_cg"
group_plot_dir = "./plots/" + group_plot_name + "/"
if not os.path.isdir(group_plot_dir):
    os.mkdir(group_plot_dir)
cg_only = True


'''
END**
------------------------------------------------------
'''

cyclopropyl_nums = []
for job in project: 
    # if the job has been polymerized to any extent
    if job.isfile('polymerize.gsd'):
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)
            ### only analyze the plot group:
            if (statepoint["chain_side_reaction_probability"]==0 or statepoint["chain_transfer_probability"]!=0):
                continue
        
        try:
            ''' initial look '''
            # get and trajectory and make sure the first frame checks out
            trajectory = gsd.hoomd.open(job.fn('polymerize.gsd'))
            frame = trajectory[0]
            # get the dummy_ids to avoid them
            dummy_id = len(frame.bonds.types)-1
            bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id] # dummy type bond 
            # get all the ids and bonds between them
            nids, counts = np.unique(np.concatenate(bonds).flatten(),return_counts=True)

            ''' analyze '''
            frame = trajectory[-1]
            bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
            count, cp_list = cyclopropyl_groups(atom_bonds(bonds))
            cyclopropyl_nums.append(count)
        except:
            print("file exist but 0 size", job.id)

# molecule size histogram at final time point
cycloprop_fig = make_plotly_fig()
cycloprop_fig.update_layout(bargap=0)
cycloprop_fig.update_xaxes(tickvals = [0,1,2,3,4,5,6,7,8,9])
cycloprop_fig.add_bar(name="beep", x=list(range(0,10,1)), y=cyclopropyl_nums,marker_color=colors)
cycloprop_fig.write_image(group_plot_dir + "cyclopropyl_groups.pdf")