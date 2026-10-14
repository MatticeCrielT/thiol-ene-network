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

def dist_pbc(x1,x0,Box):
    delta = np.abs(x1 - x0)
    delta= np.where(delta > 0.5 * Box, Box - delta, delta)
    return np.sqrt(np.sum(delta**2.0))

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


def connected_components(lists):
    R"""
    merges lists with common elements

    args: bonds from configuration (frame.bonds.group)

    returns: list of connected particles by id

    Useful for finding bonded particles in a configuration, tested for
    linear polymers with consecutive bonds (0-1-2-3-4-5, 6-7-8-9-10,..)
    and non consecutive ids ( 0-5-6-8-10, 1-4-3-2-9,...) but no other
    configuration yet. Works with ints as well as str.

    """
    neighbors = defaultdict(set)
    seen = set()
    for each in lists:
        for item in each:
            neighbors[item].update(each)
    def component(node, neighbors=neighbors, seen=seen, see=seen.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            see(node)
            nodes |= neighbors[node] - seen
            yield node
    for node in neighbors:
        if node not in seen:
            yield sorted(component(node))

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
color = iter(cm.rainbow(np.linspace(0, 1, 11)))

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

for job in project: 
    # if the job has been deformed
    if job.isfile('polymerize.gsd'):
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)

        # try:
        ''' initial look '''
        # get and trajectory and make sure the first frame checks out
        trajectory = gsd.hoomd.open(job.fn('polymerize.gsd'))
        print(len(trajectory))


        ''' analyze '''
        frames = []
        avg_diffusions = []

        # iterate through each frame
        for i,frame in enumerate(trajectory):
            if i < len(trajectory)-1:
                frames.append(i)
                pos = frame.particles.position
                next_pos = trajectory[i+1].particles.position
                diffusions = []
                for j in range(len(pos)):
                    diff = dist_pbc(pos[j],next_pos[j],frame.configuration.box[0:3])
                    diffusions.append(diff)
                avg_diffusions.append(np.average(diffusions))
        print("frames",frames)
        print("avg_diff:",avg_diffusions)
                    
            
        # ### jsons ------------------------------------------
        # json_dir = "./workspace/" + str(job) + "/json/"

        # if not os.path.isdir(json_dir):
        #     os.mkdir(json_dir)
        # # make a json for number of bonds on each atom through the trajectory
        # with open(json_dir + 'bond_histograms.json', "w") as f :
        #     json.dump(trajectory_bond_histograms, f)

        # # make a json for number of radicals through the trajectory
        # with open(json_dir + 'radical_numbers.json', "w") as f:
        #     json.dump(trajectory_radical_numbers, f)

        # # make a json for histogram of molecule sizes at last frame
        # molec_size_hist = dict()
        # for i in strand_lengths:
        #     molec_size_hist[i] = molec_size_hist.get(i, 0) + 1
        # with open(json_dir + 'molecules_histogram.json', "w") as f :
        #     json.dump(trajectory_molecule_sizes, f)
        
        # # make a json for histogram of undreacted atoms at last frame
        # with open(json_dir + 'non_reacted_histogram.json', "w") as f :
        #     json.dump(trajectory_unreacted_atoms, f)
        
        ### plots ------------------------------------------
        plot_dir = "./workspace/" + str(job) + "/plots/"
        if not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)

        fig, ax = plt.subplots(1,1,sharey=False)

        ax.plot(frames,
                avg_diffusions,
                linewidth=2)
                    # N = 5
        # ax.plot(np.convolve(unique_strains,np.ones(N)/N,mode='valid'),
        #         np.convolve(unique_true_stress_straight,np.ones(N)/N,mode='valid'),
        #         linewidth=2,label=r'$\sigma_{xx}$')
        # ax.plot(unique_strains,unique_true_stress_straight)
        # ax[1].plot(strain,true_stress_deviatoric)
        # plt.savefig(plot_dir + 'stress_strain.png')

        plt.savefig(plot_dir + 'stress_strain.png')

        # c = rgb2hex(c)
        # # molecule size histogram at final time point
        # mol_size_fig = make_plotly_fig()
        # mol_size_fig.update_layout(bargap=0)
        # molec_size_hist = OrderedDict(sorted(molec_size_hist.items()))
        # mol_size_fig.add_bar(name=i, x=[str(i) for i in list(molec_size_hist.keys())], y=list(molec_size_hist.values()),marker_color=c)
        # mol_size_fig.write_image(plot_dir + "molecule_size_histogram.pdf")
        
        # # unreacted atoms histogram at final time point
        # unreact_fig = make_plotly_fig()
        # unreact_fig.update_layout(bargap=0)
        # unreact_fig.add_bar(name=i,x=[str(i) for i in unreact_hist.keys()], y=list(unreact_hist.values()), marker_color=c)
        # unreact_fig.write_image(plot_dir + "unreacted_atoms_histogram.pdf")

        # # bond number histogram at final time point
        # bond_no_fig = make_plotly_fig()
        # bond_no_fig.update_layout(bargap=0)
        # bond_no_fig.add_bar(name=i,x=list(bond_histogram.keys()), y=list(bond_histogram.values()), marker_color=c)
        # bond_no_fig.write_image(plot_dir + "bond_number_histogram.pdf")

        # # radical consistency
        # radical_fig = make_plotly_fig(xaxis='Frames', yaxis='Number of radicals', yaxis_range=[7.75, 10.25])
        # radical_fig.update_yaxes(tickvals=[8,9,10])
        # radical_fig.add_hline(y=10, line_dash="dash", line_color="grey")
        # radical_fig.add_trace(
        #                       go.Scatter(x=list(trajectory_radical_numbers.keys()), 
        #                                  y = [i[0] for i in trajectory_radical_numbers.values()],
        #                                  mode='markers', 
        #                                  marker=dict(color=c, opacity=0.5)
        #                                  )
        #                       )
        # radical_fig.write_image(plot_dir + "radical_consistency.pdf")
        
        # # strands between crosslinks
        # if group_plot_dir != None:
        #     if not os.path.isdir(group_plot_dir + "moleculue_sizes/"):
        #         os.mkdir(group_plot_dir + "moleculue_sizes/")
        #         os.mkdir(group_plot_dir + "unreacted_atoms/")
        #         os.mkdir(group_plot_dir + "bond_numbers/")
        #         os.mkdir(group_plot_dir + "radical_consistency")
        #     mol_size_fig.write_image(group_plot_dir + "moleculue_sizes/" + "{}.pdf".format(job.id))
        #     unreact_fig.write_image(group_plot_dir + "unreacted_atoms/" + "{}.pdf".format(job.id))
        #     bond_no_fig.write_image(group_plot_dir + "bond_numbers/" + "{}.pdf".format(job.id))
        #     radical_fig.write_image(group_plot_dir + "radical_consistency/" + "{}.pdf".format(job.id))
        # except:
        #     print("file exist but 0 size", job.id)
plt.show()
        
        


