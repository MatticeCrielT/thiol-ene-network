import sys, os, re, json, io, itertools
import numpy as np
import signac
import gsd, gsd.hoomd 
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from datetime import datetime

import networkx as nx

def dangling_end_analysis(job_id):
    project = signac.get_project()
    direc = project.fn('') + 'workspace/'

    # for job_id in os.listdir(direc):
    jdir = direc + job_id
    input_file = direc + job_id + "/polymerize.gsd"
    if not os.path.isfile(input_file):
        raise FileNotFoundError("File not found")
        exit()
    # open the polymerized gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)

    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
    # make a graph - each bond is an edge
    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    # calculate gel fraction
    sizes = [len(n) for n in Gcc ]
    total = sum(sizes)
    sizes = sorted(Counter(sizes).items(), key=lambda item: item[0], reverse=True)

    new_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                new_bonds.append(each)
                continue
   
   # make cluster of biggest graph
    G2 = nx.Graph()
    G2.add_edges_from(new_bonds)

    # make a fake G3 for the first run through
    G3 = nx.Graph()
    G3.add_edges_from([[-1,1]])
    first_run_q == 0
    while G3.nodes != G2.nodes:
        if first_run_q != 0:
            G2 = G3.copy()
        else:
            first_run_q == 1
        # find all one degree beads in the gel
        one_degree_beads = [x for  x in G2.nodes() if G2.degree(x) == 1]
        # find all crosslinker beads in the gel
        crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) >= 3]
        G_new = G2.copy()
        # remove all crosslink beads
        G_new.remove_nodes_from(crosslink_beads)
        strands = nx.connected_components(G_new)
        # remove strands with only one end attached to network and store in new copy
        remove = [i for i in strands if set(strands).intersection(set(one_degree_beads))]
        remove = np.array(remove).flatten()
        G3 = G2.copy()
        G3.remove_nodes_from(remove)
    
    Gcc_new = sorted(nx.connected_components(Gnew), key=len, reverse=True)
    print("number of strands network",len(Gcc_new))
    sizes = [len(n) for n in Gcc_new]
    print("length of strands in network")
    sizes_of_strands,count = np.unique(sizes,return_counts=True)
    print(sizes_of_strands,count)    

