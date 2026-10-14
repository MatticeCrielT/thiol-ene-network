import gsd, gsd.hoomd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


trajectory= gsd.hoomd.open('polymerize.gsd')
frame = trajectory[-1]
bonds = frame.bonds.group
Ntotal = frame.particles.N

# make a graph - each bond is an edge
G = nx.Graph()
G.add_edges_from(bonds)


# find all disconnected/connected sub-networks
Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
print("number of disconnected clusters in network",len(Gcc))
sizes = [len(n) for n in Gcc ]
print("sizes of clusters in network",sizes)

# Gcc[0]=largest cluster = gel part of the system, everything else = sol


# find number of bonds ("degree") on each particle
degrees = [val for (node, val) in G.degree()]
num_bonds_on_node,count = np.unique(degrees,return_counts=True)
print("number of bonds on monomers")
print(num_bonds_on_node,count)

# strand lengths (in entrie system. If you only want the ones in the gel, you need to look at Gcc[0] only)
# make copy to modify
Gnew = G.copy()

# remove everyone that has 3 or more bonds on it, only leaving linear strands
to_be_removed = [x for  x in G.nodes() if G.degree(x) >= 3]
for x in to_be_removed:
    Gnew.remove_node(x)

Gcc_new = sorted(nx.connected_components(Gnew), key=len, reverse=True)
print("number of strands network",len(Gcc_new))
sizes = [len(n) for n in Gcc_new]
print("length of strands in network")
sizes_of_strands,count = np.unique(sizes,return_counts=True)
print(sizes_of_strands,count)


# loops/cycles - this takes some time
# only find loops shorter than 50 monomers - can increase this, but then it takes longer
# I am also not 100% sure this finds ALL loops, I don't think so...needs to be checked
# we might need to write our own version of the algorithm to get all loops
loops = list(nx.simple_cycles(G,length_bound=50))
lengths_loops = [len(l) for l in loops]
print(lengths_loops)

