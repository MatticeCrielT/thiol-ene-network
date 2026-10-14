import gsd, gsd.hoomd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


trajectory= gsd.hoomd.open('polymerize.gsd')
frame = trajectory[-1]
bonds = frame.bonds.group

# make a graph - each bond is an edge
G = nx.Graph()
G.add_edges_from(bonds)


Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
print("number of disconnected clusters in network",len(Gcc))


# only plot 9 largest clusters (or all if less than 9)
N = np.min([9,len(Gcc)])

# rendering
fig = plt.figure(1); plt.clf()
# compute a grid size that will fit all graphs on it (couple blanks likely)
nr = int(np.ceil(np.sqrt(N)))
fig, ax = plt.subplots(nr, nr, num=1)

for i in range(N):
    Gi = G.subgraph(Gcc[i])
    # this is very slow!
    pos = nx.kamada_kawai_layout(Gi)

    # compute index for the subplot, and set this subplot as current
    ix = np.unravel_index(i, ax.shape)
    plt.sca(ax[ix])
    nx.draw_networkx_nodes(Gi, pos, ax=ax[ix])
    nx.draw_networkx_edges(Gi, pos, ax=ax[ix])

    ax[ix].set_axis_off()

plt.show()


