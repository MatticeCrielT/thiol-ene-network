import gsd, gsd.hoomd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os.path

trajectory= gsd.hoomd.open('polymerize.gsd')
frame = trajectory[-1]
bonds = frame.bonds.group

# make a graph - each bond is an edge
G = nx.Graph()
G.add_edges_from(bonds)

Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
Gi = G.subgraph(Gcc[0])

# since determining positions is very slow, this is an attempt of saving the layout to file
if os.path.exists('positions.pickle'):
    # load data
    with open('positions.pickle', 'rb') as f:
        pos = pickle.load(f)
else:
    pos = nx.kamada_kawai_layout(Gi)
    # save data
    with open('positions.pickle', 'wb') as f:
        pickle.dump(pos, f)


fig = plt.figure()
fig.add_subplot(111)

# sometimes, if the network is very heterogeneous this shows something interesting...
# there are a whole bunch of other algorithms, look at networkx documentation
# https://networkx.org/documentation/stable/reference/algorithms/community.html
lpc = nx.community.louvain_communities(Gi)
community_index = {n: i for i, com in enumerate(lpc) for n in com}
node_color = [community_index[n] for n in Gi]

nx.draw(Gi,pos,with_labels=True,node_size=10)

#nx.draw_networkx_nodes(Gi, pos,node_size=10,node_color=node_color,cmap=plt.cm.jet)

# thick edge = "important" edge - ie if this edge would be removed, the network connectivity
# would change a lot (read documentation, it's much more complicated than that... but..)
edge_centrality = nx.edge_betweenness_centrality(Gi)
#edge_centrality2 = nx.edge_current_flow_betweenness_centrality(Gi)
edge_weight= [v*100 for v in edge_centrality.values()]

nx.draw_networkx_edges(Gi, pos,width=edge_weight)


plt.show()


