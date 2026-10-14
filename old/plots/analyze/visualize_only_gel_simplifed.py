import gsd, gsd.hoomd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os.path

def simplifyGraph(G):
    ''' Loop over the graph until all nodes of degree 2 have been removed and their incident edges fused '''
    g = G.copy()
    while any(degree==2 for _, degree in g.degree):
        g0 = g.copy() #<- simply changing g itself would cause error `dictionary changed size during iteration`
        for node, degree in g.degree():
            if degree==2:
                edges = g0.edges(node)
                edges = list(edges.__iter__())
                try:
                    a0,b0 = edges[0]
                    a1,b1 = edges[1]
                    e0 = a0 if a0!=node else b0
                    e1 = a1 if a1!=node else b1
                    g0.remove_node(node)
                    g0.add_edge(e0, e1)
                except:
                    g0.remove_node(node)
        g = g0
    return g

trajectory= gsd.hoomd.open('polymerize.gsd')
frame = trajectory[-1]
bonds = frame.bonds.group

# make a graph - each bond is an edge
G = nx.Graph()
G.add_edges_from(bonds)

Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
Gi = G.subgraph(Gcc[0])

# let's try to simplify the gel a little bit. Each linear strand is condensed down to
# one edge, the length/10 is the weight/color of this edge
Gnew = simplifyGraph(Gi)

weights = []
for e in Gnew.edges:
    u = nx.shortest_path(G, source=e[0], target=e[1])
    weights.append(len(u)/10.)


fig = plt.figure()
fig.add_subplot(111)

pos = nx.kamada_kawai_layout(Gnew)

nx.draw_networkx_nodes(Gnew, pos,node_size=10)
nx.draw_networkx_edges(Gnew, pos,width=weights,edge_color=weights,edge_cmap=plt.cm.rainbow)

plt.show()

