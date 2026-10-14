import sys, os, re, json, io, itertools
import numpy as np
import signac
import gsd, gsd.hoomd 
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from datetime import datetime
import matplotlib.pyplot as plt

import networkx as nx

"""
Calculates the distance between two points in a periodic box.
"""
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
    # create a default dictionary with an element returning a set if nothing in it
    neighbors = defaultdict(list) 
    # add every instance of an atom in a bond to the dictionary
    for each in bonds:
        for item in each:
            neighbors[item].extend(each)
            # remove the atom itself from the list of atoms bonded to it
            neighbors[item].remove(item)
    return neighbors

def draw_graph(G):
    pos = nx.planar_layout(G)
    labels = {}    
    for node in G.nodes():
        #set the node name as the key and the label as its value 
        labels[node] = node
    nx.draw_networkx_nodes(G, pos, node_color = 'r', node_size = 100, alpha = 1)
    nx.draw_networkx_labels(G, pos, labels, font_size = 10)
    ax = plt.gca()
    if G.is_multigraph():
        for e in G.edges:
            if e[0] != e[1]:
                ax.annotate("",
                            xy=pos[e[0]], xycoords='data',
                            xytext=pos[e[1]], textcoords='data',
                            arrowprops=dict(arrowstyle="-", color="0.5",
                                            shrinkA=5, shrinkB=5,
                                            patchA=None, patchB=None,
                                            connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
                                            ),
                                            ),
                            )
            else:
                ax.annotate("",
                            xy=pos[e[0]], xycoords='data',
                            xytext=pos[e[1]], textcoords='data',
                            arrowprops=dict(arrowstyle="-", color="0.5",
                                            shrinkA=5, shrinkB=5,
                                            patchA=None, patchB=None,
                                            connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
                                            ),
                                            ),
                            )
    else:
        for e in G.edges:
            ax.annotate("",
                        xy=pos[e[0]], xycoords='data',
                        xytext=pos[e[1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="0.5",
                                        shrinkA=5, shrinkB=5,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=0"
                                        ),
                        )
    plt.axis('off')

def directly_bonded(beads, bonds):
    direct_bonds = []
    for each in bonds:
        for n,bead in enumerate(each):
            if bead in beads:
                direct_bonds.append(each[1-n])
                break
    return direct_bonds

def classify_loop_type(cycle,bonds):
    n_multifunctional = 0
    for bead in cycle:
        n_bonds = 0
        for bond in bonds:
            if bead in bond:
                n_bonds += 1
        if n_bonds > 2:
            n_multifunctional += 1
    if n_multifunctional == 0:
        # return "primary"
        raise ValueError("Complete loop found. Are you running on the largest cluster?")
    elif n_multifunctional == 1:
        return "primary"
    elif n_multifunctional == 2:
        return "secondary"
    elif n_multifunctional == 3:
        return "tertiary"
    else:
        return "quaternary+"

def classify_loop_type_simple(cycle):
    #cycle is the list of nodes in the loop
    #for this classification to be correct, these nodes need to only include crosslinkers
    match len(cycle):
        case 1:
            return "primary"
        case 2:
            return "secondary"
        case 3:
            return "tertiary"
        case 4:
            return "quaternary+"
        case _:
            return "quaternary+"

def snap_molecule_indices(snap):
    """Find molecule index for each particle.

    Given a snapshot from a trajectory, compute clusters of bonded molecules
    and return an array of the molecule index of each particle.

    Parameters
    ----------
    snap : gsd.hoomd.Snapshot
        Trajectory snapshot.

    Returns
    -------
    numpy array (N_particles,)

    """
    system = freud.AABBQuery.from_system(snap)
    num_query_points = num_points = snap.particles.N
    query_point_indices = snap.bonds.group[:, 0]
    point_indices = snap.bonds.group[:, 1]
    # print(np.shape(query_point_indices))
    # print(np.shape(point_indices))
    # This may fix the error of the indices not being sorted, but uncertain
    # sorted_indices = np.argsort(query_point_indices)
    # query_point_indices = query_point_indices[sorted_indices]
    # point_indices = point_indices[sorted_indices]
    vectors = system.box.wrap(
        system.points[query_point_indices] - system.points[point_indices]
    )
    nlist = freud.NeighborList.from_arrays(
        num_query_points, num_points, query_point_indices, point_indices, vectors
    )
    cluster = freud.cluster.Cluster()
    cluster.compute(system=system, neighbors=nlist)
    return cluster.cluster_idx

def intermolecular_rdf(
    gsd_traj,
    A_name=None,
    B_name=None,
    r_max=None,
    r_min=0,
    bins=100,
    exclude_bonded=True,
):
    """Compute intermolecular RDF from a GSD file.

    This function calculates the radial distribution function given a GSD frame
    and the names of the particle types. By default it will calculate the RDF
    for all particles.

    It is assumed that the bonding, number of particles, and simulation box do
    not change during the simulation.

    Parameters
    ----------
    gsd_traj : gsd.hoomd.Trajectory
        The GSD trajectory.
    A_name, B_name : str
        Name(s) of particles between which to calculate the RDF (found in
        gsd.hoomd.Snapshot.particles.types)
    r_max : float
        Maximum radius of RDF. If None, half of the maximum box size is used.
        (Default value = None)
    r_min : float
        Minimum radius of RDF. (Default value = 0)
    bins : int
        Number of bins to use when calculating the RDF. (Default value = 100)
    exclude_bonded : bool
        Whether to remove particles in same molecule from the neighbor list.
        (Default value = True)

    Returns
    -------
    freud.density.RDF
    """
    import freud
    
    snap = gsd_traj[0]

    if r_max is None:
        # Use a value just less than half the maximum box length.
        r_max = np.nextafter(
            np.max(snap.configuration.box[:3]) * 0.5, 0, dtype=np.float32
        )

    rdf = freud.density.RDF(bins=bins, r_max=r_max, r_min=r_min)

    if A_name == None and B_name == None:
        type_A = [True]*snap.particles.N
        type_B = [True]*snap.particles.N
    else:
        type_A = snap.particles.typeid == snap.particles.types.index(A_name)
        type_B = snap.particles.typeid == snap.particles.types.index(B_name)

    if exclude_bonded:
        molecules = snap_molecule_indices(snap)
        molecules_A = molecules[type_A]
        molecules_B = molecules[type_B]

    for snap in gsd_traj:

        # operate on the frame
        A_pos = snap.particles.position[type_A]
        if A_name == B_name:
            B_pos = A_pos
            exclude_ii = True
        else:
            B_pos = snap.particles.position[type_B]
            exclude_ii = False

        box = snap.configuration.box
        system = (box, A_pos)
        aq = freud.locality.AABBQuery.from_system(system)
        nlist = aq.query(
            B_pos, {"r_max": r_max, "exclude_ii": exclude_ii}
        ).toNeighborList()

        if exclude_bonded:
            pre_filter = len(nlist)
            indices_A = molecules_A[nlist.point_indices]
            indices_B = molecules_B[nlist.query_point_indices]
            nlist.filter(indices_A != indices_B)
            post_filter = len(nlist)

        # rdf.compute(aq, neighbors=nlist, reset=False)
        rdf.compute(aq, reset=False)
    normalization = post_filter / pre_filter if exclude_bonded else 1
    return rdf, normalization

def reduce_network_dangling_ends(in_graph):
    '''
    Does one iteration of removing dangling ends from a network.
    Dangling ends are defined as nodes with only one edge, and must be evaluated recursively
        to account for branching dangling ends. (see remove_dangling_ends())
    Args:
        in_graph: networkx.Graph instance
    Returns:
        out_graph: networkx.Graph instance with one iteration of dangling ends removed
        dangling_ends: list of beads that were a part of dangling ends removed from the network
    '''
    dangling_ends = []
    out_graph = in_graph.copy()
    for i in in_graph.nodes():
        if in_graph.degree(i) == 1:
            out_graph.remove_node(i)
            dangling_ends.append(i)
    return out_graph, dangling_ends

def id_and_remove_dangling_and_primary(in_graph):
    """
    Args:
        in_graph: unmodified networkx.Graph instance
    Returns:
        new_graph: networkx.MuliGraph instance with dangling ends and primary loops removed
        dangling_ends: indexes and lengths of dangling ends removed from the network
        primary_loops: indexes and lengths of primary loops removed from the network
        super_dangling_ends: indexes and lengths of super dangling ends removed from the network
        super_primary_loops: indexes and lengths of super primary loops removed from the network

    Recursively:
        remove extenders
        remove dangling ends
        remove primary loops
    """
    #These will contain lists of the reduced network nodes with dangling ends
    dangling_end_nodes = []
    primary_loop_nodes = []

    dangling_ends = []
    primary_loops = []
    super_dangling_ends = []
    super_primary_loops = []

    new_graph = nx.MultiGraph()
    last_graph = remove_extenders(in_graph)
    scanlan_case_graph = last_graph.copy()

    for i in range(10000):
        #remove extenders
        new_graph = remove_extenders(last_graph)

        #remove dangling ends
        # new_graph,dangling_ends = reduce_network_dangling_ends(last_graph)
        nodes_to_remove = []
        edges_to_remove = []
        for node in new_graph.nodes():
            #find nodes with only one neighbor
            if len(list(new_graph.neighbors(node))) == 1:
                # edge = new_graph.edges(node)
                # tmp = 0
                # for k in range(len(edge)):
                #     edge_weight = edge[k]['weight']
                #     tmp += edge_weight
                # check the lists for any already associated defects
                # add the weight of the edge to the correct list (dangling end or super dangling end)
                # dangling_ends.append([in_graph.neighbors(node)[0],tmp])
                #find any other nodes that were previously truncated
                if len(list(scanlan_case_graph.neighbors(node))) == 1:
                    dangling_end_nodes.append([list(new_graph.neighbors(node))[0],node])
                else:
                    for i in list(scanlan_case_graph.neighbors(node)):
                        # check if the neighboring nodes are in the dangling end list or primary loop list
                        


                        print("REPLACE ME")



                    dangling_end_nodes.append([list(new_graph.neighbors(node))[0],node,])
                nodes_to_remove.append(node)

            # if there is a primary loop edge,
            if new_graph.has_edge(u=node,v=node):
                primary_loop_nodes.append(node)
                edges_to_remove.append((node,node))
            

        new_graph.remove_nodes_from(nodes_to_remove)

                



    
        if new_graph.edges() == last_graph.edges():
            break
        else:
            last_graph = new_graph.copy()

def id_ineffective_junctions(in_graph):
    """
    Args:
        in_graph: networkx.MultiGraph of only junctions, terminal beads, and weighted edges
    Returns:
        effective_graph: networkx.MultiGraph instance with only effective junctions and edges
        ineffective_graph: networkx.MultiGraph instance with only ineffective junctions and edges

    Recursively:
        remove extenders
        remove dangling ends
        remove primary loops
    """
    #Assume all nodes are effective
    effective_graph = in_graph.copy()
    ineffective_graph = nx.MultiGraph()
    ineffective_nodes = []
    ineffective_edges = []

    last_graph = in_graph.copy()
    new_graph = nx.MultiGraph()
    # scanlan_case_graph = last_graph.copy()

    for i in range(100000):
        # draw_graph(last_graph)
        # plt.show()
        new_graph = last_graph.copy()

        nodes_to_remove = []
        edges_to_remove = []
        for node in new_graph.nodes():
            neighbors = list(new_graph.neighbors(node))
            # if the node has only one neighbor and is not a primary loop (degree 2),
            if len(neighbors) == 1 and new_graph.degree(node) == 1:
                # effective_graph.remove_node(node)
                # ineffective_nodes.append(node)
                # ineffective_edges.append((new_graph.neighbors(node)[0],node))
                nodes_to_remove.append(node)
                edges_to_remove.append((neighbors[0],node))
                ineffective_graph.add_edge(node,neighbors[0],weight=new_graph[node][neighbors[0]][0]["weight"])
                ineffective_nodes.append(node)
                print("dangling_end: ",node,neighbors[0])
                break

            # if the node has a primary loop edge,
            if new_graph.has_edge(u=node,v=node):
                weight = new_graph[node][node][0]['weight']
                # print(new_graph[node][node][0]['weight'])
                ineffective_graph.add_edge(node,node,weight=weight)
                # ineffective_edges.append((node,node))
                edges_to_remove.append((node,node))
                ineffective_nodes.append(node) # NOTE: this node can only be categorized as ineffective because the maximum functionality is 4
                print("primary_loop: ",node)
                break

            # if the node is only an extender, then it should be in the ineffective graph, not effective
            if new_graph.degree(node) == 2:
                # to avoid fringe cases where the node is part of a primary loop, so it has degree 2, but only one neighbor
                if len(neighbors) == 2:
                    # connect the two neighbors
                    weight = new_graph[node][neighbors[0]][0]['weight'] + \
                                new_graph[node][neighbors[1]][0]['weight'] + 1
                    new_graph.add_edge(neighbors[0],neighbors[1],weight=weight)
                    ineffective_graph.add_edge(node,neighbors[0],weight=new_graph[node][neighbors[0]][0]['weight'])
                    ineffective_graph.add_edge(node,neighbors[1],weight=new_graph[node][neighbors[1]][0]['weight'])
                elif len(neighbors) == 1:
                    # connect the neighbor to itself
                    weight = new_graph[node][neighbors[0]][0]['weight'] + \
                                new_graph[node][neighbors[0]][1]['weight'] + 1
                    new_graph.add_edge(neighbors[0],neighbors[0],weight=weight)
                    ineffective_graph.add_edge(node,neighbors[0],weight=new_graph[node][neighbors[0]][0]['weight'])
                    ineffective_graph.add_edge(node,neighbors[0],weight=new_graph[node][neighbors[0]][1]['weight'])
                else:
                    raise ValueError("Node has more than 2 neighbors, but is not a primary loop")
                # remove the node
                ineffective_graph.add_node(node)
                nodes_to_remove.append(node)
                print("extender: ",node)
                ineffective_nodes.append(node)
                break

        new_graph.remove_nodes_from(nodes_to_remove)
        new_graph.remove_edges_from(edges_to_remove)
    
        # if new_graph.edges() == last_graph.edges():
        #     break
        # else:
        #     last_graph = new_graph.copy()
        if nodes_to_remove == [] and edges_to_remove == []:
            print(new_graph.edges())
            print(last_graph.edges())
            break
        else:
            last_graph = new_graph.copy()
    else: #nobreak
        raise ValueError("Dangling ends not removed after 10,000 iterations. Consider increasing the limit, or check for infinite recursion.")
    effective_graph = new_graph.copy()
    # plt.show()
    ineffective_graph = nx.MultiGraph()
    ineffective_graph.add_nodes_from(ineffective_nodes)
    ineffective_graph.add_edges_from(ineffective_edges)
    for edge in in_graph.edges(data=True):
        print(edge)
        if edge[0] or edge[1] in ineffective_nodes:
            ineffective_graph.add_edge(edge[0],edge[1],weight=edge[-1]['weight'])
    return effective_graph, ineffective_graph

#returns a MULTIGRAPH
def remove_extenders(in_graph):
    ''''
    Recursively removes all nodes with degree 2, replacing them with a single edge 
        connecting the two neighbors. Essentially generates a CROSSLINK-ONLY graph.
    Args:
        in_graph: networkx.Graph instance
    Returns:
        out_graph: networkx.MultiGraph instance with only crosslinks
    '''
    # out_graph = in_graph.copy().to_multigraph()
    out_graph = nx.MultiGraph()
    out_graph.add_edges_from(in_graph.edges(),weight=0)
    out_graph.add_nodes_from(in_graph.nodes(),weight=0)
    strands,crosslinkers = cut_into_strands(in_graph)
    for strand in strands:
        #get the nodes that are connected to crosslinkers
        #get the crosslinkers that the end nodes are connected to
        end_nodes = []
        strand_crosslinkers = []
        for node in strand:
            neighbors = list(in_graph.neighbors(node))
            for neighbor in neighbors:
                if neighbor not in strand:
                    strand_crosslinkers.append(neighbor)
            if len(neighbors) == 1:
                end_nodes.append(node)
        # remove all the nodes in the strand, replacing them with an edge connecting the crosslinkers
        out_graph.remove_nodes_from(strand)
        if len(strand_crosslinkers) == 2:
            # print("adding edge between:",strand_crosslinkers[0],strand_crosslinkers[1])
            out_graph.add_edge(strand_crosslinkers[0],strand_crosslinkers[1],weight=len(strand))
        elif len(strand_crosslinkers) == 1:
            # print("adding edge between:",strand_crosslinkers[0],end_nodes[0])
            out_graph.add_edge(strand_crosslinkers[0],end_nodes[0],weight=len(strand))
        else:
            print(f"Strand {strand} has {len(strand_crosslinkers)} ends/crosslinkers: {strand_crosslinkers}")
            raise ValueError("Strand has more than 2 ends/crosslinkers")
    return out_graph

def cut_into_strands(in_graph):
    '''
    Takes a graph and cuts it into strands which end right before any crosslinkers with 
        functionality 3 or higher

    Args:
        in_graph: networkx.Graph instance
    Returns:
        Gcc_new: list of networkx.Graph instances, each representing a strand sorted largest to smallest
        crosslink_beads: list of beads that were removed from the network
    '''
    Gnew = in_graph.copy()
    crosslink_beads = [x for  x in in_graph.nodes() if in_graph.degree(x) >= 3]
    
    for x in crosslink_beads:
        Gnew.remove_node(x)

    # strands
    Gcc_new = sorted(nx.connected_components(Gnew), key=len, reverse=True)
    return Gcc_new,crosslink_beads

def remove_dangling_ends(in_graph):
    '''
    Recursively removes dangling ends from a network. 
    Dangling ends are defined as nodes with only one edge, and must be evaluated recursively
        to account for branching dangling ends.
    Args:
        in_graph: unmodified networkx.Graph instance
    Returns:
        new_graph: networkx.Graph instance with dangling ends removed
        dangling_ends: list of beads that were a part of dangling ends removed from the network
    '''
    dangling_ends = []
    new_graph = nx.Graph()
    last_graph = in_graph.copy()
    for i in range(10000):
        new_graph, new_dangling_ends = reduce_network_dangling_ends(last_graph)
        dangling_ends.extend(new_dangling_ends)
        if new_graph.edges() == last_graph.edges():
            break
        else:
            last_graph = new_graph.copy()
    else:
        raise ValueError("Dangling ends not removed after 10,000 iterations. Consider increasing the limit, or check for infinite recursion.")
    return new_graph, dangling_ends

# To contract chains of neighbouring vertices with degree 2 into one hypernode:
#https://stackoverflow.com/questions/52313551/graph-reduction

def is_in_2D(list, element):
    """
    Checks if an element exists in a two-dimensional list.

    Args:
        list: The two-dimensional list to search within.
        element: The element to search for.

    Returns:
        The index of the sublist containing the element, or -1 if not found.
    """
    sublist_index = 0
    for sublist in list:
        if element in sublist:
            return sublist_index
        sublist_index += 1
    return -1 #Not found

def categorize_loop_edges(loops,crosslink_graph):
    '''
    Categorizes the edges in each loop into primary, secondary, tertiary, and quaternary+ loops.
    Assumes that the crosslink_graph is a graph of only crosslinks, where the nodes in 
        between have been replaced with a single edge.
    Primary loops are loops that contain only one edge.
    Secondary loops are loops that contain two edges.
    Tertiary loops are loops that contain three edges.
    Quaternary+ loops are loops that contain four or more edges.

    If an edge is shared between two types of loop, it is categorized as the lower degree type.
    
    Args:
        loops: list of lists of nodes in each loop
        crosslink_graph: networkx.Graph instance of ONLY CROSSLINKS (see function remove_extenders)
    Returns:
        primary_edges: list of lists of edges in primary loops
        secondary_edges: list of lists of edges in secondary loops
        tertiary_edges: list of lists of edges in tertiary loops
        quaternary_edges: list of lists of edges in quaternary+ loops
    '''
    #graph needs to be the graph of only crosslinks
    loop_edges = np.empty((len(loops),),dtype=object)
    loop_types = np.array([classify_loop_type_simple(l) for l in loops])
    #returns lists of the edges in each loop type
    primary_edges = []
    secondary_edges = []
    tertiary_edges = []
    quaternary_edges = []

    #Compile lists of the edges in each loop
    for i,loop in enumerate(loops):
        loop_edges[i] = []
        if len(loop) == 1:
            loop_edges[i].append(crosslink_graph.get_edge_data(u=loop[0],v=loop[0]))
        else:
            for j in range(len(loop)):
                for k in range(j+1,len(loop)):
                    bead0=loop[j]
                    bead1=loop[k]
                    loop_edges[i].append(crosslink_graph.get_edge_data(u=bead0,v=bead1))
        if loop_types[i] == "primary":
            primary_edges.append(loop_edges[i])
        elif loop_types[i] == "secondary":
            secondary_edges.append(loop_edges[i])
        elif loop_types[i] == "tertiary":
            tertiary_edges.append(loop_edges[i])
        elif loop_types[i] == "quaternary+":
            quaternary_edges.append(loop_edges[i])

    #Remove double counting of edges
    for i in range(len(primary_edges)):
        test_index = is_in_2D(element=primary_edges[i], list=secondary_edges)
        if (test_index != -1):
            secondary_edges[test_index].remove(primary_edges[i])
        test_index = is_in_2D(element=primary_edges[i], list=tertiary_edges)
        if (test_index != -1):
            tertiary_edges[test_index].remove(primary_edges[i])
        test_index = is_in_2D(element=primary_edges[i], list=quaternary_edges)
        if (test_index != -1):
            quaternary_edges[test_index].remove(primary_edges[i])
    for i in range(len(secondary_edges)):
        for j in range(len(secondary_edges[i])):
            test_index = is_in_2D(element=secondary_edges[i][j], list=tertiary_edges)
            if (test_index != -1):
                tertiary_edges[test_index].remove(secondary_edges[i][j])
            test_index = is_in_2D(element=secondary_edges[i][j], list=quaternary_edges)
            if (test_index != -1):
                quaternary_edges[test_index].remove(secondary_edges[i][j])
        # if (is_in_2D(element=secondary_edges[i], list=tertiary_edges):
        #     tertiary_edges.remove(secondary_edges[i])
        # if (is_in_2D(element=secondary_edges[i], list=quaternary_edges):
        #     quaternary_edges.remove(secondary_edges[i])
    for i in range(len(tertiary_edges)):
        test_index = is_in_2D(element=tertiary_edges[i], list=quaternary_edges)
        if test_index != -1:
            quaternary_edges[test_index].remove(tertiary_edges[i])

    for i in range(len(quaternary_edges)):
        while None in quaternary_edges[i]:
            quaternary_edges[i].remove(None)
        # if quaternary_edges[i] == []:
        #     quaternary_edges.pop(i)
        # if is_in_2D(element=None,list=quaternary_edges) != -1:
        #     print(quaternary_edges[index])
        #     quaternary_edges.pop(index)

    return primary_edges, secondary_edges, tertiary_edges, quaternary_edges

def get_length_categorized_loops(primary_edges,secondary_edges,tertiary_edges,quaternary_edges):
    '''
    Returns the length of each loop in each category of loop type.
    Args:
        primary_edges: list of lists of edges in primary loops
        secondary_edges: list of lists of edges in secondary loops
        tertiary_edges: list of lists of edges in tertiary loops
        quaternary_edges: list of lists of edges in quaternary+ loops
    Returns:
        primary_lengths: list of lengths of primary loops
        secondary_lengths: list of lengths of secondary loops
        tertiary_lengths: list of lengths of tertiary loops
        quaternary_lengths: list of lengths of quaternary+ loops
    '''
    #graph needs to be the graph of only crosslinks
    primary_lengths = np.zeros(len(primary_edges))
    secondary_lengths = np.zeros(len(secondary_edges))
    tertiary_lengths = np.zeros(len(tertiary_edges))
    quaternary_lengths = np.zeros(len(quaternary_edges))

    for i,loop_edges in enumerate(primary_edges):
        # print("primary edges:",edges)
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        primary_lengths[i] = tmp

    for i,loop_edges in enumerate(secondary_edges):
        # print("secondary edges:",edges)
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        secondary_lengths[i] = tmp

    for i,loop_edges in enumerate(tertiary_edges):
        # print("tertiary edges:",edges)
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        tertiary_lengths[i] = tmp

    for i,loop_edges in enumerate(quaternary_edges):
        tmp = 0
        for edge in loop_edges:
            for k in range(len(edge)):
                edge_weight = edge[k]['weight']
                tmp += edge_weight
        quaternary_lengths[i] = tmp

    return primary_lengths, secondary_lengths, tertiary_lengths, quaternary_lengths

#---------------------------------------------------------------------------------------
# Analysis Functions
#---------------------------------------------------------------------------------------

#Dangling end lengths do not include the beads within the gel
#
def defect_analysis(job_id,testing=False, defect_frame = -1):
    
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

    # if not analyzing the final frame
    if defect_frame != -1:
        from scripts.conversion import get_frame_at_95_conversion
        defect_frame = int(get_frame_at_95_conversion(job_id))
        print(defect_frame)


    frame = trajectory[defect_frame]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]


    Ntotal = frame.particles.N

    if testing:
        jdir = "./test_network_properties"
        bonds = [[0,1],
                 [1,1.5],[1.5,2],[1,1.6],[1.6,2],[11,2.11],[2,2.11],
                 [2,3],[3,4],[4,1],[0,5],[5,6],[6,7],[5,8],[8,9],
                    [2,10],[2,11],[10,12],[12,11],[12,13],[12,14],[13,15],[15,14],[0,16],[0,11],[6,17],[1,1]]
        print("bonds:",bonds)
        G = nx.Graph()
        G.add_edges_from(bonds)
        G2 = nx.Graph()
        G2.add_edges_from(bonds)
        frame = gsd.hoomd.Frame()
        frame.particles.position = np.array([[0,0,0],[1,0,0],[1.5,0,0],[1.6,0,0],[2,0,0],[2.11,0,0],[2,1,0],[1,1,0],[-1,0,0],[-2,0,0],[-3,0,0],[-1,1,0],[-1,2,0],
                                             [3,0,0],[2,1,0],[3,1,0],[4,1,0],[3,2,0],[4,2,0],[-1,0,0],[-2,-1,0]])
        frame.particles.velocity = np.zeros((len(frame.particles.position),3))
        frame.particles.N = len(frame.particles.position)
        frame.bonds.N = len(bonds)
        frame.bonds.group = np.array(bonds)
        frame.bonds.typeid = np.zeros(len(bonds))
        frame.bonds.types = ['bond']
        new_bonds = bonds

        # print('G2 Nodes:',G2.nodes())
        # print('G2 Edges:',G2.edges())
        # G2 = remove_extenders(G2)
        # print('G2 Nodes:',G2.nodes())
        # print('G2 Edges:',G2.edges())

        '''
        trajectory = gsd.hoomd.open('/home/bj21/simulations/thiol-ene/network_polymerization/workspace/b00fe99f07ce7c0bc08a9cf517888351/polymerize.gsd')

        frame = trajectory[-1]
        dummy_id = len(frame.bonds.types)-1
        bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
        '''

        G = nx.Graph()
        G.add_edges_from(bonds)
        # find all disconnected/connected sub-networks
        # sort the list of networks by length
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)

        # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
        # bonds in biggest cluster:
        new_bonds = []
        for each in bonds:
            for i in each:
                if i in Gcc[0]:
                    new_bonds.append(each)
                    continue
        # make cluster of biggest graph
        G2 = nx.Graph()
        G2.add_edges_from(new_bonds)
    else:
        # make a graph - each bond is an edge
        G = nx.Graph()
        G.add_edges_from(bonds)
        # find all disconnected/connected sub-networks
        # sort the list of networks by length
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
            
        # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
        # bonds in biggest cluster:
        new_bonds = []
        for each in bonds:
            for i in each:
                if i in Gcc[0]:
                    new_bonds.append(each)
                    continue
        # make cluster of biggest graph
        G2 = nx.Graph()
        G2.add_edges_from(new_bonds)
    
    print("calculating dangling ends")
    dangling_end_strands = []
    dangling_end_lengths = []

    no_danglingEnds_graph, dangling_ends = remove_dangling_ends(G2)
    # print('no_danglingEnds_graph edges', no_danglingEnds_graph.edges())
    # print('no_danglingEnds_graph nodes', no_danglingEnds_graph.nodes())
    dangling_end_graph = G2.copy()
    for x in G2.nodes():
        if x not in dangling_ends:
            dangling_end_graph.remove_node(x)
    # print("dangling ends:",dangling_end_graph.nodes())
    # print("dangling end edges:",dangling_end_graph.edges())

    for dangling_strand in nx.connected_components(dangling_end_graph):
        dangling_strand = list(dangling_strand)
        # print(dangling_strand)
        dangling_end_strands.append(dangling_strand)
        dangling_end_lengths.append(len(dangling_strand))
    unique_dangling_end_lengths,count = np.unique(dangling_end_lengths,return_counts=True)
    dangling_end_data = list(zip(unique_dangling_end_lengths,count))
    # print(dangling_end_data)

    print("calculating loops",  flush=True)
    start = datetime.now()
    crosslink_graph = remove_extenders(no_danglingEnds_graph)
    # print('crosslink_graph edges', crosslink_graph.edges())
    # print('crosslink_graph nodes', crosslink_graph.nodes())
    print('start loop analysis:',start,  flush=True)
    loops = list(nx.simple_cycles(crosslink_graph,length_bound=4))
    print('end loop analysis:',datetime.now(),  flush=True)
    print('duration:',datetime.now()-start,  flush=True)
    # lengths_loops = np.array([len(l) for l in loops])
    loop_types = np.array([classify_loop_type_simple(l) for l in loops])
    primary_edges, secondary_edges, tertiary_edges, quaternary_edges = categorize_loop_edges(loops, crosslink_graph)
    primary_lengths, secondary_lengths, tertiary_lengths, quaternary_lengths = get_length_categorized_loops(primary_edges,secondary_edges,tertiary_edges,quaternary_edges)
    # print('primary edges:',primary_edges)
    # print('secondary edges:',secondary_edges)
    # print('tertiary edges:',tertiary_edges)
    # print('quaternary edges:',quaternary_edges)
    #if testing, print out the ids of the particles participating in loops

    if (testing):
        for i in range(len(loops)):
            print("identified loop:",loops[i]," as type:",loop_types[i])

    if (testing):
        out_frame = frame
        out_frame.particles.velocity[:] = -1
        for strand in dangling_end_strands:
            for bead in strand:
                out_frame.particles.velocity[bead] = [0,0,0]
        for i,loop in enumerate(loops):
            match loop_types[i]:
                case "primary":
                    type = 1
                case "secondary":
                    type = 2
                case "tertiary":
                    type = 3
                case "quaternary+":
                    type = 4
            for bead in loop:
                out_frame.particles.velocity[bead] = [type,0,0]
        output_gsd = gsd.hoomd.open(jdir + f'/test_{os.path.basename(__file__)}'.replace('.py','.gsd'), 'w')
        output_gsd.append(out_frame)

    #Format the loop data as discrete histograms
    unique_primary_lengths,count = np.unique(primary_lengths,return_counts=True)
    primary_data = list(zip(unique_primary_lengths,count))
    print('primary lengths and counts:',primary_data)
    unique_secondary_lengths,count = np.unique(secondary_lengths,return_counts=True)
    secondary_data = list(zip(unique_secondary_lengths,count))
    print('secondary lengths and counts:',secondary_data)
    unique_tertiary_lengths,count = np.unique(tertiary_lengths,return_counts=True)
    tertiary_data = list(zip(unique_tertiary_lengths,count))
    print('tertiary lengths and counts:',tertiary_data)
    unique_quaternary_lengths,count = np.unique(quaternary_lengths,return_counts=True)
    quaternary_data = list(zip(unique_quaternary_lengths,count))
    print('quarternary lengths and counts:', quaternary_data)

    #write the data to a file
    # Here, each line is a tuple of the form (loop_size, count), with header lines 
    # indicating the type of loop data to follow until the next header

    # if not analyzing the final frame
    if defect_frame != -1:
        file_name = jdir + "/loop_counts_conversion95.txt"
    else:
        file_name = jdir + "/loop_counts.txt"

    with open(file_name, "w") as f:
        f.write("dangling_strand data\n")
        for x in dangling_end_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("primary_loop_size data\n")
        for x in primary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("secondary_loop_size data\n")
        for x in secondary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("tertiary_loop_size data\n")
        for x in tertiary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("quaternary_loop_size data\n")
        for x in quaternary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

def defect_analysis_in_progress(job_id,testing=False):
    
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


    Ntotal = frame.particles.N

    if testing:
        jdir = "./test_network_properties"
        bonds = [[0,1],
                 [1,1.5],[1.5,2],[11,2.11],[2,2.11],
                 [2,3],[3,4],[4,1],[0,5],[5,6],[6,7],[5,8],[8,9],
                    [2,10],[10,12],[12,11],[12,13],[12,14],[13,15],[15,14],[0,16],[0,11],[6,17],[8,8]]
        print("bonds:",bonds)
        G = nx.Graph()
        G.add_edges_from(bonds)
        G2 = nx.Graph()
        G2.add_edges_from(bonds)
        frame = gsd.hoomd.Frame()
        frame.particles.position = np.array([[0,0,0],[1,0,0],[1.5,0,0],[1.6,0,0],[2,0,0],[2.11,0,0],[2,1,0],[1,1,0],[-1,0,0],[-2,0,0],[-3,0,0],[-1,1,0],[-1,2,0],
                                             [3,0,0],[2,1,0],[3,1,0],[4,1,0],[3,2,0],[4,2,0],[-1,0,0],[-2,-1,0]])
        frame.particles.velocity = np.zeros((len(frame.particles.position),3))
        frame.particles.N = len(frame.particles.position)
        frame.bonds.N = len(bonds)
        frame.bonds.group = np.array(bonds)
        frame.bonds.typeid = np.zeros(len(bonds))
        frame.bonds.types = ['bond']
        new_bonds = bonds

        # print('G2 Nodes:',G2.nodes())
        # print('G2 Edges:',G2.edges())
        # G2 = remove_extenders(G2)
        # print('G2 Nodes:',G2.nodes())
        # print('G2 Edges:',G2.edges())

        '''
        trajectory = gsd.hoomd.open('/home/bj21/simulations/thiol-ene/network_polymerization/workspace/b00fe99f07ce7c0bc08a9cf517888351/polymerize.gsd')

        frame = trajectory[-1]
        dummy_id = len(frame.bonds.types)-1
        bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
        '''

        G = nx.Graph()
        G.add_edges_from(bonds)
        # find all disconnected/connected sub-networks
        # sort the list of networks by length
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)

        # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
        # bonds in biggest cluster:
        new_bonds = []
        for each in bonds:
            for i in each:
                if i in Gcc[0]:
                    new_bonds.append(each)
                    continue
        # make cluster of biggest graph
        gel_graph = nx.Graph()
        gel_graph.add_edges_from(new_bonds)
    else:
        # make a graph - each bond is an edge
        G = nx.Graph()
        G.add_edges_from(bonds)
        # find all disconnected/connected sub-networks
        # sort the list of networks by length
        Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
            
        # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
        # bonds in biggest cluster:
        new_bonds = []
        for each in bonds:
            for i in each:
                if i in Gcc[0]:
                    new_bonds.append(each)
                    continue
        # make cluster of biggest graph
        gel_graph = nx.Graph()
        gel_graph.add_edges_from(new_bonds)
    

    draw_graph(gel_graph)
    plt.show()
    gel_graph = remove_extenders(gel_graph)
    draw_graph(gel_graph)
    plt.show()

    effective_graph, ineffective_graph = id_ineffective_junctions(gel_graph)
    draw_graph(effective_graph)
    plt.show()
    draw_graph(ineffective_graph)
    plt.show()
    exit()









    print("calculating dangling ends")
    dangling_end_strands = []
    dangling_end_lengths = []

    no_danglingEnds_graph, dangling_ends = remove_dangling_ends(G2)
    # print('no_danglingEnds_graph edges', no_danglingEnds_graph.edges())
    # print('no_danglingEnds_graph nodes', no_danglingEnds_graph.nodes())
    dangling_end_graph = G2.copy()
    for x in G2.nodes():
        if x not in dangling_ends:
            dangling_end_graph.remove_node(x)
    # print("dangling ends:",dangling_end_graph.nodes())
    # print("dangling end edges:",dangling_end_graph.edges())

    for dangling_strand in nx.connected_components(dangling_end_graph):
        dangling_strand = list(dangling_strand)
        # print(dangling_strand)
        dangling_end_strands.append(dangling_strand)
        dangling_end_lengths.append(len(dangling_strand))
    unique_dangling_end_lengths,count = np.unique(dangling_end_lengths,return_counts=True)
    dangling_end_data = list(zip(unique_dangling_end_lengths,count))
    # print(dangling_end_data)

    print("calculating loops",  flush=True)
    start = datetime.now()
    crosslink_graph = remove_extenders(no_danglingEnds_graph)
    # print('crosslink_graph edges', crosslink_graph.edges())
    # print('crosslink_graph nodes', crosslink_graph.nodes())
    print('start loop analysis:',start,  flush=True)
    loops = list(nx.simple_cycles(crosslink_graph,length_bound=4),  flush=True)
    print('end loop analysis:',datetime.now(),  flush=True)
    print('duration:',datetime.now()-start,  flush=True)
    # lengths_loops = np.array([len(l) for l in loops])
    loop_types = np.array([classify_loop_type_simple(l) for l in loops])
    primary_edges, secondary_edges, tertiary_edges, quaternary_edges = categorize_loop_edges(loops, crosslink_graph)
    primary_lengths, secondary_lengths, tertiary_lengths, quaternary_lengths = get_length_categorized_loops(primary_edges,secondary_edges,tertiary_edges,quaternary_edges)
    # print('primary edges:',primary_edges)
    # print('secondary edges:',secondary_edges)
    # print('tertiary edges:',tertiary_edges)
    # print('quaternary edges:',quaternary_edges)
    #if testing, print out the ids of the particles participating in loops

    if (testing):
        for i in range(len(loops)):
            print("identified loop:",loops[i]," as type:",loop_types[i])

    if (testing):
        out_frame = frame
        out_frame.particles.velocity[:] = -1
        for strand in dangling_end_strands:
            for bead in strand:
                out_frame.particles.velocity[bead] = [0,0,0]
        for i,loop in enumerate(loops):
            match loop_types[i]:
                case "primary":
                    type = 1
                case "secondary":
                    type = 2
                case "tertiary":
                    type = 3
                case "quaternary+":
                    type = 4
            for bead in loop:
                out_frame.particles.velocity[bead] = [type,0,0]
        output_gsd = gsd.hoomd.open(jdir + f'/test_{os.path.basename(__file__)}'.replace('.py','.gsd'), 'w')
        output_gsd.append(out_frame)

    #Format the loop data as discrete histograms
    unique_primary_lengths,count = np.unique(primary_lengths,return_counts=True)
    primary_data = list(zip(unique_primary_lengths,count))
    print('primary lengths and counts:',primary_data)
    unique_secondary_lengths,count = np.unique(secondary_lengths,return_counts=True)
    secondary_data = list(zip(unique_secondary_lengths,count))
    print('secondary lengths and counts:',secondary_data)
    unique_tertiary_lengths,count = np.unique(tertiary_lengths,return_counts=True)
    tertiary_data = list(zip(unique_tertiary_lengths,count))
    print('tertiary lengths and counts:',tertiary_data)
    unique_quaternary_lengths,count = np.unique(quaternary_lengths,return_counts=True)
    quaternary_data = list(zip(unique_quaternary_lengths,count))
    print('quarternary lengths and counts:', quaternary_data)

    #write the data to a file
    # Here, each line is a tuple of the form (loop_size, count), with header lines 
    # indicating the type of loop data to follow until the next header
    with open(jdir + "/loop_counts.txt", "w") as f:
        f.write("dangling_strand data\n")
        for x in dangling_end_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("primary_loop_size data\n")
        for x in primary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("secondary_loop_size data\n")
        for x in secondary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("tertiary_loop_size data\n")
        for x in tertiary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)

        f.write("quaternary_loop_size data\n")
        for x in quaternary_data:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)


def xlink_rdf_analysis(job_id):
    print(
    '''
    This function is intended to calculate the rdf of crosslinks in the system.
    However, running this on only 5 frames results in a lot of noise. Thus it needs 
    multiple frames to be run on. For the RDF analysis, there cannot be any new crosslinks
    being formed over the course of this trajectory. This analysis may not be reliable
    enough without additional simulation.
    ''')
    testing = True


    project = signac.get_project()
    direc = project.fn('') + 'workspace/'

    # for job_id in os.listdir(direc):
    jdir = direc + job_id
    input_file = direc + job_id + "/polymerize.gsd"
    if not os.path.isfile(input_file):
        raise FileNotFoundError("File not found")
        exit()
    # open the polymerized gsd file's trajectory (list of frames)
    # trajectory = gsd.hoomd.open(input_file)
    trajectory = gsd.hoomd.open(input_file, mode='r')
    print(np.shape(trajectory[-1].particles.position))
    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]

    # make a graph - each bond is an edge
    G = nx.MultiGraph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        
    # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
    # bonds in biggest cluster:
    new_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                new_bonds.append(each)
                continue
    # make cluster of biggest graph
    G2 = nx.Graph()
    G2.add_edges_from(new_bonds)
    #get only the crosslink_beads
    ene_crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) == 3]
    thiol_crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) == 4]
    print("# of ene crosslinks:",len(ene_crosslink_beads))
    # print('ene_crosslink_beads:',ene_crosslink_beads)
    print("# of thiol crosslinks:",len(thiol_crosslink_beads))
    # print('thiol_crosslink_beads:',thiol_crosslink_beads)
    # exit()
    crosslinks_only_traj = []
    # make the new frames
    for frame in trajectory[-5:-1]:
        crosslinks_only_frame = gsd.hoomd.Frame()
        crosslinks_only_frame.particles.position = frame.particles.position[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.typeid = frame.particles.typeid[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.orientation = frame.particles.orientation[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.mass = frame.particles.mass[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.charge = frame.particles.charge[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.diameter = frame.particles.diameter[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.body = frame.particles.body[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.moment_inertia = frame.particles.moment_inertia[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.velocity = frame.particles.velocity[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.angmom = frame.particles.angmom[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.image = frame.particles.image[ene_crosslink_beads + thiol_crosslink_beads]
        crosslinks_only_frame.particles.N = len(ene_crosslink_beads) + len(thiol_crosslink_beads)
        crosslinks_only_frame.particles.types = frame.particles.types
        crosslinks_only_frame.configuration = frame.configuration

        crosslinks_only_traj.append(crosslinks_only_frame)
    #get the rdfs
    # 'C' = tetrafunctional thiol crosslinker center, 'Carbon' = chain growthed ene
    job = project.open_job(id=job_id)
    if job.sp.get('chain_side_reaction_probability') != 0.0:
        rdf_ene, norm_ene = intermolecular_rdf(crosslinks_only_traj, 
                                            A_name='Carbon', B_name='Carbon', 
                                            r_max=10.0, r_min=0.25, bins=500, 
                                            exclude_bonded=False,)
        with open(jdir + "/ene_ene_rdf.txt", "w") as f:
            f.write("r g(r)\n")
            for r, g in zip(rdf_ene.bin_centers, rdf_ene.rdf*norm_ene):
                f.write(f"{r} {g}\n")
    if job.sp.get("crosslinker_percent") != 0.0:
        rdf_thiol, norm_thiol = intermolecular_rdf(crosslinks_only_traj, 
                                                A_name='C', B_name='C', 
                                                r_max=10.0, r_min=0.25, bins=500, 
                                                exclude_bonded=False,)
        with open(jdir + "/thiol_thiol_rdf.txt", "w") as f:
            f.write("r g(r)\n")
            for r, g in zip(rdf_thiol.bin_centers, rdf_thiol.rdf*norm_thiol):
                f.write(f"{r} {g}\n")
    if job.sp.get("crosslinker_percent") != 0.0 and job.sp.get('chain_side_reaction_probability') != 0.0:
        rdf_thiol_ene, norm_thiol_ene = intermolecular_rdf(crosslinks_only_traj, 
                                                        A_name='C', B_name='Carbon', 
                                                        r_max=10.0, r_min=0.25, bins=500, 
                                                        exclude_bonded=False,)
        with open(jdir + "/thiol_ene_rdf.txt", "w") as f:
            f.write("r g(r)\n")
            for r, g in zip(rdf_thiol_ene.bin_centers, rdf_thiol_ene.rdf*norm_thiol_ene):
                f.write(f"{r} {g}\n")

    # save the rdf data into txts


def strand_lengths_analysis(job_id, defect_frame = -1):
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
    if defect_frame != -1:
        from scripts.conversion import get_frame_at_95_conversion
        defect_frame = int(get_frame_at_95_conversion(job_id))
        print(defect_frame)
    frame = trajectory[defect_frame]
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
    print(sizes[0][0])
    gf = sizes[0][0]/total

    if defect_frame == -1:
        file = "/gel_fraction.txt"
    else:
        file = "/conv95_gel_fraction.txt"

    with open(jdir + file, "w") as f:
        f.write(str(gf))

    new_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                new_bonds.append(each)
                continue
    # make cluster of biggest graph
    G2 = nx.Graph()
    G2.add_edges_from(new_bonds)
    Gnew = G2.copy()
    # remove everyone that has 3 or more bonds on it, only leaving linear strands
    crosslink_beads = [x for  x in G2.nodes() if G2.degree(x) >= 3]
    for x in crosslink_beads:
        Gnew.remove_node(x)
    
    # strands
    Gcc_new = sorted(nx.connected_components(Gnew), key=len, reverse=True)
    print("number of strands network",len(Gcc_new))
    sizes = [len(n) for n in Gcc_new]
    print("length of strands in network")
    sizes_of_strands,count = np.unique(sizes,return_counts=True)
    print(sizes_of_strands,count)

    strand_count = list(zip(sizes_of_strands,count))
    print(strand_count)
    if defect_frame == -1:
        file = "/strand_sizes.txt"
    else:
        file = "/conv95_strand_sizes.txt"
    with open(jdir + file, "w") as f:
        f.write("strand_size count\n")
        for x in strand_count:
            line = str(x[0]) + " " + str(x[1]) + "\n"
            f.write(line)
    
    '''
    crosslink_density = len(crosslink_beads)/sum([i[0]*i[1] for i in strand_count])
    with open(jdir + "/crosslink_density_1.txt", "w") as f:
        f.write(str(crosslink_density))'''

def calculate_crosslinking_density(job_id, analysis_frame = -1):
    project = signac.get_project()
    job = project.open_job(id=job_id)

    polymerize_gsd = job.fn('polymerize.gsd')
    trajectory = gsd.hoomd.open(polymerize_gsd)
    frame = trajectory[analysis_frame]
    box = frame.configuration.box
    volume = box[0]*box[1]*box[2]

    direc = project.fn('') + 'workspace/'
    if analysis_frame == -1: 
        strand_direc = direc + job_id + "/contract_bonds_analysis/"
    else:
        strand_direc = direc + job_id + "/percolation_contract_bonds_analysis/"
    effective_strands = np.loadtxt(job.fn(strand_direc + "effective_strand_hist.txt"), delimiter=",")
    total_effective = np.sum(effective_strands[:, 0])
    cld = total_effective/volume
    print(cld)
    if analysis_frame == -1:
        job.doc["crosslinking_density"] = cld
    else:
        job.doc["perc_crosslinking_density"] = cld





def contract_bonds_analysis(job_id, file="contract_bonds.gsd"):
    project = signac.get_project()
    direc = project.fn('') + 'workspace/'

    # for job_id in os.listdir(direc):
    jdir = direc + job_id
    input_file = direc + job_id + "/" + file

    if not os.path.isfile(input_file):
        raise FileNotFoundError("File not found")
        exit()
    # open the contract_bonds gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)

    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
    
    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        
    # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
    # bonds in biggest cluster:
    gel_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                gel_bonds.append(each)
                continue
    # make cluster of biggest graph
    gel_G = nx.Graph()
    gel_G.add_edges_from(gel_bonds)

    # get all bonds that are elastically ineffective
    ineffective_bonds = []
        # find the cutoff based on the maximum bond length of known ineffective bonds
    HS_SH_id = frame.bonds.types.index('-SH')
    C_C_double_id = frame.bonds.types.index('CC_Double')
    known_ineffective_bonds = frame.bonds.group[np.where((frame.bonds.typeid == HS_SH_id) | (frame.bonds.typeid == C_C_double_id))[0]]
    max_ineffective_length = 0
    for bond in known_ineffective_bonds:
        length = dist_pbc(frame.particles.position[bond[0]],\
                          frame.particles.position[bond[1]],\
                          Box=frame.configuration.box[0:3])
        if length > max_ineffective_length:
            max_ineffective_length = length
        # Use the maximum ineffective bond length as the cutoff
    for bond in gel_G.edges():
        if dist_pbc(frame.particles.position[bond[0]],\
                    frame.particles.position[bond[1]],\
                    Box=frame.configuration.box[0:3]) <= max_ineffective_length:
                    # and bond not in ineffective_bonds:
            # print("ineffective bond:",bond)
            ineffective_bonds.append(bond)

    maybe_ineffective_particles, ineffective_particle_bonds = np.unique(ineffective_bonds,return_counts=True)
    print("len(maybe_ineffective_particles):",len(maybe_ineffective_particles))
    print("ineffective_particle_bonds:",ineffective_particle_bonds)
    ineffective_particles = []
    for i in range(len(maybe_ineffective_particles)):
        # if the particle has the same number of ineffective bonds as the number of bonds on it, it is truly inefffective
        if ineffective_particle_bonds[i] == gel_G.degree(maybe_ineffective_particles[i]):
            # ineffective_particles.remove(ineffective_particles[i])
            ineffective_particles.append(maybe_ineffective_particles[i])
        # else:
        #     print("effective particle:",maybe_ineffective_particles[i])
    print("len(ineffective_particles):",len(ineffective_particles))
    print("elastically effective portion: ",1-len(ineffective_particles)/len(gel_G.nodes()))

    # make a copy graph of elastically effective network
    gel_G_effective = gel_G.copy()
    for particle in ineffective_particles:
        gel_G_effective.remove_node(particle)

    # make a copy graph of ineffective strands
    gel_ineffective_G = nx.Graph()
    gel_ineffective_G.add_edges_from(ineffective_bonds)
    nodes_to_remove = []
    for node in gel_ineffective_G.nodes():
        if node not in ineffective_particles:
            nodes_to_remove.append(node)
    for node in nodes_to_remove:
        gel_ineffective_G.remove_node(node)

    #identify the strands in effective graph
    # remove everyone that has 3 or more bonds on it, only leaving linear strands
    crosslink_beads = [x for  x in gel_G_effective.nodes() if gel_G_effective.degree(x) >= 3]
    G_tmp = gel_G_effective.copy()
    for x in crosslink_beads:
        G_tmp.remove_node(x)
    strands_effective = sorted(nx.connected_components(G_tmp), key=len, reverse=True)
    sizes_effective = [len(n) for n in strands_effective]
    effective_strand_hist, effective_strand_binEdges = np.histogram(sizes_effective,bins=np.arange(0,np.max(sizes_effective)+1))
    effective_strand_hist = effective_strand_hist.astype(float)
    print("effective strands:",effective_strand_hist)
    print("effective strands bin edges:",effective_strand_binEdges)

    #identify the strands in ineffective graph
    strands_ineffective = sorted(nx.connected_components(gel_ineffective_G), key=len, reverse=True)
    sizes_ineffective = [len(n) for n in strands_ineffective]
    ineffective_strand_hist, ineffective_strand_binEdges = np.histogram(sizes_ineffective,bins=np.arange(0,np.max(sizes_ineffective)+1))
    ineffective_strand_hist = ineffective_strand_hist.astype(float)
    print("ineffective strands:",ineffective_strand_hist)
    print("ineffective strands bin edges:",ineffective_strand_binEdges)

    # output a frame with the ineffective and effective strands colored
    testing = False
    if testing:
        out_frame = frame
        out_frame.particles.velocity[:] = -1
        for i in range(len(strands_ineffective)):
            for bead in strands_ineffective[i]:
                out_frame.particles.velocity[bead] = [i,0,0]
        for i in range(len(strands_effective)):
            for bead in strands_effective[i]:
                out_frame.particles.velocity[bead] = [0,i,0]
        output_gsd = gsd.hoomd.open(jdir + f'/test_contract_bonds_analysis.gsd', 'w')
        output_gsd.append(out_frame)

    # save the data
    if file == "contract_bonds.gsd":
        w_direc = jdir + "/contract_bonds_analysis/"
    else:
        w_direc = jdir + "/percolation_contract_bonds_analysis/"
    if not os.path.exists(w_direc):
        os.mkdir(w_direc)

    np.savetxt(w_direc + "ineffective_strand_hist.txt", \
               np.vstack((ineffective_strand_hist,np.arange(len(ineffective_strand_hist)))).T, \
                fmt='%d', delimiter=',', newline='\n')
    
    np.savetxt(w_direc + "effective_strand_hist.txt", \
                np.vstack((effective_strand_hist,np.arange(len(effective_strand_hist)))).T, \
                 fmt='%d', delimiter=',', newline='\n')
    
    with open(w_direc + "overall.txt", "w") as f:
        f.write("n_effective_particles\tn_ineffective_particles\n")
        f.write(f"{len(gel_G.nodes())-len(ineffective_particles)}\t{len(ineffective_particles)}\n")

def scanlan_case_analysis_on_contract_bonds(job_id, file="contract_bonds.gsd"):
    """
    Analyze the contract_bonds.gsd file to find the effective network from the 
    perspective of the crosslinks and scanlan case criterion.

    Outcome:
        Summary crosslink properties in the contract_bonds_analysis/scanlan_case_analysis.txt
        *The average functionality of the effective crosslinks does NOT include the 
        bridging crosslinks.
        
        The following crosslink attributes in file 
            contract_bonds_analysis/crosslink_properties.txt:
        - types (Tetra-S or Carbon)
        - effective functionality (number of effective bonds on it)
        - scanlan classification (effective, bridging, ineffective)
    """
    project = signac.get_project()
    direc = project.fn('') + 'workspace/'

    # for job_id in os.listdir(direc):
    jdir = direc + job_id
    input_file = direc + job_id + "/" + file
    if not os.path.isfile(input_file):
        raise FileNotFoundError("File not found")
        exit()
    # open the contract_bonds gsd file's trajectory (list of frames)
    trajectory = gsd.hoomd.open(input_file)

    frame = trajectory[-1]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]

    # find the effective bond cutoff based on the maximum bond length of known ineffective bonds
    HS_SH_id = frame.bonds.types.index('-SH')
    C_C_double_id = frame.bonds.types.index('CC_Double')
    known_ineffective_bonds = frame.bonds.group[np.where((frame.bonds.typeid == HS_SH_id) | (frame.bonds.typeid == C_C_double_id))[0]]
    max_ineffective_length = 0
    for bond in known_ineffective_bonds:
        length = dist_pbc(frame.particles.position[bond[0]],\
                          frame.particles.position[bond[1]],\
                          Box=frame.configuration.box[0:3])
        if length > max_ineffective_length:
            max_ineffective_length = length
        # Use the maximum ineffective bond length as the cutoff
    
    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        
    # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
    # bonds in biggest cluster:
    gel_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                gel_bonds.append(each)
                continue
    # make cluster of biggest graph
    gel_G = nx.Graph()
    gel_G.add_edges_from(gel_bonds)

    # categorize crosslinks as effective, bridging, or ineffective
    effective_crosslinks = []
    bridging_crosslinks = []
    ineffective_crosslinks = []
    crosslink_beads = [x for  x in gel_G.nodes() if gel_G.degree(x) >= 3]
    functionalities = []
    pd_types = []
    pd_effective_functionalities = []
    pd_scanlan_classifications = []
    error_crosslinks = []
    for x in crosslink_beads:
        #count how many effective bonds are on it
        count = 0
        for bond in gel_G.edges(x):
            if dist_pbc(frame.particles.position[bond[0]],\
                        frame.particles.position[bond[1]],\
                        Box=frame.configuration.box[0:3]) > max_ineffective_length:
                count += 1
        if count == 0:
            ineffective_crosslinks.append(x)
        elif count == 1:
            print("ERROR: Impossible crosslink found:",x)
            print("This crosslink has only one effective bond on it. This is impossible.")
            print("Bond lengths: ",)
            for bond in gel_G.edges(x):
                print(dist_pbc(frame.particles.position[bond[0]],\
                        frame.particles.position[bond[1]],\
                        Box=frame.configuration.box[0:3]))
            error_crosslinks.append(x)
            print("For now, this crosslink will be considered ineffective.")
            ineffective_crosslinks.append(x)
            # exit()
        elif count == 2:
            bridging_crosslinks.append(x)
        elif count >= 3:
            effective_crosslinks.append(x)
            functionalities.append(count)

        # Save crosslink properties to lists for data frame
        match frame.particles.typeid[x]:
            case 5:
                pd_types.append("Carbon")
            case 6:
                pd_types.append("Tetra-S")
            case _:
                print("ERROR: Unknown crosslink type found:",frame.particles.typeid[x], "at bead:",x)
                raise ValueError("Unknown crosslink type")
        pd_effective_functionalities.append(count)
        if count == 2:
            pd_scanlan_classifications.append("bridging")
        elif count >= 3:
            pd_scanlan_classifications.append("effective")
        else:
            pd_scanlan_classifications.append("ineffective")

    print("Number of error crosslinks:", len(error_crosslinks))
    print("number of crosslinks:",len(crosslink_beads))
    print("effective crosslinks:",len(effective_crosslinks)," ",len(effective_crosslinks)/len(crosslink_beads))
    print("bridging crosslinks:",len(bridging_crosslinks)," ",len(bridging_crosslinks)/len(crosslink_beads))
    print("ineffective crosslinks:",len(ineffective_crosslinks)," ",len(ineffective_crosslinks)/len(crosslink_beads))
    print("average functionality of effective crosslinks:",np.average(functionalities))

    # save the data
    if file == "contract_bonds.gsd":
        w_direc = jdir + "/contract_bonds_analysis/"
    else:
        w_direc = jdir + "/percolation_contract_bonds_analysis/"
    # save the data
    with open(w_direc + "/scanlan_case_analysis.txt", "w") as f:
        f.write("n_crosslinkers\tn_effective_crosslinkers\tn_bridging_crosslinkers\tn_ineffective_crosslinkers\tavg_effective_Functionality\n")
        f.write(f"{len(crosslink_beads)}\t{len(effective_crosslinks)}\t{len(bridging_crosslinks)}\t{len(ineffective_crosslinks)}\t{np.average(functionalities)}\n")

    # make a data frame of the crosslink properties
    # crosslink_df = pd.DataFrame({
    #     'crosslink_id': crosslink_beads,
    #     'type': pd_types,
    #     'effective_functionality': pd_effective_functionalities,
    #     'scanlan_classification': pd_scanlan_classifications
    # })

    # save the data to a txt file
    with open(w_direc + "/crosslink_properties.txt", "w") as f:
        f.write("crosslink_id\ttype\teffective_functionality\tscanlan_classification\n")
        for i in range(len(crosslink_beads)):
            f.write(f"{crosslink_beads[i]}\t{pd_types[i]}\t{pd_effective_functionalities[i]}\t{pd_scanlan_classifications[i]}\n")

def crosslinker_heterogeneity_by_VV(job_id):
    """
    Analyze the final network to check for voronoi volume heterogeneity of EFFECTIVE 
    crosslinkers. (requires the scanlan_case_analysis_on_contract_bonds function to be 
    run first)
    
    Output:
        - voronoi volume distribution of all elastically effective crosslinkers
        - voronoi volume distribution of elastically effective tetrathiol crosslinkers
        - voronoi volume distribution of elastically effective ene crosslinkers
        # These were discarded due to the complexity necessary for unknown returns
        # - distribution of distance from center of cell to all effective crosslinkers 
        #         ("Centrality")
        # - distribution of distance from center of cell to effective tetrathiol crosslinkers
        # - distribution of distance from center of cell to effective ene crosslinkers
    """
    # get the project and job
    project = signac.get_project()
    job = project.open_job(id=job_id)

    #read the crosslink properties from the file
    crosslink_properties = np.genfromtxt(job.fn('contract_bonds_analysis/crosslink_properties.txt'),skip_header=1,\
        dtype='str', delimiter='\t')
    # crosslink_id	type	effective_functionality	scanlan_classification

    # read the frame from the gsd file
    input_file = job.fn("polymerize.gsd")
    frame = gsd.hoomd.open(input_file)[-1]
    pos = frame.particles.position

    # get the crosslinker positions based on the crosslink properties
    effective_ene_crosslinkers_pos = []
    effective_thiol_crosslinkers_pos = []
    for i in range(len(crosslink_properties)):
        if crosslink_properties[i][1] == "Tetra-S" and crosslink_properties[i][3] == "effective":
            effective_thiol_crosslinkers_pos.append(pos[int(crosslink_properties[i][0])])
        elif crosslink_properties[i][1] == "Carbon" and crosslink_properties[i][3] == "effective":
            effective_ene_crosslinkers_pos.append(pos[int(crosslink_properties[i][0])])
    effective_ene_crosslinkers_pos = np.array(effective_ene_crosslinkers_pos).reshape(-1,3)
    effective_thiol_crosslinkers_pos = np.array(effective_thiol_crosslinkers_pos).reshape(-1,3)

    # set up for voronoi volume calculation
    import freud
    box = freud.box.Box(Lx=frame.configuration.box[0], \
                        Ly=frame.configuration.box[1], \
                        Lz=frame.configuration.box[2])
    voro = freud.locality.Voronoi()
    # get the voronoi volume of ALL the effective crosslinkers
    if job.sp['crosslinker_percent'] > 0 and job.sp['chain_side_reaction_probability'] > 0.0:
        q_pos = np.concatenate((effective_ene_crosslinkers_pos, effective_thiol_crosslinkers_pos), axis=0)
        voro.compute((box, q_pos))
        volumes = voro.volumes
        polytopes = voro.polytopes
        np.savetxt(job.fn('contract_bonds_analysis/voronoi_volumes_all.txt'), volumes, fmt='%f')
    else:
        np.savetxt(job.fn('contract_bonds_analysis/voronoi_volumes_all.txt'), [], fmt='%f')
    # get the voronoi volume of the effective tetrathiol crosslinkers
    if job.sp['crosslinker_percent'] > 0:
        q_pos = effective_thiol_crosslinkers_pos
        voro.compute((box, q_pos))
        volumes = voro.volumes
        polytopes = voro.polytopes
        np.savetxt(job.fn('contract_bonds_analysis/voronoi_volumes_thiol.txt'), volumes, fmt='%f')
    else:
        np.savetxt(job.fn('contract_bonds_analysis/voronoi_volumes_thiol.txt'), [], fmt='%f')
    # get the voronoi volume of the effective ene crosslinkers
    if job.sp['chain_side_reaction_probability'] > 0.0:
        q_pos = effective_ene_crosslinkers_pos
        voro.compute((box, q_pos))
        volumes = voro.volumes
        polytopes = voro.polytopes
        np.savetxt(job.fn('contract_bonds_analysis/voronoi_volumes_ene.txt'), volumes, fmt='%f')
    else:
        np.savetxt(job.fn('contract_bonds_analysis/voronoi_volumes_ene.txt'), [], fmt='%f')

    # centrality = np.zeros(len(volumes))
    # for i in range(len(polytopes)):
    #     CM = calc_pbc_CM(polytopes[i],box)
    #     # calculate the distance from the center of mass to the points
    #     centrality[i] = np.sqrt()

def primary_loop_types_analysis(job_id,analyzed_frame = -1):
    """
    Analyze the primary loops to categorize them into different types based on 
    their formation mechanism.
    """
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

    # if not analyzing the final frame
    if analyzed_frame != -1:
        from scripts.conversion import get_frame_at_95_conversion
        analyzed_frame = int(get_frame_at_95_conversion(job_id))
        print(analyzed_frame)


    frame = trajectory[analyzed_frame]
    dummy_id = len(frame.bonds.types)-1
    bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
    type_ids = frame.particles.typeid


    Ntotal = frame.particles.N

    # make a graph - each bond is an edge
    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
        
    # strand lengths (in entire system. If you only want the ones in the gel, you need to look at Gcc[0] only)
    # bonds in biggest cluster:
    new_bonds = []
    for each in bonds:
        for i in each:
            if i in Gcc[0]:
                new_bonds.append(each)
                continue
    # make cluster of biggest graph
    G2 = nx.Graph()
    G2.add_edges_from(new_bonds)

    print("calculating dangling ends")
    dangling_end_strands = []
    dangling_end_lengths = []

    no_danglingEnds_graph, dangling_ends = remove_dangling_ends(G2)
    # print('no_danglingEnds_graph edges', no_danglingEnds_graph.edges())
    # print('no_danglingEnds_graph nodes', no_danglingEnds_graph.nodes())
    dangling_end_graph = G2.copy()
    for x in G2.nodes():
        if x not in dangling_ends:
            dangling_end_graph.remove_node(x)
    # print("dangling ends:",dangling_end_graph.nodes())
    # print("dangling end edges:",dangling_end_graph.edges())

    for dangling_strand in nx.connected_components(dangling_end_graph):
        dangling_strand = list(dangling_strand)
        # print(dangling_strand)
        dangling_end_strands.append(dangling_strand)
        dangling_end_lengths.append(len(dangling_strand))
    unique_dangling_end_lengths,count = np.unique(dangling_end_lengths,return_counts=True)
    dangling_end_data = list(zip(unique_dangling_end_lengths,count))
    # print(dangling_end_data)

    print("calculating loops",  flush=True)
    start = datetime.now()
    crosslink_graph = remove_extenders(no_danglingEnds_graph)
    # print('crosslink_graph edges', crosslink_graph.edges())
    # print('crosslink_graph nodes', crosslink_graph.nodes())
    print('start loop analysis:',start,  flush=True)
    loops = list(nx.simple_cycles(crosslink_graph,length_bound=4))
    # print("found loops",loops)
    print('end loop analysis:',datetime.now(),  flush=True)
    print('duration:',datetime.now()-start,  flush=True)
    # lengths_loops = np.array([len(l) for l in loops])
    loop_types = np.array([classify_loop_type_simple(l) for l in loops])
    # print("loops[loop_types=='primary']:",loops[np.where(loop_types=='primary')[0]])
    tetrathiol_primary_loops = []
    ene_primary_loops = []
    for loop in loops:
        if len(loop) == 1:
            if type_ids[loop[0]] == 6:
                print("found primary loop of type C:", loop)
                tetrathiol_primary_loops.append(loop)
            elif type_ids[loop[0]] == 5:
                print("found primary loop of type Carbon:", loop)
                ene_primary_loops.append(loop)
            else:
                print("found primary loop of unknown type:",type_ids[loop[0]])
                raise ValueError("Unknown primary loop type")
    primary_edges, secondary_edges, tertiary_edges, quaternary_edges = categorize_loop_edges(loops, crosslink_graph)
    primary_lengths, secondary_lengths, tertiary_lengths, quaternary_lengths = get_length_categorized_loops(primary_edges,secondary_edges,tertiary_edges,quaternary_edges)
    # print("primary loop edges and lengths:",primary_edges[0:100], primary_lengths[0:100])
    print("# of primary loops:",len(primary_lengths))
    print('# of primary loops of type C', len(tetrathiol_primary_loops))
    print('# of primary loops of type Carbon', len(ene_primary_loops))

    # save the data
    output_file_name = jdir + "/primary_loop_types.txt"

    with open(output_file_name, "w") as f:
        f.write("N_Tetrathiol_primary_loops\tN_ene_primary_loops\n")
        line = str(len(tetrathiol_primary_loops)) + "\t" + str(len(ene_primary_loops)) + "\n"
        f.write(line)
    
def main(job_id):
    xlink_rdf_analysis(job_id)
    exit()
    defect_analysis(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 network_properties.py job_id")
        sys.exit(1)
    main(sys.argv[1])
