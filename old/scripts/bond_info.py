import sys, os, re
import numpy as np
import signac
import json
import gsd, gsd.hoomd 
import networkx as nx 
from collections import defaultdict



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

def connected_components(lists):
    """
    merges lists with common elements

    args: bonds from configuration (frame.bonds.group)

    returns: list of connected particles by id

    Useful for finding bonded particles in a configuration, tested for
    linear polymers with consecutive bonds (0-1-2-3-4-5, 6-7-8-9-10,..)
    and non consecutive ids ( 0-5-6-8-10, 1-4-3-2-9,...) but no other
    configuration yet. Works with ints as well as str.

    """
    # create a default dictionary with an element returning a set if nothing in it
    neighbors = defaultdict(list)
    seen = set()

    # for every element in 
    for each in lists:
        for item in each:
            neighbors[item].extend(each)
            neighbors[item].remove(item)
    return neighbors
    """
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
            print("hi")
    """

# our signac project
project = signac.get_project()
# look through each POLYMERIZED job in the project
for job in project: 
    if job.isfile('polymerize.gsd'):

        # only look at our base cases for now...
        with open(job.fn('signac_statepoint.json')) as f:
            statepoint = json.load(f)
            '''
            if not (statepoint["crosslinker_percent"]==50.0 and statepoint["chain_side_reaction_probability"]==0):
                continue
            '''
            
            if (statepoint["chain_side_reaction_probability"]==0):
                continue

        print(job)
        # open the polymerized gsd file's trajectory (list of frames)
        trajectory = gsd.hoomd.open(job.fn('polymerize.gsd'))
        # get the dummy_bond_id and number of bonds from the initial frame
        frame = trajectory[0]
        dummy_id = len(frame.bonds.types)-1
        
        # go through each frame of the trajectory
        trajectory_bond_histograms = {}
        trajectory_alternating_q = {}
        trajectory_radical_numbers = {}
        for i, frame in enumerate(trajectory):
            # get all non dummy bonds, and atom types
            bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
            all_bonds = atom_bonds(bonds)
            type_map = frame.particles.typeid

            # count the number of radicals in the frame
            radical_number = list(type_map).count(2) + list(type_map).count(3)
            # count the atoms with (1,2,3,4, or 5) bonds
            bond_histogram = defaultdict(int, {k:0 for k in (1,2,3,4) })
            # record the atoms that have more than 3 bonds and what they are bonded to
            multibonded_atoms = defaultdict(dict)
            # check that for atoms bonded to two atoms that it is alternating atoms - holds any that are not alt.
            ## this only works for the dubmbells
            alternating_q = {}
            
            
            for atom, neighbors in all_bonds.items():
                bond_histogram[int(len(neighbors))] += 1   
                if len(neighbors) == 2:
                    neighbor_types = [type_map[i] for i in neighbors]
                    # if there are two neighbors check if it's alternating 
                    if not ((0 in neighbor_types or 2 in neighbor_types or 4 in neighbor_types) and (1 in neighbor_types or 3 in neighbor_types or 5 in neighbor_types)):
                        alternating_q[int(atom)] = neighbor_types
                        
                if len(neighbors) > 2:
                    multibonded_atoms[int(len(neighbors))][atom] = neighbors
            # record in overall trajectory lists
            trajectory_bond_histograms[i] = bond_histogram
            #trajectory_alternating_q[i] = alternating_q
            trajectory_radical_numbers[i] = (radical_number, list(type_map).count(2), list(type_map).count(3))

        # write into files
        json_dir = "./workspace/" + str(job) + "/json/"
        if not os.path.isdir(json_dir):
            os.mkdir(json_dir)
        with open(json_dir + 'bond_hitsograms.json', "w") as f :
            json.dump(trajectory_bond_histograms, f)
        #print(trajectory_alternating_q)
        #with open(job.fn('alternating_q.json'), "w") as f:
        #    json.dump(trajectory_alternating_q, f)
        with open(json_dir + 'radical_numbers.json', "w") as f:
            json.dump(trajectory_radical_numbers, f)
        print(list(trajectory_radical_numbers.values())[-1])