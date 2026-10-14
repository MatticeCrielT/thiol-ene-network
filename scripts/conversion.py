import sys, os, re
import numpy as np
import signac
import json
import gsd, gsd.hoomd 
from collections import defaultdict, Counter
from scripts.extract import connected_components

'''
VIA GSD.HOOMD
for every frame of the polymerization trajectory
by going through each string of connected atoms, calculate conversion, 
the largest molecule by bead size, and the average molecule bead size 
'''
def conversion_molecule_sizes(N_0, trajectory, dummy_id):
    traj_conversion = []
    thiol_ene_conversion = []
    # through each frame of the trajectory
    for i,frame in enumerate(trajectory):
        # find the types of each atom in the 
        type_map = list(frame.particles.typeid)
        if i == 0:
            init_functionalities = type_map.count(1)/2

        bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
        connections = connected_components(bonds)
        strand_lengths = []
        n_molecules = 0

        for strand in connections:
            strand_lengths.append(len(strand))
            n_molecules +=1
        
        counts = dict(Counter(strand_lengths))
        sum_num = 0
        sum_den = 0
        for l in list(counts.keys()):
            sum_num += float(l)**2 * counts[l]
            sum_den += float(l) * counts[l]

        m_w = sum_num/sum_den

        # calculate the extent of the reaction based off of molecule formation
        # also store the largest molecule size and average molecule size
        extent_of_reaction = (N_0 - n_molecules)/N_0
        traj_conversion.append([i,extent_of_reaction, np.max(strand_lengths), np.average(strand_lengths), m_w])

        # calculate functionality conversion
        thiols = type_map.count(0) 
        enes = type_map.count(1)/2
        thiol_ene_conversion.append([i, 1 - thiols/init_functionalities, 1 - enes/init_functionalities])

    return np.array(traj_conversion), np.array(thiol_ene_conversion)

project = signac.get_project()


'''
runs on polymerize file in the workspace, adds on to trajectory conversions that are already recorded
calculates conversion by 
(1) molecule -> trajectory_conversion.txt
(2) functionality -> trajectory_thiol_ene_conversion.txt
'''
def calculate_conversion(job_id):
    # find workspace directory
    direc = project.fn('') + "workspace/"
    jdir = direc + job_id
    polymerize_gsd = jdir + "/polymerize.gsd"
    if not os.path.isfile(polymerize_gsd):
         return "oops no polymer file"
    print(polymerize_gsd)

    trajectory = gsd.hoomd.open(polymerize_gsd)
    print(len(trajectory))
    # get the dummy_bond_id and number of bonds from the initial frame
    initial_frame = trajectory[0]
    dummy_id = len(initial_frame.bonds.types)-1
    print(dummy_id)
    bonds = initial_frame.bonds.group[initial_frame.bonds.typeid!=dummy_id] # dummy type bond 
    # find all strands in the box
    all_strands = connected_components(bonds)
    # count the number of molecules
    n_molecules = 0
    for strand in all_strands:
            n_molecules +=1
    N_0 = n_molecules

    # get molecule size based trajectory conversion
    t_c, t_te = conversion_molecule_sizes(N_0, trajectory, dummy_id)
    print(t_c)
    trajectory_conversion_txt = jdir + "/trajectory_conversion.txt"
    thiolene_conversion_txt = jdir + "/trajectory_thiol_ene_conversion.txt"
    
    # molecule conversion
    # if there is already data stored, add onto the old conversion data (do not restart frame numbering)
    # this will allow us to ignore if we truncate our files
    if os.path.isfile(trajectory_conversion_txt):
        data = np.loadtxt(trajectory_conversion_txt, delimiter=' ', skiprows=1)
        # check if we already are fully reacted
        cutoff = 0
        if float(data[-1,1]) - float(data[-5,1]) <= cutoff:
            print("already reacted")
            return None
        # counting if we have reached steady state but not 5 frames:
        old_conversion_values = [d[1] for d in data]
        # ignore past repeated frames
        repeats = int(old_conversion_values.count(data[-1,1]))
        # get the new values and frames
        new_convs = [[float(i[1]),i[2],i[3], i[4]] for i in t_c if i[1] >= data[-1][1]][repeats:]
        t_c_new = np.array([[int(n + 1 + int(data[-1][0]))] + i for n,i in enumerate(new_convs)])
        if t_c_new.size > 0:
            t_c = np.concatenate((data, t_c_new))
    
    # dump the molecule conversions into the txt
    header = "frame conversion largest_molecule average_molecule weight_average_molecule"
    np.savetxt(trajectory_conversion_txt, t_c, header=header, fmt='%i %.16f %.16f %.16f %.16f', comments='')

    # thiol and ene functionality conversion
    # if there is already data stored, add onto the old conversion data (do not restart frame numbering)
    # this will allow us to ignore if we truncate our files
    
    if os.path.isfile(thiolene_conversion_txt):
        data = np.loadtxt(thiolene_conversion_txt, delimiter=' ', skiprows=1)
        # compare ene conversion of last file to current
        new_convs = [[float(i[1]), float(i[2])] for i in t_te if i[2] >= data[-1][2]]
        header = "frame thiol_conv ene_conv"
        if new_convs == []:
            np.savetxt(thiolene_conversion_txt, t_te, header=header, fmt='%i %.16f %.16f', comments='')
        else:
            # add the new conversions to the corrected 
            # counting final repeats:
            new_convs = new_convs[repeats:]
            old_conversion_values = [d[1] for d in data]
            t_te_new = np.array([[int(n + 1 + int(data[-1][0])), i[0], i[1]] for n,i in enumerate(new_convs)])
            if t_te_new.size > 0:
                t_te = np.concatenate((data, t_te_new))
            # dump the thiol-ene conversions into the txt
            np.savetxt(thiolene_conversion_txt, t_te, header=header, fmt='%i %.16f %.16f', comments='')

'''
runs on polymerize file in the workspace, adds on to conversion text files that are already recorded
calculates weight average molecular weight of molecules during the polymerization
records as the fifth column of trajectory_converison.txt
'''
def calculate_mw(job_id):
    # find workspace directory
    direc = project.fn('') + "workspace/"
    jdir = direc + job_id
    polymerize_gsd = jdir + "/polymerize.gsd"
    if not os.path.isfile(polymerize_gsd):
         return "oops no polymer file"
    print(polymerize_gsd)

    trajectory = gsd.hoomd.open(polymerize_gsd)
    print(len(trajectory))

    # get the dummy_bond_id and number of bonds from the initial frame
    initial_frame = trajectory[0]
    dummy_id = len(initial_frame.bonds.types)-1

    bonds = initial_frame.bonds.group[initial_frame.bonds.typeid!=dummy_id] # dummy type bond 

    m_w = []
    for i,frame in enumerate(trajectory):
        # get all the molecules in the system to analyze their lengths
        bonds = frame.bonds.group[frame.bonds.typeid!=dummy_id]
        connections = connected_components(bonds)
        strand_lengths = []
        for strand in connections:
            strand_lengths.append(len(strand))

        # calculate and store weight-average molecular weight
        counts = dict(Counter(strand_lengths))
        sum_num = 0
        sum_den = 0
        for l in list(counts.keys()):
            sum_num += float(l)**2 * counts[l]
            sum_den += float(l) * counts[l]

        m_w.append(sum_num/sum_den)
    
    trajectory_conversion_txt = jdir + "/trajectory_conversion.txt"
    if os.path.isfile(trajectory_conversion_txt):
        data = np.loadtxt(trajectory_conversion_txt, delimiter=' ', skiprows=1)
        new_data = np.array([[d[0],d[1],d[2],d[3], m_w[n]] for n,d in enumerate(data)])
        # dump into conversion file
        header = "frame conversion largest_molecule average_molecule weight_average_molecule"
        np.savetxt(trajectory_conversion_txt, new_data, header=header, fmt='%i %.16f %.16f %.16f %.16f', comments='')
    else:
        print("no conversion file")

def get_frame_at_95_conversion(job_id):
    # find workspace directory
    direc = project.fn('') + "workspace/"
    jdir = direc + job_id
    trajectory_conversion_txt = jdir + "/trajectory_conversion.txt"
    if not os.path.isfile(trajectory_conversion_txt):
        return "oops no conversion file"
    
    data = np.loadtxt(trajectory_conversion_txt, delimiter=' ', skiprows=1)
    # find the frame where conversion is just above 0.95, starting from the end
    for i in np.arange(start=-1, stop=-len(data)-1, step=-1):
        if data[i, 1] <= 0.95:
            return i + 1 + len(data)
    return False

def get_percolation_frame(job_id):
    # find job directory
    project = signac.get_project()
    job = project.open_job(id=job_id)
    # get percolation file
    percolation_txt = job.fn("percolation_data.txt")
    data = np.loadtxt(percolation_txt, delimiter=' ', skiprows=1)
    # first frame in which the largest cluster spans the box in all 3 directions
    for i in data:
        if i[1] > 2:
            job.doc['percolation_frame'] = i[0]
            break

