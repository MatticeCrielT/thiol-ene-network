import signac
import numpy as np
import matplotlib.pyplot as plt
import gsd.hoomd
import networkx as nx

project = signac.get_project()

# find the jobs with [TE,MM] mechanism and Xthiol [12.5, 100]
mechanisms = [0.0,0.15]
crosslinker_percents = [12.5, 100]

for mechanism in mechanisms:
    for crosslinker_percent in crosslinker_percents:
        job_list = project.find_jobs({'chain_side_reaction_probability':mechanism, 'crosslinker_percent':crosslinker_percent,\
                                      'angle_constant':100.0, 'polymerize_period':100})
        print(f"Mechanism: {mechanism}, Crosslinker Percent: {crosslinker_percent}, Number of jobs: {len(job_list)}")

        n_replicates = 1
        # n_replicates = len(job_list)
        avg_n_thiol_crosslinkers = 0
        avg_n_ene_crosslinkers = 0

        for job in job_list:
            if job.sp["replica_index"] != 0:
                continue
            # print(f"Processing job {job.id}")
            # read percolation data to find the frame 3D percolation occurs
            percolation_data = np.loadtxt(job.fn('percolation_data.txt'),skiprows=1)
            percolated_rows = percolation_data[percolation_data[:,1] == 3]
            percolated_frame = int(percolated_rows[0,0])  # first frame where 3D percolation occurs
            print(f"Job {job.id}: 3D percolation occurs at frame {percolated_frame}")

            # load the percolated frame of the gsd file
            with gsd.hoomd.open(name=job.fn('polymerize.gsd'), mode='r') as traj:
                frame = traj[percolated_frame]
                typeids = frame.particles.typeid
                bonds = frame.bonds.group

            # build graph from bonds
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

            # down select to only the nodes with functionality > 2
            crosslinker_nodes = [n for n, d in gel_G.degree() if d > 2]
            # print(f"Job {job.id}: Number of nodes in gel went from: {len(gel_G.nodes())} to {len(crosslinker_nodes)}")

            n_thiol_crosslinkers = 0
            n_ene_crosslinkers = 0
            thiol_crosslinker_typeid = frame.particles.types.index('C')
            ene_crosslinker_typeid = frame.particles.types.index('Carbon')
            for node in crosslinker_nodes:
                if typeids[node] == thiol_crosslinker_typeid:
                    n_thiol_crosslinkers += 1
                elif typeids[node] == ene_crosslinker_typeid:
                    n_ene_crosslinkers += 1
                else:
                    print("Unknown crosslinker typeid:", typeids[node])
                    exit()
            # print(f"Job {job.id}: Number of thiol crosslinkers in gel: {n_thiol_crosslinkers}, Number of ene crosslinkers in gel: {n_ene_crosslinkers}")
            avg_n_thiol_crosslinkers += n_thiol_crosslinkers / n_replicates
            avg_n_ene_crosslinkers += n_ene_crosslinkers / n_replicates
        
        # jobs have all been iterated over
        print(f"Mechanism: {mechanism}, Crosslinker Percent: {crosslinker_percent}, Average number of thiol crosslinkers in gel: {avg_n_thiol_crosslinkers}, Average number of ene crosslinkers in gel: {avg_n_ene_crosslinkers}")