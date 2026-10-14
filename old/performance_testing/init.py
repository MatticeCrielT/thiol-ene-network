#!/opt/miniconda3/bin/python
"""Initialize the project's data space."""
import numpy as np
import signac
import itertools


def grid(gridspec):
    """Yields the Cartesian product of a `dict` of iterables.

    The input ``gridspec`` is a dictionary whose keys correspond to
    parameter names. Each key is associated with an iterable of the
    values that parameter could take on. The result is a sequence of
    dictionaries where each dictionary has one of the unique combinations
    of the parameter values.
    """
    for values in itertools.product(*gridspec.values()):
        yield dict(zip(gridspec.keys(), values))


def main():
    project = signac.init_project()
    # define all parameters of interest in a dictionary
    # this script can be changed and simply executed again to extend the parameter space
    # new parameters can also be added to the workspace after creating it

    # change small blob density
    statepoint_grid = {
        "replica_index": [0,1],
        "density": [0.8475],
        "temperature":[0.9],
        "crosslinker_percent":[10.0],#0.0,10.0,25.0,50.0,75.0],
        "N_monomers":[5000,10000],
        "monomer_size":[0],
        "extender_size":[0],
        "radical_percent":[1.0],
        "chain_side_reaction_probability": [0.333],#0.0,0.5,0.6,0.9],
        "chain_transfer_probability": [0.333],
        "thiol_reaction_probability": [0.5],
        # "polymerize_period": [100],#,1,1000000000],
        "polymerize_period": [1,10,100,1000,10000,100000],#,1,1000000000],
        "r_cut_reaction":[1.1],
        "angle_constant":[100],
        # "polymerization_method":["custom_action_GPU","custom_action_GPU_bulk"],
        "polymerization_method":["cpu_local_snapshot","custom_action_GPU","custom_action_CPU"],
        }

    for sp in grid(statepoint_grid):
        # only create the statepoint with 10,000 monomers with polymerize_period 100 and on cpu_local_snapshot
        if(sp.get("N_monomers") == 10000):
            if (sp.get("polymerize_period") != 100) or (sp.get("polymerization_method") != "cpu_local_snapshot") or (sp.get("replica_index") != 0):
                continue

        # open the job and initialize
        job = project.open_job(sp).init()
        if job.sp['monomer_size']>0 and job.sp['extender_size']>0:
            size_monomer =  job.sp['monomer_size']+6 
            size_extender =  job.sp['extender_size']+1 
            reactive_beads = job.sp['N_monomers']*2 
            crosslinker_fraction = job.sp['crosslinker_percent']/100.
            N_crosslinker = int(np.round(reactive_beads*crosslinker_fraction/4.))
            N_extenders =  int(np.round((reactive_beads-N_crosslinker*4)/2.))
            N_particles =  N_crosslinker*13 + N_extenders*size_extender+ job.sp['N_monomers']*size_monomer
        else: 
            size_monomer =  2
            size_extender =  2
            reactive_beads = job.sp['N_monomers']*2 
            crosslinker_fraction = job.sp['crosslinker_percent']/100.
            N_crosslinker = int(np.round(reactive_beads*crosslinker_fraction/4.))
            N_extenders =  int(np.round((reactive_beads-N_crosslinker*4)/2.))
            N_particles =  N_crosslinker*5 + N_extenders*size_extender+ job.sp['N_monomers']*size_monomer

        job.doc['N_extenders'] = N_extenders
        job.doc['N_crosslinkers'] = N_crosslinker
        job.doc['N_total']=  N_particles
        job.doc['L'] =  (N_particles/job.sp['density'])**(1/3.0)
        
        # prevent overwriting of existing jobs 
        try:
            job.doc['reacted_monomers'] 
        except:
            job.doc['reacted_monomers'] = 0

        job.doc['average_polymerization_time'] = 0
        job.doc['average_period_time'] = 0

        try:
            job.doc['avg_polymerization_time']
        except:
            job.doc['avg_polymerization_time'] = 0
        try:
            job.doc['avg_integration_time']
        except:
            job.doc['avg_integration_time'] = 0
        try:
            job.doc['integration_tps']
        except:
            job.doc['integration_tps'] = 0

        print(f"initializing state point with id {job.id}, N_extenders {N_extenders}, N_crosslinkers {N_crosslinker}, N_monomers {job.sp['N_monomers']}")


if __name__ == "__main__":
    main()
