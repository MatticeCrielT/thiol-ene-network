#!/opt/miniconda3/bin/python
"""Initialize the project's data space."""
import numpy as np
import signac
import itertools
from scripts.system import molecule_counts


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
        "replica_index": [0,1,2,3,4],
        "density": [0.81],
        "temperature":[1],
        "crosslinker_percent":[50],#[0, 12.5, 25, 37.5, 50, 62.5, 75, 87.5, 100],
        "N_monomers":[10000],
        "monomer_size":[0],     # spacer beads between the two C=C of the di-ene (0 = E=E-E=E)
        "extender_size":[0],    # spacer beads between the two SH of the dithiol (0 = S-S)
        "radical_percent":[1.0],
        "chain_side_reaction_probability": [0.0,0.15],#0.0,0.5,0.6,0.9],
        "chain_transfer_probability": [0.1],
        "thiol_reaction_probability": [1],
        "polymerize_period": [50,500],#100
        "r_cut_reaction":[1.1],
        "angle_constant":[100],#[10,100],
        "polymerization_method":["cpu_local_snapshot"],
        }

    for sp in grid(statepoint_grid):
        # check if the sp is for one of the non-standard polymerization periods. If so, limit the other statepoints
        if sp.get("polymerize_period") != 100:
            if sp.get("replica_index") != 0 or sp.get("crosslinker_percent") != 50 or sp.get("angle_constant") != 100:
                continue

        # open the job and initialize
        job = project.open_job(sp).init()
        # same counting as the simulation itself (scripts/system.py)
        counts = molecule_counts(job.sp['N_monomers'], job.sp['crosslinker_percent'],
                                 job.sp['monomer_size'], job.sp['extender_size'])
        N_crosslinker = counts['N_crosslinkers']
        N_extenders = counts['N_extenders']
        N_particles = counts['N_particles']

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
