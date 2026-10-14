#!/opt/miniconda3/bin/python
"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python project.py --help
"""
import flow
from flow import FlowProject
import os
import subprocess

class MyProject(FlowProject):
    pass

class CampusCluster(flow.environment.DefaultSlurmEnvironment):
    hostname_pattern = r".*.campuscluster\.illinois\.edu$"
    #template = "campuscluster.sh"

# Labels
@MyProject.label
def equilibrated(job):
    return job.isfile("equi.gsd")

@MyProject.label
def reacted(job):
    # return job.doc['reacted_monomers']>0.9
    if job.isfile("polymerize.gsd"):
        return os.path.getsize(job.fn("polymerize.gsd")) > 500
    else:
        return False

@MyProject.label
def reacting(job):
    return job.isfile("polymerize.gsd")


@MyProject.post(equilibrated)
@MyProject.operation
def equilibrate(job):
    from scripts.simulate import Simulator
    # if job.sp["polymerization_method"] == "cpu_local_snapshot":
    #     from scripts.simulate import Simulator
    # elif job.sp["polymerization_method"] == "custom_action_GPU" or job.sp["polymerization_method"] == "custom_action_CPU" or job.sp["polymerization_method"] == "custom_action_GPU_bulk":
    #     from scripts.customAction.simulate import Simulator
    # else:
    #     print("ERROR polymerization_method not recognized")
    #     exit(2)
    sinit = Simulator(job)
    sinit.equilibrate()
    print("equilibrated",job.id)


@MyProject.pre(equilibrated)
@MyProject.post(reacted)
@MyProject.operation
def polymerize(job):
    from scripts.simulate import Simulator
    # if job.sp["polymerization_method"] == "cpu_local_snapshot":
    #     from scripts.simulate import Simulator
    # elif job.sp["polymerization_method"] == "custom_action_GPU" or job.sp["polymerization_method"] == "custom_action_CPU" or job.sp["polymerization_method"] == "custom_action_GPU_bulk":
    #     from scripts.customAction.simulate import Simulator
    # else:
    #     print("ERROR polymerization_method not recognized")
    #     exit(2)
    sinit = Simulator(job)
    sinit.polymerize()
    print("reacted",job.id)

# Analysis



if __name__ == "__main__":
    MyProject().main()
