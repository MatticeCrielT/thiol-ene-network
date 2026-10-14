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
import numpy as np

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
    try:
        conversion_data = np.genfromtxt(job.fn('trajectory_conversion.txt'), comments="#", delimiter=" ")
    except FileNotFoundError:
        return False
    cutoff = 0
    if float(conversion_data[-1,1]) - float(conversion_data[-5,1]) <= cutoff:
        return True
    else:
        return False

@MyProject.label
def reacting(job):
    return job.isfile("polymerize.gsd")

@MyProject.label
def deformed(job):
    return job.isfile("deform.gsd") and os.path.getsize(job.fn('deform.log')) > 100

@MyProject.label
def deformed_gel(job):
    return job.isfile("deform_gel.gsd") and os.path.getsize(job.fn('deform_gel.log')) > 100
    
@MyProject.label
def contracted_bonds(job):
    return job.isfile("contract_bonds.gsd")

@MyProject.label
def diffused(job):
    """
    Only two jobs need to run this operation
    """
    if job.sp['chain_side_reaction_probability'] == 0 and \
            job.sp['crosslinker_percent'] == 50 and \
            job.sp['replica_index'] == 0:
        # only two jobs need to run, one for each angle constant
        return job.isfile("diffuse.gsd")
    else:
        return True

#-----------------------
# Analysis labels
#-----------------------

# def percolation_analyzed(job):
#     return job.isfile('percolation_data.txt')
@MyProject.label
def basic_strand_analyzed(job):
    return job.isfile('strand_sizes.txt')

@MyProject.label
def defect_analyzed(job):
    return job.isfile('loop_counts.txt')

@MyProject.label
def xlink_rdf_analyzed(job):
    #If the system cannot form any crosslinks, don't even try to make an RDF of them
    if job.sp['chain_side_reaction_probability'] != 0 :
        one = job.isfile('ene_ene_rdf.txt')
    else:
        one = True
    if job.sp['crosslinker_percent'] != 0:
        two = job.isfile('thiol_thiol_rdf.txt')
    else:
        two = True
    if job.sp['crosslinker_percent'] != 0 and job.sp['chain_side_reaction_probability'] != 0:
        three = job.isfile('thiol_ene_rdf.txt')
    else:
        three = True
    return one and two and three

# @MyProject.label
# def stress_strain_analyzed(job):
#     return job.isfile('stress_strain.txt')

@MyProject.label
def stress_strain_gel_analyzed(job):
    return job.isfile('stress_strain_gel.txt')


@MyProject.label
def modulus_analyzed(job):
    return job.isfile('youngs_modulus.txt')

@MyProject.label
def modulus_gel_analyzed(job):
    return job.isfile('youngs_modulus_gel.txt')

@MyProject.label
def gelation_conversion_2(job):
    try:
        return job.doc['gelation_conversion_2_2'][0] >= 0.0
    except KeyError:
        return False

@MyProject.label
def crosslinking_density_calculated(job):
    try:
        return False #job.doc['crosslinking_density'] >= 0.0
    except KeyError:
        return False
    
@MyProject.label
def contract_bonds_analyzed(job):
    if job.sp['crosslinker_percent'] == 0 and job.sp['chain_side_reaction_probability'] == 0:
        return True
    else: 
        return job.isfile('contract_bonds_analysis/effective_strand_hist.txt') and \
        job.isfile('contract_bonds_analysis/ineffective_strand_hist.txt') and \
        job.isfile('contract_bonds_analysis/overall.txt') and \
        job.isfile('contract_bonds_analysis/scanlan_case_analysis.txt') and \
        job.isfile('contract_bonds_analysis/crosslink_properties.txt')

@MyProject.label
def vv_analyzed(job):
    return job.isfile('contract_bonds_analysis/voronoi_volumes_all.txt') and \
            job.isfile('contract_bonds_analysis/voronoi_volumes_thiol.txt') and \
            job.isfile('contract_bonds_analysis/voronoi_volumes_ene.txt')

@MyProject.label
def diffusion_analyzed(job):
    """
    Only two jobs need to run this operation
    """
    if job.sp['chain_side_reaction_probability'] == 0 and \
            job.sp['crosslinker_percent'] == 75 and \
            job.sp['replica_index'] == 0:
        return job.isfile(f'diffusion_length_at_delta_timesteps_{job.sp["polymerize_period"]}.txt')
    else:
        return True

@MyProject.label
def primary_loop_types_analyzed(job):
    return job.isfile('primary_loop_types.txt')

#-----------------------
# Simulation operations
#-----------------------
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

@MyProject.pre(reacted)
@MyProject.post(deformed)
@MyProject.operation
def deform(job):
    import scripts.deform
    scripts.deform.main(job,job.fn('polymerize.gsd'),-1)
    print('deformed: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(deformed_gel)
@MyProject.operation
def deform_gel(job):
    import scripts.deform_gel;
    scripts.deform_gel.main(job,job.fn('polymerize.gsd'),-1)
    print('deformed gel: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(contracted_bonds)
@MyProject.operation
def contract_bonds(job):
    from scripts.simulate import Simulator
    sinit = Simulator(job)
    sinit.contract_bonds()
    print("contracted bonds: ",job.id)

@MyProject.pre(equilibrated)
@MyProject.post(diffused)
@MyProject.operation
def diffuse(job):
    from scripts.simulate import Simulator
    sinit = Simulator(job)
    sinit.diffuse()
    print("diffused: ",job.id)

#-----------------------
# Analysis operations
#-----------------------
# @MyProject.pre(reacting)
# @MyProject.post(trajectory_conversion_analyzed)
# @MyProject.operation
# def analyze_trajectory_conversion(job):
#     os.system('python3 ./scripts/trajectory_conversion.py')

@MyProject.pre(reacting)
@MyProject.post(reacted)
@MyProject.operation
def calculate_conversion(job):
    import scripts.conversion
    print("now conducting trajectory conversion analysis")
    scripts.conversion.calculate_conversion(job.id)
    print('analyzed the trajectory conversion of ', job.id)

@MyProject.pre(reacted)
@MyProject.post(defect_analyzed)
@MyProject.operation
def defect_analysis(job):
    import scripts.network_properties
    scripts.network_properties.defect_analysis(job.id)
    print('analyzed defectivity: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(xlink_rdf_analyzed)
@MyProject.operation
def xlink_rdf_analysis(job):
    import scripts.network_properties
    scripts.network_properties.xlink_rdf_analysis(job.id)
    print('analyzed crosslink rdfs: ',job.id)


@MyProject.pre(reacted)
@MyProject.post(basic_strand_analyzed)
@MyProject.operation
def strand_analysis(job):
    import scripts.network_properties
    scripts.network_properties.strand_lengths_analysis(job.id)
    print('analyzed strands: ',job.id)


# @MyProject.pre(deformed)
# @MyProject.post(stress_strain_analyzed)
# @MyProject.operation
# def deformation_analysis(job):
#     import scripts.plot_deforms
#     scripts.plot_deforms.analyze_stress_strain(job.id)
#     print('analyzed stress-strain: ',job.id)

@MyProject.pre(deformed_gel)
@MyProject.post(stress_strain_gel_analyzed)
@MyProject.operation
def deformation_gel_analysis(job):
    import scripts.plot_deforms
    scripts.plot_deforms.analyze_stress_strain(job.id)
    print('analyzed stress-strain: ',job.id)
    

# @MyProject.pre(stress_strain_analyzed)
# @MyProject.post(modulus_analyzed)
# @MyProject.operation
# def modulus_analysis(job):
#     import scripts.plot_deforms
#     scripts.plot_deforms.analyze_modulus(job.id)
#     print('analyzed modulus: ',job.id)

@MyProject.pre(stress_strain_gel_analyzed)
@MyProject.post(modulus_gel_analyzed)
@MyProject.operation
def modulus_gel_analysis(job):
    import scripts.plot_deforms
    scripts.plot_deforms.analyze_modulus(job.id)
    print('analyzed modulus: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(gelation_conversion_2)
@MyProject.operation
def gelation_2_analysis(job):
    import scripts.gelation2
    scripts.gelation2.gelation_analysis_via_reducedMw(job.id)
    print('analyzed gelation via reduced Mw method: ',job.id)

@MyProject.pre(contracted_bonds)
@MyProject.post(contract_bonds_analyzed)
@MyProject.operation
def contract_bonds_analysis(job):
    import scripts.network_properties
    # scripts.network_properties.calculate_crosslinking_density(job.id)
    scripts.network_properties.contract_bonds_analysis(job.id)
    scripts.network_properties.scanlan_case_analysis_on_contract_bonds(job.id)
    print('analyzed contract bonds: ',job.id)

@MyProject.pre(contract_bonds_analyzed)
@MyProject.post(crosslinking_density_calculated)
@MyProject.operation
def calculate_cld(job):
    import scripts.network_properties
    scripts.network_properties.calculate_crosslinking_density(job.id)
    print('calculated cld: ',job.id)

@MyProject.pre(reacted and contract_bonds_analyzed)
@MyProject.post(vv_analyzed)
@MyProject.operation
def vv_analysis(job):
    import scripts.network_properties
    scripts.network_properties.crosslinker_heterogeneity_by_VV(job.id)
    print('analyzed voronoi volumes of job: ',job.id)

@MyProject.pre(diffused)
@MyProject.post(diffusion_analyzed)
@MyProject.operation
def diffusion_analysis(job):
    from scripts.diffusion import DiffusionAnalyzer
    analyzer = DiffusionAnalyzer(job, dt=0.005)
    # if not job.isfile('MSD.txt'):
    #     print('MSD.txt does not exist, calculating MSD')
    #     analyzer.calculate_MSD()
    # else:
    #     print('MSD.txt already exists, skipping MSD calculation')
    # analyzer.calculate_diffusion_coeff()
    analyzer.calculate_diffusion_length_at_delta_timestep(timesteps=job.sp['polymerize_period'])
    print('analyzed diffusion: ',job.id)

'''
assess percolated frame
'''

@MyProject.label
def percolation_analyzed(job):
    return job.isfile("percolation_data.txt")

@MyProject.label
def percolation_frame_recorded(job):
    try:
        return job.doc['percolation_frame'] > 0
    except KeyError:
        return False

@MyProject.pre(diffused)
@MyProject.post(percolation_frame_recorded)
@MyProject.operation
def get_percolation_frame(job):
    import scripts.conversion
    scripts.conversion.get_percolation_frame(job.id)

@MyProject.label
def percolation_contracted_bonds(job):
    return job.isfile("percolation_contract_bonds.gsd")

@MyProject.label
def percolation_contracted_bonds_analyzed(job):
    if job.sp['crosslinker_percent'] == 0 and job.sp['chain_side_reaction_probability'] == 0:
        return True
    else: 
        return job.isfile('percolation_contract_bonds_analysis/effective_strand_hist.txt') and \
        job.isfile('percolation_contract_bonds_analysis/ineffective_strand_hist.txt') and \
        job.isfile('percolation_contract_bonds_analysis/overall.txt') and \
        job.isfile('percolation_contract_bonds_analysis/scanlan_case_analysis.txt') and \
        job.isfile('percolation_contract_bonds_analysis/crosslink_properties.txt')

@MyProject.pre(percolation_frame_recorded)
@MyProject.post(percolation_contracted_bonds)
@MyProject.operation
def percolation_contract_bonds(job):
    from scripts.simulate import Simulator
    sinit = Simulator(job, contract_frame = int(job.doc['percolation_frame']))
    sinit.contract_bonds(contract_frame = int(job.doc['percolation_frame']))
    print("contracted bonds of percolation frame: ",job.id)


@MyProject.pre(percolation_contracted_bonds)
@MyProject.post(percolation_contracted_bonds_analyzed)
@MyProject.operation
def percolation_contract_bonds_analysis(job):
    import scripts.network_properties
    # scripts.network_properties.calculate_crosslinking_density(job.id)
    scripts.network_properties.contract_bonds_analysis(job.id, file="percolation_contract_bonds.gsd")
    scripts.network_properties.scanlan_case_analysis_on_contract_bonds(job.id, file="percolation_contract_bonds.gsd")
    print('analyzed contract bonds: ',job.id)


@MyProject.label
def percolation_crosslinking_density_calculated(job):
    try:
        return job.doc['perc_crosslinking_density'] >= 0.0
    except KeyError:
        return False

@MyProject.pre(percolation_contracted_bonds_analyzed)
@MyProject.post(percolation_crosslinking_density_calculated)
@MyProject.operation
def calculate_percolation_cld(job):
    import scripts.network_properties
    scripts.network_properties.calculate_crosslinking_density(job.id, frame=int(job.doc['percolation_frame']))
    print('calculated cld: ',job.id)

@MyProject.label
def conv95_defect_analyzed(job):
    return job.isfile("loop_counts_conversion95.txt")

@MyProject.label
def conv95_strands_analyzed(job):
    return job.isfile("conv95_strand_sizes.txt")

@MyProject.pre(reacted)
@MyProject.post(conv95_defect_analyzed)
@MyProject.operation
def conv95_defect_analysis(job):
    import scripts.network_properties
    scripts.network_properties.defect_analysis(job.id, defect_frame=95)
    print('analyzed defectivity: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(conv95_strands_analyzed)
@MyProject.operation
def conv95_strand_analysis(job):
    import scripts.network_properties
    scripts.network_properties.strand_lengths_analysis(job.id, defect_frame=95)
    print('analyzed strands: ',job.id)

@MyProject.pre(reacted)
@MyProject.post(primary_loop_types_analyzed)
@MyProject.operation
def primary_loop_types_analysis(job):
    import scripts.network_properties
    scripts.network_properties.primary_loop_types_analysis(job.id)
    print('analyzed primary loop types: ',job.id)

if __name__ == "__main__":
    MyProject().main()
