#The pressure tensor will NOT be output to the deform.log file. In hoomd 4, it must be 
# accessed via a gsd file. It will be output as pressures.gsd. Use these lines to read:
#   import gsd.hoomd
#   logs = gsd.hoomd.read_log(***Path to pressures.gsd***)
#   pressure_tensors = logs["log/md/compute/ThermodynamicQuantities/pressure_tensor"]

import numpy as np
import sys
import gsd,gsd.hoomd
import sys, os 
import hoomd
# from timeit import default_timer as timer
# from scripts.system import System
from pathlib import Path
import signac

class BoxWriter(hoomd.custom.Action):
    # def __init__(self):
    #     super().attach(simulation)

    def attach(self,simulation):
        super().attach(simulation)

    @hoomd.logging.log(category='scalar')
    def Lx(self):
        snap = self._state.get_snapshot()
        if snap.communicator.rank == 0:
            return snap.configuration.box[0]
        
    @hoomd.logging.log(category='scalar')
    def Ly(self):
        snap = self._state.get_snapshot()
        if snap.communicator.rank == 0:
            return snap.configuration.box[1]
        
    @hoomd.logging.log(category='scalar')
    def Lz(self):
        snap = self._state.get_snapshot()
        if snap.communicator.rank == 0:
            return snap.configuration.box[2]

    def act(self, timestep):
        pass

class NVTCustomBoxVariant(hoomd.variant.box.BoxVariant):
    def __init__(self,initial_box,t_start,t_ramp,final_lam):
        hoomd.variant.box.BoxVariant.__init__(self)
        self.initial_box = np.asarray(initial_box)
        self.deform_time = t_ramp
        self.t_start = t_start
        self.final_lam = final_lam

    def __call__(self, timestep):
        if timestep > self.t_start + self.deform_time:
            lam = self.final_lam
        elif timestep < self.t_start:
            lam = 1
        else:
            lam = (timestep - self.t_start)/self.deform_time*(self.final_lam-1) + 1 # 1 is the initial lam
        return [self.initial_box[0]*lam, self.initial_box[1]/np.sqrt(lam), self.initial_box[2]/np.sqrt(lam), 
                self.initial_box[3], self.initial_box[4], self.initial_box[5]]
'''
These functions were used to run the equilibration before the deformation, but strange and
 violent temperature fluctuations were observed, so now the equilibration is done in the main function
def run_equilibration(job,gsd_file_path, frame_number):
    parent_path = str(Path(gsd_file_path).parent.absolute())

    try:
        device = hoomd.device.GPU()
    except:
        device = hoomd.device.CPU()

    sim = hoomd.Simulation(device=device, seed=1)
    sim.create_state_from_gsd(filename=gsd_file_path, frame=int(frame_number))

    #-----------------------Force fields------------------------
    particle_types = sim.state.particle_types
    cell = hoomd.md.nlist.Cell(buffer=0.4)

    lj = hoomd.md.pair.LJ(nlist=cell)
    lj.params[(particle_types, particle_types)] = dict(epsilon=1.0,sigma=1.0)
    lj.r_cut[particle_types, particle_types] = 2**(1/6) #2.5

    lj.params[(particle_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
    lj.r_cut[particle_types, 'Dummy'] = 0
    lj.mode = 'shift'

    fene = hoomd.md.bond.FENEWCA()
    fene.params[sim.state.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
    fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

    if len(sim.state.angle_types)==0:
        FJ_system=True
    else:
        FJ_system=False

    if FJ_system==False:
        angle_constant = job.sp['angle_constant']
        cosinesq = hoomd.md.angle.CosineSquared()
        cosinesq.params[sim.state.angle_types] = dict(k=angle_constant, t0=np.pi*110/180)# https://www.sciencedirect.com/science/article/pii/S0032386110003642?ref=cra_js_challenge&fr=RR-1
        cosinesq.params['Dummy'] = dict(k=0.0001, t0=np.pi)  # k>0 to make warning go away (should not do anything)
    #-----------------------Force fields------------------------
    #----------------------Integration parameters-------------------------

    kT = 1 # Could read, but just assuming job.sp["temperature"]
    print(f"Could read, but assuming temperature = {kT}")
    dt = 0.005 # Also taken from Brandon's deformation parameter scan
    types_to_integrate =  hoomd.filter.Type(sim.state.particle_types[:-1]) # everything but "Dummy" particles

    integrator = hoomd.md.Integrator(dt=dt)
    
    if FJ_system==False:
        integrator.forces = [lj,fene,cosinesq]
    else:
        integrator.forces = [lj,fene]

    mttk = hoomd.md.methods.thermostats.MTTK(kT=kT,tau=100*dt)
    nvt = hoomd.md.methods.ConstantVolume(filter=types_to_integrate,thermostat=mttk)
    integrator.methods = [nvt]

    sim.operations.integrator = integrator
    #------------------------Configure output---------------------------
    logger = hoomd.logging.Logger(categories=['scalar','string'])
    logger.add(sim, quantities=['timestep','tps'])

    thermo = hoomd.md.compute.ThermodynamicQuantities(types_to_integrate)
    sim.operations.computes.append(thermo)
    logger.add(thermo, ['kinetic_temperature','kinetic_energy', 'potential_energy','pressure','volume'])
    
    box_writer_action = BoxWriter()
    logger.add(box_writer_action, ['Lx','Ly','Lz'])
    box_writer = hoomd.write.CustomWriter(action=box_writer_action,
                                             trigger=hoomd.trigger.Periodic(10))
    sim.operations.writers.append(box_writer)

    log_writer = hoomd.write.Table(output=open(parent_path + '/equi_WCA.log','w'),
                                    trigger=hoomd.trigger.Periodic(100),
                                    logger=logger,)
    sim.operations.writers.append(log_writer)

    logger = hoomd.logging.Logger(categories=['sequence'])
    logger.add(thermo, ['pressure_tensor'])
    gsd_writer = hoomd.write.GSD(filename=parent_path + '/equi_WCA.gsd',
                                    trigger=hoomd.trigger.Periodic(100000),
                                    dynamic=['property','momentum','topology','attribute'],
                                    mode='wb',
                                    logger=logger)
    sim.operations.writers.append(gsd_writer)

    #--------------------------Run the simulation----------------------

    sim.state.thermalize_particle_momenta(kT=kT,filter=types_to_integrate)

    # Run a short equilibration
    sim.run(1e5,write_at_start=True)

    print("Equilibration done! Initial state ready for deformation.")
    del sim

def run_feneOnly_equilibration(job,gsd_file_path, frame_number):
    parent_path = str(Path(gsd_file_path).parent.absolute())

    try:
        device = hoomd.device.GPU()
    except:
        device = hoomd.device.CPU()

    sim = hoomd.Simulation(device=device, seed=1)
    sim.create_state_from_gsd(filename=gsd_file_path, frame=int(frame_number))

    #-----------------------Force fields------------------------
    # particle_types = sim.state.particle_types
    # cell = hoomd.md.nlist.Cell(buffer=0.4)

    # lj = hoomd.md.pair.LJ(nlist=cell)
    # lj.params[(particle_types, particle_types)] = dict(epsilon=1.0,sigma=1.0)
    # lj.r_cut[particle_types, particle_types] = 2**(1/6) #2.5

    # lj.params[(particle_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
    # lj.r_cut[particle_types, 'Dummy'] = 0
    # lj.mode = 'shift'

    fene = hoomd.md.bond.FENEWCA()
    fene.params[sim.state.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
    fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params[sim.state.bond_types] = dict(k=30,r0=0.97)
    harmonic.params['Dummy'] = dict(k=0,r0=0.97)

    # if len(sim.state.angle_types)==0:
    #     FJ_system=True
    # else:
    #     FJ_system=False

    # if FJ_system==False:
    #     angle_constant = job.sp['angle_constant']
    #     cosinesq = hoomd.md.angle.CosineSquared()
    #     cosinesq.params[sim.state.angle_types] = dict(k=angle_constant, t0=np.pi*110/180)# https://www.sciencedirect.com/science/article/pii/S0032386110003642?ref=cra_js_challenge&fr=RR-1
    #     cosinesq.params['Dummy'] = dict(k=0.0001, t0=np.pi)  # k>0 to make warning go away (should not do anything)
    #-----------------------Force fields------------------------
    #----------------------Integration parameters-------------------------

    kT = 1 # Could read, but just assuming job.sp["temperature"]
    print(f"Could read, but assuming temperature = {kT}")
    dt = 0.005 # Also taken from Brandon's deformation parameter scan
    types_to_integrate =  hoomd.filter.Type(sim.state.particle_types[:-1]) # everything but "Dummy" particles

    integrator = hoomd.md.Integrator(dt=dt)
    
    integrator.forces = [harmonic]

    mttk = hoomd.md.methods.thermostats.MTTK(kT=kT,tau=100*dt)
    nvt = hoomd.md.methods.ConstantVolume(filter=types_to_integrate,thermostat=mttk)
    integrator.methods = [nvt]

    sim.operations.integrator = integrator
    #------------------------Configure output---------------------------
    logger = hoomd.logging.Logger(categories=['scalar','string'])
    logger.add(sim, quantities=['timestep','tps'])

    thermo = hoomd.md.compute.ThermodynamicQuantities(types_to_integrate)
    sim.operations.computes.append(thermo)
    logger.add(thermo, ['kinetic_temperature','kinetic_energy', 'potential_energy','pressure','volume'])
    
    box_writer_action = BoxWriter()
    logger.add(box_writer_action, ['Lx','Ly','Lz'])
    box_writer = hoomd.write.CustomWriter(action=box_writer_action,
                                             trigger=hoomd.trigger.Periodic(10))
    sim.operations.writers.append(box_writer)

    log_writer = hoomd.write.Table(output=open(parent_path + '/equi_feneOnly.log','w'),
                                    trigger=hoomd.trigger.Periodic(100),
                                    logger=logger,)
    sim.operations.writers.append(log_writer)

    logger = hoomd.logging.Logger(categories=['sequence'])
    logger.add(thermo, ['pressure_tensor'])
    gsd_writer = hoomd.write.GSD(filename=parent_path + '/equi_feneOnly.gsd',
                                    trigger=hoomd.trigger.Periodic(99999),
                                    dynamic=['property','momentum','topology','attribute'],
                                    mode='wb',
                                    logger=logger)
    sim.operations.writers.append(gsd_writer)

    #--------------------------Run the simulation----------------------

    sim.state.thermalize_particle_momenta(kT=kT,filter=types_to_integrate)

    # Run a short equilibration
    sim.run(1e5,write_at_start=True)

    print("Equilibration done! Initial state ready for deformation.")
    del sim
'''

def isolate_gel(job,input_gsd,frame_number):
    """
    Isolate the gel from the input gsd file and save it to a new gsd file.
    """
    import networkx as nx
    parent_path = str(Path(input_gsd).parent.absolute())
    output_gsd = parent_path + '/isolated_gel.gsd'

    traj = gsd.hoomd.open(input_gsd, mode='r')
    in_frame = traj[int(frame_number)]
    bonds = in_frame.bonds.group

    G = nx.Graph()
    G.add_edges_from(bonds)
    # find all disconnected/connected sub-networks
    # sort the list of networks by length
    Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
    # get the largest connected component
    gel_graph = Gcc[0]
    gel_nodes = list(gel_graph)
    print(f"Isolated {len(gel_nodes)} gel particles from {len(in_frame.particles.position)} total particles.")
    
    return gel_nodes


def main(job,gsd_file_path,frame_number):
    WCA_deformation = False # If False, will use LJ potential with r_cut=2.5
    feneOnly_deformation = True # If True, will not use pair potential, only FENEWCA
    if feneOnly_deformation:
        print("Using no pair deformation")
        WCA_deformation = False
    else:
        print("Using pair deformation")
        if WCA_deformation:
            print("Using WCA deformation")
        else:
            print("Using LJ deformation")

    parent_path = str(Path(gsd_file_path).parent.absolute())

    # Isolate the gel from the input gsd file
    gel_nodes = isolate_gel(job,gsd_file_path,frame_number)

    try:
        device = hoomd.device.GPU()
    except:
        device = hoomd.device.CPU()

    if feneOnly_deformation:
        print("No pair deformation, using FENEWCA only")
        # run_feneOnly_equilibration(job,gsd_file_path, frame_number)
        # sim = hoomd.Simulation(device=device, seed=1)
        # sim.create_state_from_gsd(filename=parent_path + '/equi_feneOnly.gsd', frame=-1)
        sim = hoomd.Simulation(device=device, seed=1)
        sim.create_state_from_gsd(filename=gsd_file_path, frame=int(frame_number))
    else:
        if WCA_deformation:
            run_equilibration(job,gsd_file_path, frame_number)
            sim = hoomd.Simulation(device=device, seed=1)
            sim.create_state_from_gsd(filename=parent_path + '/equi_WCA.gsd', frame=-1)
        else:
            print("Using existing gsd file for deformation")
            # If not WCA deformation, just read the gsd file directly
            sim = hoomd.Simulation(device=device, seed=1)
            sim.create_state_from_gsd(filename=gsd_file_path, frame=int(frame_number))

    # ------------------- Box resizer -------------------
    # Parameters taken from Brandon's previous parameter scans for NPT deformation of 
    # triblock system, except step size scaled down by a fact or 10. 1000 -> 100
    step_size = 1000
    delta_lam = 0.001
    # NVT deformation
    lam = 3.0
    deform_time = (lam - 1) / delta_lam * step_size
    traj = gsd.hoomd.open(gsd_file_path,mode='r')
    initial_box = traj[-1].configuration.box
    initial_box = np.asarray(initial_box)
    # initial_box = np.asarray(sim.state.box)

    deform_box_variant = NVTCustomBoxVariant(initial_box=initial_box,
                                              t_start=sim.timestep+100000,
                                              t_ramp=deform_time,
                                              final_lam=lam)

    box_resize = hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(step_size),
                                        box=deform_box_variant,)
    sim.operations.updaters.append(box_resize)

    #-----------------------Force fields------------------------
    particle_types = sim.state.particle_types
    cell = hoomd.md.nlist.Cell(buffer=0.4)

    lj = hoomd.md.pair.LJ(nlist=cell)
    lj.params[(particle_types, particle_types)] = dict(epsilon=1.0,sigma=1.0)
    if WCA_deformation:
        lj.r_cut[particle_types, particle_types] = 2**(1/6)
    else:
        lj.r_cut[particle_types, particle_types] = 2.5

    lj.params[(particle_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
    lj.r_cut[particle_types, 'Dummy'] = 0
    lj.mode = 'shift'

    fene = hoomd.md.bond.FENEWCA()
    fene.params[sim.state.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
    fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

    if len(sim.state.angle_types)==0:
        FJ_system=True
    else:
        FJ_system=False

    if FJ_system==False:
        angle_constant = job.sp['angle_constant']
        cosinesq = hoomd.md.angle.CosineSquared()
        cosinesq.params[sim.state.angle_types] = dict(k=angle_constant, t0=np.pi*110/180)# https://www.sciencedirect.com/science/article/pii/S0032386110003642?ref=cra_js_challenge&fr=RR-1
        cosinesq.params['Dummy'] = dict(k=0.0001, t0=np.pi)  # k>0 to make warning go away (should not do anything)


    #----------------------Integration parameters-------------------------

    # kT = 1 # Could read, but just assuming job.sp["temperature"]
    # print(f"Could read, but assuming temperature = {kT}")
    kT = job.sp['temperature']
    dt = 0.005 # Also taken from Brandon's deformation parameter scan
    # types_to_integrate =  hoomd.filter.Type(sim.state.particle_types[:-1]) # everything but "Dummy" particles
    types_to_integrate = hoomd.filter.Tags(gel_nodes) # filter for the gel particles

    integrator = hoomd.md.Integrator(dt=dt)
    
    if feneOnly_deformation:
        integrator.forces = [fene]
        if FJ_system==False:
            integrator.forces.append(cosinesq)
    else:
        integrator.forces = [lj,fene]
        if FJ_system==False:
            integrator.forces.append(cosinesq)

    mttk = hoomd.md.methods.thermostats.MTTK(kT=kT,tau=100*dt)
    nvt = hoomd.md.methods.ConstantVolume(filter=types_to_integrate,thermostat=mttk)
    integrator.methods = [nvt]

    sim.operations.integrator = integrator

    #------------------------Configure output---------------------------

    logger = hoomd.logging.Logger(categories=['scalar','string'])
    logger.add(sim, quantities=['timestep','tps'])

    thermo = hoomd.md.compute.ThermodynamicQuantities(types_to_integrate)
    sim.operations.computes.append(thermo)
    logger.add(thermo, ['kinetic_temperature','kinetic_energy', 'potential_energy','pressure','volume'])
    
    box_writer_action = BoxWriter()
    logger.add(box_writer_action, ['Lx','Ly','Lz'])
    box_writer = hoomd.write.CustomWriter(action=box_writer_action,
                                             trigger=hoomd.trigger.Periodic(10))
    sim.operations.writers.append(box_writer)

    #--------------------------Run the simulation----------------------

    sim.state.thermalize_particle_momenta(kT=kT,filter=types_to_integrate)

    sim.run(100000,write_at_start=True)
    print("Thermalization done! Starting deformation...")

    log_writer = hoomd.write.Table(output=open(parent_path + '/deform_gel.log','w'),
                                    trigger=hoomd.trigger.Periodic(10),
                                    logger=logger,)
    sim.operations.writers.append(log_writer)

    # add logging for pressure tensor and gsd file output
    logger = hoomd.logging.Logger(categories=['sequence'])
    logger.add(thermo, ['pressure_tensor'])
    gsd_writer = hoomd.write.GSD(filename=parent_path + '/deform_gel.gsd',
                                    # trigger=hoomd.trigger.Periodic(int(5)),
                                    trigger=hoomd.trigger.Periodic(int(deform_time/4/10)),
                                    dynamic=['property','momentum','topology','attribute'],
                                    mode='wb',
                                    logger=logger)
    sim.operations.writers.append(gsd_writer)

    #add pressure to a new gsd writer
    logger_pressure = hoomd.logging.Logger(categories = ['scalar','sequence'])
    logger_pressure.add(thermo, ['pressure_tensor'])
    gsd_writer_pressure = hoomd.write.GSD(filename=parent_path+"/pressures_gel.gsd",
                                          filter = hoomd.filter.Null(),
                                          mode = "wb", 
                                          trigger=hoomd.trigger.Periodic(10),
                                          logger=logger_pressure)
    sim.operations.writers.append(gsd_writer_pressure)

    sim.run(int(deform_time/4),write_at_start=True)

    # sim.run(deform_time+1)


    print("Deformation done! Written gsd file to: ",parent_path + "/deform_gel.gsd")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 deform_gel.py gsd_file_path frame_number")
        sys.exit(1)
    main(str(sys.argv[1]), str(sys.argv[2]))