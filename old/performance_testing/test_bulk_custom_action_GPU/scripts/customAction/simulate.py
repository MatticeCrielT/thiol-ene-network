from scripts.customAction.system import System
import numpy as np
import gsd,gsd.hoomd
import sys, os 
sys.path.insert(0,'/Users/statt/Programs/hoomd_fork_4.0.1/')
import hoomd
import scripts.customAction.polymerize as polymerize
from timeit import default_timer as timer

import datetime
class Status:
    def __init__(self, simulation):
        self.simulation = simulation

    @property
    def seconds_remaining(self):
        try:
            return (
                self.simulation.final_timestep - self.simulation.timestep
            ) / self.simulation.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))

class Simulator():
    def __init__(self,job):
        
        self.job = job
        self.kT = job.sp['temperature']
        self.rho = job.sp['density']
        self.crosslinker_percent = job.sp['crosslinker_percent']
        self.N_monomers = job.sp['N_monomers']
        self.monomer_size = job.sp['monomer_size']
        self.extender_size = job.sp['extender_size']
        self.angle_constant = job.sp['angle_constant']

        self.radical_number_percent = job.sp['radical_percent']
        self.r_cut_reaction = job.sp['r_cut_reaction']
        
        self.chain_transfer_probability= job.sp['chain_transfer_probability']
        self.chain_side_reaction_probability= job.sp['chain_side_reaction_probability']

        self.polymerize_period = job.sp['polymerize_period'] 

        # file names 
        self.equi_gsd_file = job.fn('equi.gsd')
        self.polymerize_gsd_file = job.fn('polymerize.gsd')
        self.output_txt =job.fn('polymerize.txt')
        self.output_run_txt =job.fn('run.txt')
        self.run_gsd_file = job.fn('run.gsd')

    def equilibrate(self): 
        S = System()
        frame = S.create_initial_configuration(density=self.rho,
                                            N_monomers = self.N_monomers,
                                            monomer_size=self.monomer_size,
                                            extender_size=self.extender_size,
                                            crosslinker=self.crosslinker_percent)
        

        with gsd.hoomd.open(name=self.equi_gsd_file, mode='w') as f:
            
            f.append(frame)

        # try:
        #     cpu = hoomd.device.GPU()
        # except:
        #     cpu = hoomd.device.CPU()
        if self.job.sp["polymerization_method"] == "custom_action_GPU" or self.job.sp["polymerization_method"] == "custom_action_GPU_bulk":
            device = hoomd.device.GPU(notice_level=3)
        else:
            device = hoomd.device.CPU(notice_level=3)
        print(self.job.sp["polymerization_method"]," is being run on ",device)

        sim = hoomd.Simulation(device=device, seed=1)
        sim.create_state_from_gsd(filename=self.equi_gsd_file)

        if len(sim.state.angle_types)==0:
            simple_system=True
        else:
            simple_system=False

        integrator = hoomd.md.Integrator(dt=0.005)
        cell = hoomd.md.nlist.Cell(buffer=0.4)

        lj = hoomd.md.pair.LJ(nlist=cell)

        lj.params[(S.particles_types, S.particles_types)] = dict(epsilon=1.0,sigma=1.0)
        lj.r_cut[S.particles_types, S.particles_types] = 2.5

        lj.params[(S.particles_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
        lj.r_cut[S.particles_types, 'Dummy'] = 0
        lj.mode = 'shift'

        fene = hoomd.md.bond.FENEWCA()
        fene.params[S.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
        fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

        gaussian = hoomd.md.pair.DPDConservative(nlist=cell)
        gaussian.params[(S.particles_types,S.particles_types)] = dict(A=50.0)
        gaussian.r_cut[S.particles_types,S.particles_types] = 1.0

        gaussian.params[(S.particles_types, 'Dummy')] = dict(A=0.0)
        gaussian.r_cut[S.particles_types, 'Dummy'] = 0
        gaussian.mode = 'none'

        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params[S.bond_types] = dict(k=100.0, r0=0.96)
        harmonic.params['Dummy'] = dict(k=0.0, r0=0.0)
        
        if simple_system==False:
            harmonic_a = hoomd.md.angle.Harmonic()
            harmonic_a.params[S.angle_types] = dict(k=self.angle_constant,    t0=np.pi)
            harmonic_a.params['Dummy'] = dict(k=0.0001, t0=np.pi)  # k>0 to make warning go away (should not do anything)

        types_to_integrate =  hoomd.filter.Type(S.particles_types[:-1]) # everything but "Dummy" particles

        print("Simulation set up, starting energy minimization...")

        nve = hoomd.md.methods.ConstantVolume(filter=types_to_integrate)
        fire = hoomd.md.minimize.FIRE(dt=0.05,
                                    force_tol=1e-2,
                                    angmom_tol=1e-2,
                                    energy_tol=1e-4)


        fire.methods.append(nve)
        fire.forces = [gaussian, harmonic]
        sim.operations.integrator = fire

        sim.run(1000)

        fire.forces = [gaussian,fene]
        sim.run(1000)

        fire.forces = [lj,fene]

        sim.run(1000)
        if simple_system==False:
            fire.forces = [lj,fene,harmonic_a]

            sim.run(1000)

        fire.forces = []
        del fire

        print("Simulation initialized, starting equilibration...")

        # Nose-Hoover
        #nvt = hoomd.md.methods.ConstantVolume(filter=types_to_integrate,
        #    thermostat=hoomd.md.methods.thermostats.MTTK(kT=kT,tau=0.005*100))
        #integrator.methods.append(nvt)

        langevin = hoomd.md.methods.Langevin(filter=types_to_integrate, kT=self.kT,default_gamma=0.1)

        integrator.methods.append(langevin)
        
        if simple_system==False:
            integrator.forces = [lj,fene,harmonic_a]
        else:
            integrator.forces = [lj,fene]

        sim.operations.integrator = integrator
        sim.state.thermalize_particle_momenta(filter=types_to_integrate, kT=self.kT)

        sim.run(5000)
        hoomd.write.GSD.write(state=sim.state, mode='wb', filename=self.equi_gsd_file)
        print("writing gsd file to:",self.equi_gsd_file)

    def polymerize(self):
        
        S = System()

        # try:
        #     cpu = hoomd.device.GPU(notice_level=5)
        # except:
        #     cpu = hoomd.device.CPU(notice_level=5)
        if self.job.sp["polymerization_method"] == "custom_action_GPU" or self.job.sp["polymerization_method"] == "custom_action_GPU_bulk":
            device = hoomd.device.GPU(notice_level=3)
        else:
            device = hoomd.device.CPU(notice_level=3)
        print(f"job_id: {self.job.id}",self.job.sp["polymerization_method"]," is being run on ",device)

        sim = hoomd.Simulation(device=device, seed=1)

        if os.path.isfile(self.polymerize_gsd_file) and os.path.getsize(self.polymerize_gsd_file)>1e4:
            sim.create_state_from_gsd(filename=self.polymerize_gsd_file)
        else: 
            sim.create_state_from_gsd(filename=self.equi_gsd_file)
        
        if len(sim.state.angle_types)==0:
            simple_system=True
        else:
            simple_system=False

        integrator = hoomd.md.Integrator(dt=0.005)
        sim.operations.integrator = integrator 

        cell = hoomd.md.nlist.Cell(buffer=0.4)

        lj = hoomd.md.pair.LJ(nlist=cell)
        lj.params[(S.particles_types, S.particles_types)] = dict(epsilon=1.0,sigma=1.0)
        lj.r_cut[S.particles_types, S.particles_types] = 2.5

        lj.params[(S.particles_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
        lj.r_cut[S.particles_types, 'Dummy'] = 0
        lj.mode = 'shift'

        fene = hoomd.md.bond.FENEWCA()
        fene.params[S.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
        fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

        if simple_system==False:
            harmonic_a = hoomd.md.angle.Harmonic()
            harmonic_a.params[S.angle_types] = dict(k=self.angle_constant,    t0=np.pi)
            harmonic_a.params['Dummy'] = dict(k=0.0001, t0=np.pi)  # k>0 to make warning go away (should not do anything)

            integrator.forces = [lj,fene,harmonic_a]
        else:
            integrator.forces = [lj,fene]

        types_to_integrate =  hoomd.filter.Type(S.particles_types[:-1])

       

        if os.path.isfile(self.polymerize_gsd_file) and os.path.getsize(self.polymerize_gsd_file)>1e4:
            # already has radicals in polymerize.gsd 
            print("Continue reaction...")
        else:
            print("Starting reaction...")
            # flip some sulfurs to be radical, i.e. "initialization by photoinitiator"
            snapshot = sim.state.get_snapshot()
            S.flip_radicals_on(snapshot,self.radical_number_percent)
            sim.state.set_snapshot(snapshot)
        
        npt = hoomd.md.methods.ConstantPressure(
            filter=types_to_integrate,
            tauS=1000*sim.operations.integrator.dt,
            gamma=2/(1000*sim.operations.integrator.dt),
            S=0.0,
            couple="xyz",
            rescale_all=True,
            thermostat=hoomd.md.methods.thermostats.MTTK(kT=self.kT,tau=sim.operations.integrator.dt*100))
        
        sim.state.thermalize_particle_momenta(filter=types_to_integrate, kT=self.kT)
        #zero_momentum = hoomd.md.update.ZeroMomentum( hoomd.trigger.On(1000))
        #sim.operations.updaters.append(zero_momentum)
        
        sim.operations.integrator.methods.append(npt)

        # Define and add the GSD operation.
        gsd_writer = hoomd.write.GSD(filename=self.polymerize_gsd_file,
                                    trigger=hoomd.trigger.Periodic(100*self.polymerize_period),
                                    dynamic=['property','momentum','topology','attribute'],
                                    mode='ab')
        sim.operations.writers.append(gsd_writer)
        print("writing gsd file to:",self.polymerize_gsd_file)
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=types_to_integrate)
        sim.operations.computes.append(thermodynamic_properties)
       
        status = Status(sim)
        logger = hoomd.logging.Logger(categories=['scalar','string'])
       
        logger.add(sim, quantities=['timestep','tps'])
        logger.add(thermodynamic_properties, quantities=['kinetic_temperature','pressure','kinetic_energy','potential_energy','volume'])

        table_stdout = hoomd.write.Table(trigger=hoomd.trigger.Periodic(self.polymerize_period,100),logger=logger)
        table_file = hoomd.write.Table(trigger=hoomd.trigger.Periodic(self.polymerize_period,100),logger=logger, output=open(self.output_txt,'a'))
        sim.operations.writers.append(table_stdout)
        sim.operations.writers.append(table_file)

        if self.job.sp["polymerization_method"] == "custom_action_GPU" or self.job.sp["polymerization_method"] == "custom_action_CPU" :
            #define and add the actions and operations to keep track of extent of reaction
            extent_of_reaction_action = polymerize.calc_extent_of_reaction(
                                        simulator=self) 
            extent_of_reaction_operation = hoomd.update.CustomUpdater(
                action=extent_of_reaction_action, trigger=self.polymerize_period
            )
            sim.operations += extent_of_reaction_operation

            #define and add the polymerization actions and operations
            bond_formation_action = polymerize.bond_formation(
                                        probability=1.0,r_cut=self.r_cut_reaction,
                                        simple_system=simple_system) 
                #thiol reaction probability defaults to 1.0 in system's propagate reactions
            bond_formation_operation = hoomd.update.CustomUpdater(
                action=bond_formation_action, trigger=self.polymerize_period
            )
            sim.operations += bond_formation_operation

            # polymerize.chain_growth(self.)
            chain_growth_action = polymerize.chain_growth(
                                        probability=self.chain_side_reaction_probability,
                                        r_cut=self.r_cut_reaction,simple_system=simple_system)
            chain_growth_operation = hoomd.update.CustomUpdater(
                action=chain_growth_action, trigger=self.polymerize_period
            )
            sim.operations += chain_growth_operation

            # polymerize.chain_transfer(self.)
            chain_transfer_action = polymerize.chain_transfer(
                                            probability=self.chain_transfer_probability,
                                            r_cut=self.r_cut_reaction,simple_system=simple_system) 
            chain_transfer_operation = hoomd.update.CustomUpdater(
                action=chain_transfer_action, trigger=self.polymerize_period
            )
            sim.operations += chain_transfer_operation

            '''
            termination_action = polymerize.termination_reactions(
                                    probability=1.0,r_cut=self.r_cut_reaction,simple_system=simple_system)
                #thiol reaction probability defaults to 1.0 in system's propagate reactions
            termination_operation = hoomd.update.CustomUpdater(
                action=termination_action, trigger=self.polymerize_period
            )
            sim.operations += termination_operation
            '''
        elif self.job.sp["polymerization_method"] == "custom_action_GPU_bulk":
            bulk_polymerization_action = polymerize.bulk_polymerize(
                                    r_cut=self.r_cut_reaction,simple_system=simple_system,simulator=self)
            bulk_polymerization_operation = hoomd.update.CustomUpdater(
                action=bulk_polymerization_action, trigger=self.polymerize_period
            )
            sim.operations += bulk_polymerization_operation


        polymerization_times = []
        integration_times = []
        integration_tps = []
        tps = []

        simulation_start = timer()
        
        for i in range(10000):
            period_start = timer()

            polymerization_start = timer()
            sim.run(1)
            polymerization_end = timer()
            polymerization_times.append(polymerization_end - polymerization_start)

            integration_start = timer()
            sim.run(self.polymerize_period-1)
            integration_end = timer()
            integration_times.append(integration_end - integration_start)

            integration_tps.append((self.polymerize_period-1)/integration_times[-1])

            gsd_writer.flush()

            period_end = timer()
            tps.append(period_end - period_start)

            if (period_end - simulation_start > 30*60):
                self.job.doc["TPS"] = np.average(tps)
                self.job.doc["avg_polymerization_time"] = np.average(polymerization_times)
                self.job.doc["avg_integration_time"] = np.average(integration_times)
                self.job.doc["integration_tps"] = np.average(integration_tps)
                exit()

    # def run(self):

    #     S = System()
        
    #     try:
    #         cpu = hoomd.device.GPU()
    #     except:
    #         cpu = hoomd.device.CPU()
        
    #     sim = hoomd.Simulation(device=cpu, seed=1)
    #     sim.create_state_from_gsd(filename=self.polymerize_gsd_file)

    #     if len(sim.state.angle_types)==0:
    #         simple_system=True
    #     else:
    #         simple_system=False

    #     integrator = hoomd.md.Integrator(dt=0.005)
    #     sim.operations.integrator = integrator

    #     cell = hoomd.md.nlist.Cell(buffer=0.4)

    #     lj = hoomd.md.pair.LJ(nlist=cell)
    #     lj.params[(S.particles_types, S.particles_types)] = dict(epsilon=1.0,sigma=1.0)
    #     lj.r_cut[S.particles_types, S.particles_types] = 2.5

    #     lj.params[(S.particles_types, 'Dummy')] = dict(epsilon=0.0,sigma=0.0)
    #     lj.r_cut[S.particles_types, 'Dummy'] = 0
    #     lj.mode = 'shift'

    #     fene = hoomd.md.bond.FENEWCA()
    #     fene.params[S.bond_types] = dict(k=30,r0=1.5,epsilon=1.0, sigma=1.0, delta=0.0)
    #     fene.params['Dummy'] = dict(k=0,r0=1.5,epsilon=0.0, sigma=1.0, delta=0.0)

    #     if simple_system==False:
    #         harmonic_a = hoomd.md.angle.Harmonic()
    #         harmonic_a.params[S.angle_types] = dict(k=self.angle_constant,    t0=np.pi)
    #         harmonic_a.params['Dummy'] = dict(k=1e-6, t0=np.pi)  # k>0 to make warning go away (should not do anything)

    #         integrator.forces = [lj,fene,harmonic_a]
    #     else: 
    #         integrator.forces = [lj,fene]

    #     types_to_integrate =  hoomd.filter.Type(S.particles_types[:-1])

    #     npt = hoomd.md.methods.ConstantPressure(
    #         filter=types_to_integrate,
    #         tauS=1000*sim.operations.integrator.dt,
    #         gamma=2/(1000*sim.operations.integrator.dt),
    #         S=0.0,
    #         couple="xyz",
    #         rescale_all=True,
    #         thermostat=hoomd.md.methods.thermostats.MTTK(kT=self.kT,tau=sim.operations.integrator.dt*100))
        
    #     sim.state.thermalize_particle_momenta(filter=types_to_integrate, kT=self.kT)
    #     #zero_momentum = hoomd.md.update.ZeroMomentum( hoomd.trigger.On(1000))
    #     #sim.operations.updaters.append(zero_momentum)
        
    #     sim.operations.integrator.methods.append(npt)

    #     # Define and add the GSD operation.
    #     gsd_writer = hoomd.write.GSD(filename=self.run_gsd_file,
    #                                 trigger=hoomd.trigger.Periodic(1000),
    #                                 dynamic=['property','momentum','topology','attribute'],
    #                                 mode='wb')
    #     sim.operations.writers.append(gsd_writer)

    #     thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=types_to_integrate)
    #     sim.operations.computes.append(thermodynamic_properties)
       
    #     logger = hoomd.logging.Logger(categories=['scalar'])
       
    #     logger.add(sim, quantities=['timestep'])
    #     logger.add(thermodynamic_properties, quantities=['kinetic_temperature','pressure','kinetic_energy','potential_energy','volume'])

    #     table_stdout = hoomd.write.Table(trigger=hoomd.trigger.Periodic(1000),logger=logger)
    #     table_file = hoomd.write.Table(trigger=hoomd.trigger.Periodic(1000),logger=logger, output=open(self.output_run_txt,'w'))
    #     sim.operations.writers.append(table_stdout)
    #     sim.operations.writers.append(table_file)

    #     sim.run(10000000)



