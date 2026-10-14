"""
This file contains classes defined as custom actions for use in the generation of
thiol-ene networks. To use these custom actions, include the path to this file, and
then call 
`import polymerize`
These custom actions inherit from Hoomd's custom actions, which were introduced in
Hoomdv4
"""

import hoomd
import numpy as np
import freud

from timeit import default_timer as timer

class calc_extent_of_reaction(hoomd.custom.Action):
    '''
    This class calculates the extent of reaction of the system, using the number of
        unreacted ene-group particles divided by the number of potential bonds, equal to
        the number of initial monomers
    '''
    def __init__(self,simulator):
        self.simulator = simulator

    def attach(self, simulation):
        super().attach(simulation)

    def act(self,timestep):
        snapshot = self._state.get_snapshot()

        # dummy_id = len(sim.state.bond_types)-1
        # bonds = snapshot.bonds.group[snapshot.bonds.typeid!=dummy_id] # remove dummy type bond 
        # nids, counts = np.unique(np.concatenate(bonds).flatten(),return_counts=True)
        # bonded = nids[counts>=2]
        # ids = np.arange(len(snapshot.particles.tag))
        # idx = snapshot.particles.rtag[ids]
        ids = np.arange(len(snapshot.particles.charge))
        idx = np.arange(len(snapshot.particles.charge))
        particle_ids = snapshot.particles.typeid[idx]

        unreacted_enes = len(ids[particle_ids==1])/2

        #reacted_monomers = bonded[snapshot.bonds.typeid[bonded]==1]
        # self.job.doc['reacted_monomers'] = 1 - unreacted_enes/(self.N_monomers*2)
        self.simulator.job.doc['reacted_monomers'] = 1 - unreacted_enes/(self.simulator.N_monomers*2)



class bond_formation(hoomd.custom.Action):
    '''
    thiol-ene propogation reaction: thiol radicals add terminally to ene group
    thiol radical (type 2) adds to terminal ene (type 1) that has only one bond
    >-<* + o=o-o=o = >-<o*=o-o=o
    0-2  + 1=1-1=1 = 0-0-3-1-1=1
    '''
    def __init__(self,probability,r_cut,FJ_system):
        ## define chemical tags ----
        self.particles_types = ['Thiol','Ene_C','RSulfur','RCarbon', 'Sulfur', 'Carbon','C','D','Dummy']
        ## ---------------[  0  , 1  ,     2     ,  3 ,     4     ,   5]
        self.bond_types = ['-SH','New','CC_Double','-S-','CC_Single','Dummy']
        self.angle_types = ['PE','Dummy']

        # needs to be the correct integer position of 'dummy' in the arrays above 
        self.dummy_type = len(self.particles_types)-1
        self.dummy_type_bond = len(self.bond_types)-1
        self.dummy_type_angle = len(self.angle_types)-1

        self.thiol_type = 0
        self.ene_type = 1
        self.radical_thiol = 2
        self.radical_carbon = 3
        self.sulfur = 4
        self.carbon = 5
        self.spacer_type1 = 6
        self.spacer_type2 = 7

        self.thiol_reaction_probability = probability
        self.rcut_neigh = r_cut

        self.FJ_system = FJ_system
        
    def attach(self, simulation):
        super().attach(simulation)

    def act(self, timestep):
        start = timer()

        snapshot = self._state.get_snapshot()

        end = timer()
        print("bond_propagation_get_snapshot: ", end - start)

        bond_propagation_start = timer()

        if self.thiol_reaction_probability>0:
            # Read the simulation state
            '''
            ids = np.arange(len(snapshot.particles.tag))
            idx = snapshot.particles.rtag[ids]
            '''
            ids = np.arange(len(snapshot.particles.charge))
            idx = np.arange(len(snapshot.particles.charge))

            particle_ids = snapshot.particles.typeid[idx]
            positions = snapshot.particles.position[idx]
            '''
            box = freud.Box.from_box(snapshot.global_box)
            '''
            box = freud.Box.from_box(snapshot.configuration.box)

            # Propagation/bond formation step 
            # find the positions of all thiol radicals
            radicals_positions = positions[particle_ids==self.radical_thiol]
            
            # find the ids of thiol-radicals and enes
            radicals_ids = ids[particle_ids==self.radical_thiol]
            polymer_ids = ids[particle_ids==self.ene_type]  # all Enes

            # bond type = 2 carbon double bond. 
            # bond type = 4 carbon single bond - can't react in bond formation 

            # count and find all unreacted c=c
            bond_ids,bond_counts = np.unique((snapshot.bonds.group[snapshot.bonds.typeid==2]).flatten(),return_counts=True)
            # find all particles with only 1 non-H covalent bond
            only_one_bond_ids = bond_ids[bond_counts==1]
            
            # the ene groups that are candidates for rxn only have one bond
            polymer_candidates_ids = np.intersect1d(only_one_bond_ids,polymer_ids)
            polymer_candidates_positions = positions[polymer_candidates_ids]

            # only if there are radicals and reactable ene groups
            if len(radicals_ids)>0 and len(polymer_candidates_ids)>0:
                # get the radical positions
                aq = freud.locality.AABBQuery(box, radicals_positions)
                # find the positions of polymer candidates near radical
                nlist = aq.query(polymer_candidates_positions, {'r_max': self.rcut_neigh, 'exclude_ii':True}).toNeighborList()
                all_pairs = nlist[:]

                # get the radical and polymer ids 
                actual_ids_polymers = polymer_candidates_ids[all_pairs[:,0]]
                actual_ids_radicals = radicals_ids[all_pairs[:,1]]
                actual_ids = np.vstack((actual_ids_radicals,actual_ids_polymers)).T

                # for each thiol radical a, there is a chance to react
                for a in np.unique(actual_ids_radicals):
                    # find the neighbors of the thiol radical a that are radical thiols or polymer candidates
                    neigh_a = actual_ids[actual_ids[:,0]==a][:,1]
                    neigh_a = neigh_a.astype(int)

                    # get the types of the neighbors
                    neigh_a_types = particle_ids[neigh_a]
                    # take the neighbors that are not radical thiols or radical carbons
                    neigh_a = neigh_a[(neigh_a_types!=self.radical_thiol)|(neigh_a_types!=self.radical_carbon)]
                    bond_ids,bond_counts = np.unique((snapshot.bonds.group).flatten(),return_counts=True)
                    # the neighbors also must have only 1 bond (terminal carbon)
                    neigh_a = neigh_a[bond_counts[neigh_a]==1]
                    
                    # roll to react
                    react = np.random.uniform()

                    # if reacts
                    if len(neigh_a)>0 and react <= self.thiol_reaction_probability:  
                        # chose a random terminal ene neighbor 
                        b = np.random.choice(neigh_a)  

                        # find b and c it's bonded to b=c-|
                        bonds_on_b_1 = snapshot.bonds.group[snapshot.bonds.group[:,0]==b]
                        bonds_on_b_2 = snapshot.bonds.group[snapshot.bonds.group[:,1]==b]
                        bonds_on_b = np.unique(np.vstack((bonds_on_b_1,bonds_on_b_2)).flatten())
                        bonds_on_b = bonds_on_b[bonds_on_b!=b]

                        # find a and the atoms it's bonded to
                        bonds_on_a_1 = snapshot.bonds.group[snapshot.bonds.group[:,0]==a]
                        bonds_on_a_2 = snapshot.bonds.group[snapshot.bonds.group[:,1]==a]
                        bonds_on_a = np.unique(np.vstack((bonds_on_a_1,bonds_on_a_2)).flatten())
                        bonds_on_a = bonds_on_a[bonds_on_a!=a]
                        
                        # if not FJ_system then include angles
                        if self.FJ_system == False: 
                            U = np.where(snapshot.angles.typeid==self.dummy_type_angle)[0][0]
                            # good
                            snapshot.angles.group[U]=[bonds_on_b[0],b,a]  # angle 
                            snapshot.angles.typeid[U]=0
                            # good
                            snapshot.angles.group[U+1]=[bonds_on_a[0],a,b] #angle 
                            snapshot.angles.typeid[U+1]=0

                        # flip types of bonds 
                        snapshot.bonds.typeid[snapshot.bonds.group[:,0]==a]=3
                        snapshot.bonds.typeid[snapshot.bonds.group[:,1]==a]=3
                        # c-c single bond now
                        snapshot.bonds.typeid[snapshot.bonds.group[:,0]==b]=4
                        snapshot.bonds.typeid[snapshot.bonds.group[:,1]==b]=4

                        # make new bond 
                        U = np.where(snapshot.bonds.typeid==self.dummy_type_bond)[0][0]
                        snapshot.bonds.group[U]=[a,b]  # propagation 
                        snapshot.bonds.typeid[U]=1

                        snapshot.particles.typeid[idx[a]]=self.sulfur  # flip back to Thiol
                        snapshot.particles.typeid[idx[b]]=self.carbon   
                        snapshot.particles.typeid[idx[bonds_on_b[0]]]=self.radical_carbon   # turn the B into a radical 
                        
                        particle_ids[a]=self.sulfur
                        particle_ids[b]=self.carbon
                        particle_ids[bonds_on_b[0]]=self.radical_carbon
                        #flipping_rads.append(bonds_on_b[0])
                        actual_ids = np.vstack([i if i[1]!=b else [None,None] for i in actual_ids])
        
        start = timer()
        print("bond_propagation_logic: ",start - bond_propagation_start)

        self._state.set_snapshot(snapshot)

        end = timer()

        print("bond_propagation_set_snapshot: ", end-start)

class chain_transfer(hoomd.custom.Action):
    '''
    chain transfer of thiol-ene reaction: H from SH abstracts to carbon radical
    effectively radical jumps from carbon to sulfur
    ''' 
    def __init__(self,probability,r_cut,FJ_system):
        ## define chemical tags ----
        self.particles_types = ['Thiol','Ene_C','RSulfur','RCarbon', 'Sulfur', 'Carbon','C','D','Dummy']
        ## ---------------[  0  , 1  ,     2     ,  3 ,     4     ,   5]
        self.bond_types = ['-SH','New','CC_Double','-S-','CC_Single','Dummy']
        self.angle_types = ['PE','Dummy']

        # needs to be the correct integer position of 'dummy' in the arrays above 
        self.dummy_type = len(self.particles_types)-1
        self.dummy_type_bond = len(self.bond_types)-1
        self.dummy_type_angle = len(self.angle_types)-1

        self.thiol_type = 0
        self.ene_type = 1
        self.radical_thiol = 2
        self.radical_carbon = 3
        self.sulfur = 4
        self.carbon = 5
        self.spacer_type1 = 6
        self.spacer_type2 = 7

        self.chain_transfer_probability = probability
        self.rcut_neigh = r_cut

        self.FJ_system = FJ_system
        
    def attach(self, simulation):
        super().attach(simulation)
    
    def act(self, timestep):
        snapshot = self._state.get_snapshot()
        # transfer carbon radical to new thiol 
        if self.chain_transfer_probability>0:
            # Read the simulation state
            '''
            ids = np.arange(len(snapshot.particles.tag))
            idx = snapshot.particles.rtag[ids]
            '''
            ids = np.arange(len(snapshot.particles.charge))
            idx = np.arange(len(snapshot.particles.charge))

            particle_ids = snapshot.particles.typeid[idx]
            positions = snapshot.particles.position[idx]
            '''
            box = freud.Box.from_box(snapshot.global_box)
            '''
            box = freud.Box.from_box(snapshot.configuration.box)

            # get carbon radicals
            radicals_positions = positions[particle_ids==self.radical_carbon]
            radicals_ids = ids[particle_ids==self.radical_carbon]

            # find thiols that are candidates to move the radical to
            polymer_ids = ids[particle_ids==self.thiol_type]  # all Thiols
            bond_ids,bond_counts = np.unique((snapshot.bonds.group).flatten(),return_counts=True)
            only_one_bond_ids = bond_ids[bond_counts==1]
            polymer_candidates_ids = np.intersect1d(only_one_bond_ids,polymer_ids)
            polymer_candidates_positions = positions[polymer_candidates_ids]

            if len(radicals_ids)>0:
                # find the local candidates for each radical 
                aq = freud.locality.AABBQuery(box, radicals_positions)
                nlist = aq.query(polymer_candidates_positions, {'r_max': self.rcut_neigh, 'exclude_ii':True}).toNeighborList()
                all_pairs = nlist[:]

                actual_ids_polymers = polymer_candidates_ids[all_pairs[:,0]]
                actual_ids_radicals = radicals_ids[all_pairs[:,1]]
                actual_ids = np.vstack((actual_ids_radicals,actual_ids_polymers)).T  

                for a in np.unique(actual_ids_radicals):  # chain transfer 
                    
                    neigh_a = actual_ids[actual_ids[:,0]==a][:,1]
                    neigh_a = neigh_a.astype(int)
                    neigh_a_types = particle_ids[neigh_a]
                    
                    neigh_a = neigh_a[(neigh_a_types!=self.radical_thiol)|(neigh_a_types!=self.radical_carbon)]
                    bond_ids,bond_counts = np.unique((snapshot.bonds.group).flatten(),return_counts=True)
                    neigh_a = neigh_a[bond_counts[neigh_a]==1]
                    # roll to transfer
                    transfer = np.random.uniform()

                    # if it transfers
                    if len(neigh_a)>0 and transfer <= self.chain_transfer_probability:  
                        b = np.random.choice(neigh_a)
                        snapshot.particles.typeid[idx[a]]=self.carbon  # carbon is no longer reactive
                        snapshot.particles.typeid[idx[b]]=self.radical_thiol  # chain transfer 

                        # transfer the radical
                        particle_ids[a]=self.carbon
                        particle_ids[b]=self.radical_thiol

                        # remove the radical thiol as a possible candidate
                        actual_ids = np.vstack([i if i[1]!=b else [None,None] for i in actual_ids])
            
        self._state.set_snapshot(snapshot)

class chain_growth(hoomd.custom.Action):
    '''
    args: self, snapshot, box, positions, paricle_ids, ids, idx
    chain growth of ene carbons: a ene radical propogates from carbon to carbon double bond
    ene radical (type 2) adds to terminal ene (type 1) that has only one bond
    r-o-o*-r' + o=o-r'' -> r-o-o-r'
                               |
                               o-o*-r''
    r-5-4-r' + 2-2-r''  -> r-5-5-r'
                               |
                               5-4-r''
    '''
    def __init__(self,probability,r_cut,FJ_system):
        ## define chemical tags ----
        self.particles_types = ['Thiol','Ene_C','RSulfur','RCarbon', 'Sulfur', 'Carbon','C','D','Dummy']
        ## ---------------[  0  , 1  ,     2     ,  3 ,     4     ,   5]
        self.bond_types = ['-SH','New','CC_Double','-S-','CC_Single','Dummy']
        self.angle_types = ['PE','Dummy']

        # needs to be the correct integer position of 'dummy' in the arrays above 
        self.dummy_type = len(self.particles_types)-1
        self.dummy_type_bond = len(self.bond_types)-1
        self.dummy_type_angle = len(self.angle_types)-1

        self.thiol_type = 0
        self.ene_type = 1
        self.radical_thiol = 2
        self.radical_carbon = 3
        self.sulfur = 4
        self.carbon = 5
        self.spacer_type1 = 6
        self.spacer_type2 = 7

        self.chain_side_reaction_probability = probability
        self.rcut_neigh = r_cut

        self.FJ_system = FJ_system

    def attach(self, simulation):
        super().attach(simulation)

    def act(self, timestep):
        snapshot = self._state.get_snapshot()
        # chain side reaction radical carbon with ene reaction                
        if self.chain_side_reaction_probability>0:
            # Read the simulation state
            '''
            ids = np.arange(len(snapshot.particles.tag))
            idx = snapshot.particles.rtag[ids]
            '''
            ids = np.arange(len(snapshot.particles.charge))
            idx = np.arange(len(snapshot.particles.charge))

            particle_ids = snapshot.particles.typeid[idx]
            positions = snapshot.particles.position[idx]
            '''
            box = freud.Box.from_box(snapshot.global_box)
            '''
            box = freud.Box.from_box(snapshot.configuration.box)

            # find all radical carbons and get their ids
            radicals_positions = positions[particle_ids==self.radical_carbon]
            radicals_ids = ids[particle_ids==self.radical_carbon]
            '''
            if:
            remove any candidates that are on the same monomer
            '''

            # find reaction candidates
            ## find all ids of enes
            polymer_ids = ids[particle_ids==self.ene_type]   # all carbon/enes
            ## find the enes that have only one bond (which is a double bond)
            bond_ids,bond_counts = np.unique((snapshot.bonds.group[snapshot.bonds.typeid==2]).flatten(),return_counts=True)
            only_one_bond_ids = bond_ids[bond_counts==1]
            polymer_candidates_ids = np.intersect1d(only_one_bond_ids,polymer_ids)
            polymer_candidates_positions = positions[polymer_candidates_ids]

            # if there are carbon radicals to react
            if len(radicals_ids)>0:
                # find and store the local candidates around each radical
                aq = freud.locality.AABBQuery(box, radicals_positions)
                nlist = aq.query(polymer_candidates_positions, {'r_max': self.rcut_neigh, 'exclude_ii':True}).toNeighborList()
                all_pairs = nlist[:]

                actual_ids_polymers = polymer_candidates_ids[all_pairs[:,0]]
                actual_ids_radicals = radicals_ids[all_pairs[:,1]]
                actual_ids = np.vstack((actual_ids_radicals,actual_ids_polymers)).T  
                # for each carbon radical a
                for a in np.unique(actual_ids_radicals):  # radical carbon -ene bond formation 
                    # take the neighbor candidates around the radicals
                    neigh_a = actual_ids[actual_ids[:,0]==a][:,1]
                    neigh_a = neigh_a.astype(int)
                    neigh_a_types = particle_ids[neigh_a]
                    
                    neigh_a = neigh_a[(neigh_a_types!=self.radical_thiol)|(neigh_a_types!=self.radical_carbon)]
                    bond_ids,bond_counts = np.unique((snapshot.bonds.group).flatten(),return_counts=True)
                    neigh_a = neigh_a[bond_counts[neigh_a]==1]

                    # roll to react
                    react = np.random.uniform()
                    # if reacts
                    if len(neigh_a)>0 and react <= self.chain_side_reaction_probability: 
                        # choose one of the reactive candidates b
                        b = np.random.choice(neigh_a)  
                    
                        # find b and the c it's bonded to b=c-|
                        bonds_on_b_1 = snapshot.bonds.group[snapshot.bonds.group[:,0]==b]
                        bonds_on_b_2 = snapshot.bonds.group[snapshot.bonds.group[:,1]==b]
                        bonds_on_b = np.unique(np.vstack((bonds_on_b_1,bonds_on_b_2)).flatten())
                        bonds_on_b = bonds_on_b[bonds_on_b!=b]

                        # find a and the atoms it's bonded to 
                        ''' this search maybe unnecessary for this rxn?'''
                        bonds_on_a_1 = snapshot.bonds.group[snapshot.bonds.group[:,0]==a]
                        bonds_on_a_2 = snapshot.bonds.group[snapshot.bonds.group[:,1]==a]
                        bonds_on_a = np.unique(np.vstack((bonds_on_a_1,bonds_on_a_2)).flatten())
                        bonds_on_a = bonds_on_a[bonds_on_a!=a]
                        
                        # if not FJ system then include angles
                        if self.FJ_system == False: 
                            U = np.where(snapshot.angles.typeid==self.dummy_type_angle)[0][0]
                            
                            n_angles_added = 0
                            for next_neighbor in bonds_on_b:
                                snapshot.angles.group[U+n_angles_added]=[bonds_on_b[n_angles_added],b,a]  # angle 
                                snapshot.angles.typeid[U+n_angles_added]=0
                                n_angles_added += 1
                            for next_neighbor in bonds_on_a:
                                snapshot.angles.group[U+n_angles_added]=[bonds_on_a[n_angles_added-len(bonds_on_b)],a,b] #angle 
                                snapshot.angles.typeid[U+n_angles_added]=0
                                n_angles_added += 1
                        
                        # flip types of bonds 
                        ''' unnecessary?
                        snapshot.bonds.typeid[snapshot.bonds.group[:,0]==a]=4
                        snapshot.bonds.typeid[snapshot.bonds.group[:,1]==a]=4
                        '''
                        snapshot.bonds.typeid[snapshot.bonds.group[:,0]==b]=4
                        snapshot.bonds.typeid[snapshot.bonds.group[:,1]==b]=4


                        U = np.where(snapshot.bonds.typeid==self.dummy_type_bond)[0][0]
                        snapshot.bonds.group[U]=[a,b]  # propagation 
                        snapshot.bonds.typeid[U]=1

                        # flip to non-reactive carbons
                        snapshot.particles.typeid[idx[a]]=self.carbon 
                        snapshot.particles.typeid[idx[b]]=self.carbon   
                        snapshot.particles.typeid[idx[bonds_on_b[0]]]=self.radical_carbon   # turn the B into a radical 
                        
                        particle_ids[a]=self.carbon
                        particle_ids[b]=self.carbon
                        particle_ids[bonds_on_b[0]]=self.radical_carbon

                        # remove the reacted ones as candidate possibilities from the other radicals
                        actual_ids = np.vstack([i if i[1]!=b else [None,None] for i in actual_ids])

        self._state.set_snapshot(snapshot)

class termination_reactions(hoomd.custom.Action):
    '''
    termination reactions by which there is a net loss of radicals.
    Occurs when two radicals meet and form a new bond between them, and lose both radicals
    The types of radicals do not matter
    '''
    def __init__(self,probability,r_cut,FJ_system):
        ## define chemical tags ----
        self.particles_types = ['Thiol','Ene_C','RSulfur','RCarbon', 'Sulfur', 'Carbon','C','D','Dummy']
        ## ---------------[  0  , 1  ,     2     ,  3 ,     4     ,   5]
        self.bond_types = ['-SH','New','CC_Double','-S-','CC_Single','Dummy']
        self.angle_types = ['PE','Dummy']

        # needs to be the correct integer position of 'dummy' in the arrays above 
        self.dummy_type = len(self.particles_types)-1
        self.dummy_type_bond = len(self.bond_types)-1
        self.dummy_type_angle = len(self.angle_types)-1

        self.thiol_type = 0
        self.ene_type = 1
        self.radical_thiol = 2
        self.radical_carbon = 3
        self.sulfur = 4
        self.carbon = 5
        self.spacer_type1 = 6
        self.spacer_type2 = 7

        self.termination_probability = probability
        self.rcut_neigh = r_cut

        self.FJ_system = FJ_system

    def attach(self, simulation):
        super().attach(simulation)
    
    def act(self, timestep):
        snapshot = self._state.get_snapshot()

        # Read the simulation state
        '''
        ids = np.arange(len(snapshot.particles.tag))
        idx = snapshot.particles.rtag[ids]
        '''
        ids = np.arange(len(snapshot.particles.charge))
        idx = np.arange(len(snapshot.particles.charge))

        particle_ids = snapshot.particles.typeid[idx]
        positions = snapshot.particles.position[idx]
        '''
        box = freud.Box.from_box(snapshot.global_box)
        '''
        box = freud.Box.from_box(snapshot.configuration.box)

        # radical- radical reaction 
        # TODO: does it matter where the thiol group is? i.e. its bonding state? 
        radicals_positions = positions[(particle_ids==self.radical_thiol) | (particle_ids==self.radical_carbon)]
        radicals_ids = ids[(particle_ids==self.radical_thiol) | (particle_ids==self.radical_carbon)]

        if len(radicals_ids)>0:
            aq = freud.locality.AABBQuery(box, radicals_positions)
            nlist = aq.query(radicals_positions, {'r_max': self.rcut_neigh, 'exclude_ii':True}).toNeighborList()
            all_pairs = nlist[:]

            actual_ids_polymers = radicals_ids[all_pairs[:,0]]
            actual_ids_radicals = radicals_ids[all_pairs[:,1]]
            actual_ids = np.vstack((actual_ids_radicals,actual_ids_polymers)).T
            # sorting and making it unique to not make same bonds again 
            actual_ids = np.sort(actual_ids,axis=1)
            actual_ids = np.unique(actual_ids,axis=0)
            

            for a in np.unique(actual_ids[:,0]):

                neigh_a = actual_ids[actual_ids[:,0]==a][:,1]

                if len(neigh_a)>0:
                    b = np.random.choice(neigh_a)  
                   
                    bonds_on_b_1 = snapshot.bonds.group[snapshot.bonds.group[:,0]==b]
                    bonds_on_b_2 = snapshot.bonds.group[snapshot.bonds.group[:,1]==b]
                    bonds_on_b = np.unique(np.vstack((bonds_on_b_1,bonds_on_b_2)).flatten())
                    bonds_on_b = bonds_on_b[bonds_on_b!=b]

                    bonds_on_a_1 = snapshot.bonds.group[snapshot.bonds.group[:,0]==a]
                    bonds_on_a_2 = snapshot.bonds.group[snapshot.bonds.group[:,1]==a]
                    bonds_on_a = np.unique(np.vstack((bonds_on_a_1,bonds_on_a_2)).flatten())
                    bonds_on_a = bonds_on_a[bonds_on_a!=a]
                    
                    if self.FJ_system == False: 
                        U = np.where(snapshot.angles.typeid==self.dummy_type_angle)[0][0]
                        # good
                        snapshot.angles.group[U]=[bonds_on_b[0],b,a]  # angle 
                        snapshot.angles.typeid[U]=0
                        # good
                        snapshot.angles.group[U+1]=[bonds_on_a[0],a,b] #angle 
                        snapshot.angles.typeid[U+1]=0

                    # good 
                    U = np.where(snapshot.bonds.typeid==self.dummy_type_bond)[0][0]
                    snapshot.bonds.group[U]=[a,b]  # propagation 
                    snapshot.bonds.typeid[U]=1

                    print("radical anhiliation!")
                    print(a,b)
                    print(particle_ids[a],particle_ids[b])

                    # radicals anhiliate each other - both flip back
                    if  snapshot.particles.typeid[idx[a]]==self.radical_thiol:
                        snapshot.particles.typeid[idx[a]]=self.thiol_type  
                        particle_ids[a]=self.thiol_type

                    if  snapshot.particles.typeid[idx[b]]==self.radical_thiol:
                        snapshot.particles.typeid[idx[b]]=self.thiol_type   
                        particle_ids[b]=self.thiol_type 

                    if  snapshot.particles.typeid[idx[a]]==self.radical_carbon:
                        snapshot.particles.typeid[idx[a]]=self.ene_type  
                        particle_ids[a]=self.ene_type
                        
                    if  snapshot.particles.typeid[idx[b]]==self.radical_carbon:
                        snapshot.particles.typeid[idx[b]]=self.ene_type   
                        particle_ids[b]=self.ene_type 

                    print(particle_ids[a],particle_ids[b])
        self._state.set_snapshot(snapshot)
        


