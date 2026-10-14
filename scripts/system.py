import numpy as np
import gsd,gsd.hoomd
from collections import defaultdict
import freud


def molecule_counts(N_monomers, crosslinker_percent, monomer_size=0, extender_size=0):
    '''
    Number of molecules and beads for an equal-stoichiometry thiol-ene mixture.
    Each ene monomer has 2 C=C groups, so 2*N_monomers thiols are needed.
    crosslinker_percent of those thiols sit on tetrathiols (4 each), the rest on dithiols (2 each).
    monomer_size / extender_size = number of spacer beads in the di-ene / dithiol.
    Used by init.py and create_initial_configuration, so both always agree.
    '''
    N_monomers = int(N_monomers)
    reactive_beads = N_monomers*2
    N_crosslinkers = int(np.round(reactive_beads*crosslinker_percent/100./4.))
    N_extenders = int(np.round((reactive_beads-N_crosslinkers*4)/2.))
    beads_per_monomer = 4 + int(monomer_size)
    beads_per_extender = 2 + int(extender_size)
    beads_per_crosslinker = 5
    return dict(N_monomers=N_monomers, N_crosslinkers=N_crosslinkers, N_extenders=N_extenders,
                beads_per_monomer=beads_per_monomer, beads_per_extender=beads_per_extender,
                beads_per_crosslinker=beads_per_crosslinker,
                N_particles=N_crosslinkers*beads_per_crosslinker + N_extenders*beads_per_extender
                            + N_monomers*beads_per_monomer)


# make the system class
class System:
    def __init__(self):
        ## define chemical tags ----
        self.particles_types = ['Thiol','Ene_C','RSulfur','RCarbon', 'Sulfur', 'Carbon','C','D','Dummy']
        ## ---------------      [  0  , 1  ,     2     ,  3 ,     4     ,   5]
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

        # overwritten by "create_initial_configuration"
        self.N_particles = 0 
        self.N_dummy_bonds = 0 
        self.N_dummy_angles = self.N_dummy_bonds*2 
        
        self.size_monomer = 0 
        self.size_extender = 0 
        self.size_crosslinker = 0

        # overwritten by "propagate_reaction" 
        self.chain_transfer_probability = 0.0
        self.thiol_reaction_probability = 0.0
        self.chain_side_reaction_probability = 0.0
        self.rcut_neigh = 0

        self.FJ_system = False
        self.unspaced_system = True

    '''
    create the initial configuration in the box
    args: density, crosslinker percent, number of ene monomers,
          monomer_size  = spacer beads between the two C=C groups of the di-ene
          extender_size = spacer beads between the two thiols of the dithiol
          (0 and 0 give the simple molecules: E=E-E=E, S-S and a 4-arm tetrathiol)
    returns: initial gsd frame

    Spacer beads are type D, which never reacts. The tetrathiol is always
    4 thiols bonded to one central C bead.
    '''
    def create_initial_configuration(self,density,crosslinker,N_monomers, monomer_size=0, extender_size=0):
        counts = molecule_counts(N_monomers, crosslinker, monomer_size, extender_size)
        self.size_monomer = counts['beads_per_monomer']
        self.size_extender = counts['beads_per_extender']
        self.size_crosslinker = counts['beads_per_crosslinker']
        self.unspaced_system = (monomer_size == 0 and extender_size == 0)

        ### in this system, equal-stoichiometric ene-thiol
        self.N_monomers = counts['N_monomers']
        self.N_crosslinker = counts['N_crosslinkers']
        self.N_extenders = counts['N_extenders']
        reactive_beads = self.N_monomers*2

        print("N Monomers A = ",self.N_monomers)
        print("N Monomers B = ", self.N_extenders)
        print("N crosslinkers B =", self.N_crosslinker)
        print("Balance: 2*%s + 4*%s = 2*%s"%(self.N_extenders,self.N_crosslinker,self.N_monomers))

        # calculate the overall tallies in the box
        self.N_particles = counts['N_particles']
        self.L =  (self.N_particles/density)**(1/3.0)
        print("Box length = %1.2f"%(self.L))
        print("total particles = %d"%(self.N_particles))
        print("number density = %1.3f"%(self.N_particles/(self.L**3)))
        # spare bonds and angles that reactions turn into real ones (see propagate_reaction)
        self.N_dummy_bonds = int(np.round(reactive_beads*1.05))
        self.N_dummy_angles = self.N_dummy_bonds*3

        # one template per molecule: bead types, bonds, bond types, angles (indices within the molecule)
        # tetrathiol crosslinker: 4 thiols around one central C bead, arms swing freely (no angles)
        crosslinker_template = (
            [self.thiol_type]*4 + [self.spacer_type1],
            [[0,4],[1,4],[2,4],[3,4]],
            [0]*4,
            [])
        # dithiol extender: S - (extender_size x D) - S
        n = self.size_extender
        extender_template = (
            [self.thiol_type] + [self.spacer_type2]*extender_size + [self.thiol_type],
            [[i,i+1] for i in range(n-1)],
            [0]*(n-1),
            [[i,i+1,i+2] for i in range(n-2)] if self.FJ_system == False else [])
        # di-ene monomer: E=E - (monomer_size x D) - E=E, terminal double bonds, internal single bonds
        m = self.size_monomer
        monomer_template = (
            [self.ene_type]*2 + [self.spacer_type2]*monomer_size + [self.ene_type]*2,
            [[i,i+1] for i in range(m-1)],
            [2] + [4]*(m-3) + [2],
            [[i,i+1,i+2] for i in range(m-2)] if self.FJ_system == False else [])

        all_positions = []
        all_types = []
        all_bonds = []
        all_bonds_types = []
        all_angles = []
        N_current = 0
        for (types, bonds, bond_types, angles), count in [(crosslinker_template, self.N_crosslinker),
                                                           (extender_template, self.N_extenders),
                                                           (monomer_template, self.N_monomers)]:
            size = len(types)
            for i in range(count):
                # beads start within 0.5 of a random centre; equilibrate pulls them apart
                p = np.random.uniform(-self.L/2., +self.L/2., size=(1,3))
                offsets = self.sample_spherical(int(size))*0.5
                position = np.repeat(p,size,axis=0) + offsets.T
                start = N_current + size*i
                all_positions.extend(position)
                all_types.extend(types)
                all_bonds.extend([[a+start, b+start] for a, b in bonds])
                all_bonds_types.extend(bond_types)
                all_angles.extend([[a+start, b+start, c+start] for a, b, c in angles])
            N_current += size*count

        # format all molecular information
        all_positions = np.asarray(all_positions)
        all_positions = self.wrap_pbc(all_positions)
        all_types = np.asarray(all_types)
        all_bonds = np.asarray(all_bonds)
        all_bonds_types = np.asarray(all_bonds_types)
        all_angles = np.asarray(all_angles)

        # write all information into gsd file
        frame = gsd.hoomd.Frame()
        frame.particles.N = self.N_particles+4

        frame.particles.position = np.zeros(shape=(self.N_particles+4,3))
        frame.particles.position[:self.N_particles] = all_positions
        frame.particles.typeid = np.hstack((all_types, np.asarray([self.dummy_type,self.dummy_type,self.dummy_type,self.dummy_type])))

        frame.particles.types = self.particles_types
        frame.bonds.types = self.bond_types
        frame.angles.types = self.angle_types

        # throw the dummy particles in a corner of the box
        eps=1e-6
        Lo = self.L/2.0-eps
        frame.particles.position[-4:] = [Lo,Lo,Lo]

        frame.configuration.box = [self.L, self.L, self.L, 0, 0, 0]

        bonds = np.asarray((all_bonds.flatten()).reshape(-1,2))
        dummy_bonds = np.tile([self.N_particles,self.N_particles+1],(self.N_dummy_bonds,1))
        frame.bonds.group = np.vstack((bonds,dummy_bonds))

        frame.bonds.typeid = np.hstack((all_bonds_types,len(dummy_bonds)*[self.dummy_type_bond]))
        frame.bonds.N = len(bonds)+len(dummy_bonds)

        # if including angles, specify them in the gsd file
        if self.FJ_system == False:
            angles = np.asarray((all_angles.flatten()).reshape(-1,3))
            dummy_angles = np.tile([self.N_particles,self.N_particles+1,self.N_particles+2],(self.N_dummy_angles,1))
            frame.angles.group = np.vstack((angles,dummy_angles))
            frame.angles.typeid = np.hstack((len(angles)*[0],len(dummy_angles)*[self.dummy_type_angle]))
            frame.angles.N = len(angles)+len(dummy_angles)
        return frame
    
    
    
    
    '''
    "turn on the light"
    args: self, snapchot, %photonitiator (radicals)
    initiate the initial radicals on the thiols in the system 
    '''
    def flip_radicals_on(self,snapshot,radical_number_percent):
        # TODO: do we need to ensure only one thiol per monomer is selected or can both ends be reactive?
        '''
        # for only chain growth: starting with carbon radicals
        
        ids = np.arange(snapshot.particles.N) 
        bond_ids,bond_counts = np.unique((np.concatenate((snapshot.bonds.group[snapshot.bonds.typeid==2], snapshot.bonds.group[snapshot.bonds.typeid==4]),axis=0)).flatten(),return_counts=True)
        two_bond_ids = bond_ids[bond_counts==2]
        all_enes =  ids[snapshot.particles.typeid==self.ene_type]
        internal_enes = np.intersect1d(two_bond_ids,all_enes)
        n_radicals = int(np.round(radical_number_percent*len(internal_enes)/100.0))
        radicals = np.random.choice(internal_enes, n_radicals, replace=False)
        snapshot.particles.typeid[radicals] = self.radical_carbon
        for b in radicals:
            snapshot.bonds.typeid[snapshot.bonds.group[:,0]==b]=4
            snapshot.bonds.typeid[snapshot.bonds.group[:,1]==b]=4
        
        # normal
        # turn some A into radical_thiol
        '''
        ids = np.arange(snapshot.particles.N) 
        all_thiols =  ids[snapshot.particles.typeid==self.thiol_type]
        n_radicals = int(np.round(radical_number_percent*len(all_thiols)/100.0))
        radicals = np.random.choice(all_thiols, n_radicals, replace=False)
        snapshot.particles.typeid[radicals] = self.radical_thiol
        
        


    '''
    args:self, snapshot, box, particle_poitions, particle_ids, ids, idx
    thiol-ene propogation reaction: thiol radicals add terminally to ene group
    thiol radical (type 2) adds to terminal ene (type 1) that has only one bond
    >-<* + o=o-o=o = >-<o*=o-o=o
    0-2  + 1=1-1=1 = 0-0-3-1-1=1
    '''
    def bond_formation(self,snapshot,box,positions,particle_ids,ids,idx):

        if self.thiol_reaction_probability>0:
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
                    neigh_a = neigh_a[(neigh_a_types!=self.radical_thiol)&(neigh_a_types!=self.radical_carbon)]
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
                    
    
    '''
    args: self, snapshot, box, positions, paricle_ids, ids, idx
    chain transfer of thiol-ene reaction: H from SH abstracts to carbon radical
    effectively radical jumps from carbon to sulfur
    '''    
    def chain_transfer(self,snapshot,box,positions,particle_ids,ids,idx):
    
        # transfer carbon radical to new thiol 
        if self.chain_transfer_probability>0:
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
                    
                    neigh_a = neigh_a[(neigh_a_types!=self.radical_thiol)&(neigh_a_types!=self.radical_carbon)]
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
    def chain_growth(self,snapshot,box,positions,particle_ids,ids,idx):

        # chain side reaction radical carbon with ene reaction                
        if self.chain_side_reaction_probability>0:
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
                    
                    neigh_a = neigh_a[(neigh_a_types!=self.radical_thiol)&(neigh_a_types!=self.radical_carbon)]
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

                            '''
                            snapshot.angles.group[U]=[bonds_on_b[0],b,a]  # angle 
                            snapshot.angles.typeid[U]=1
                           
                            snapshot.angles.group[U+1]=[bonds_on_a[0],a,b] #angle 
                            snapshot.angles.typeid[U+1]=1
                            '''

                        
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
    '''
    '''
    def termination_reactions(self,snapshot,box,positions,particle_ids,ids,idx):    

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
    '''
    '''
    def propagate_reaction(self,
                           snapshot,r_cut=1.1,
                           chain_transfer_probability=1.0,
                           chain_side_reaction_probability=0.0,
                           thiol_reaction_probability=1.0):
        self.chain_transfer_probability=chain_transfer_probability
        self.chain_side_reaction_probability=chain_side_reaction_probability
        self.thiol_reaction_probability = thiol_reaction_probability 
        self.rcut_neigh = r_cut 

        box = freud.Box.from_box(snapshot.global_box)
        ids = np.arange(len(snapshot.particles.tag))

        idx = snapshot.particles.rtag[ids]
        particle_ids = snapshot.particles.typeid[idx]

        positions = snapshot.particles.position[idx]

        # init_radical_num = len(ids[particle_ids==self.radical_carbon]) + len(ids[particle_ids==self.radical_thiol])

        choices = ['propagation','chain_transfer','chain_growth']
        while len(choices) > 0:
            reaction = np.random.choice(choices)
            match reaction:
                case 'propagation':
                    # radical radical_thiol with Ene = Propagation/bond formation step 
                    rc_0 = ids[particle_ids==self.radical_carbon]
                    rt_0 = ids[particle_ids==self.radical_thiol]
                    self.bond_formation(snapshot,box,positions,particle_ids,ids,idx)
                    choices.remove('propagation')
                case 'chain_transfer':
                    # reaction radical_carbon with Thiol  - chain transfer step
                    self.chain_transfer(snapshot,box,positions,particle_ids,ids,idx)
                    ct_radical_num = len(ids[particle_ids==self.radical_carbon]) + len(ids[particle_ids==self.radical_thiol])

                    rc_ct = ids[particle_ids==self.radical_carbon]
                    rt_ct = ids[particle_ids==self.radical_thiol]
                    choices.remove('chain_transfer')
                case 'chain_growth':
                    #  competing chain growth radical carbond with ene reaction 
                    self.chain_growth(snapshot,box,positions,particle_ids,ids,idx)
                    cg_radical_num = len(ids[particle_ids==self.radical_carbon]) + len(ids[particle_ids==self.radical_thiol])
                    c_cg = ids[particle_ids==self.radical_carbon]
                    rt_cg = ids[particle_ids==self.radical_thiol]
                    # if (cg_radical_num < ct_radical_num):
                    #     print("chain growth error! initial radicals: " + str(init_radical_num) + ", fin radicals: " + str(cg_radical_num))
                    #     print("rc b4 : " + str(rc_ct) + "rt b4 : " + str(rt_ct))
                    #     print("rc : " + str(rc_cg) + "rt : " + str(rt_cg))
                    #     print("ct")
                    #     print("ct")
                    #     print("ct")
                    choices.remove('chain_growth')
                case _:
                    print("Error in propagate_reactions. Unknown reaction requested")

        # termination reactions, i.e radical-radical interactions 
        # self.termination_reactions(snapshot,box,positions,particle_ids,ids,idx)


    def sample_spherical(self,npoints, ndim=3):
        """
        Draw npoints random numbers on a sphere in ndim. 

        """
        vec = np.random.randn(ndim, npoints)
        vec /= np.linalg.norm(vec, axis=0)
        return vec
    

    def connected_components(self,lists):
        """
        merges lists with common elements

        args: bonds from configuration (frame.bonds.group)

        returns: list of connected particles by id

        Useful for finding bonded particles in a configuration, tested for
        linear polymers with consecutive bonds (0-1-2-3-4-5, 6-7-8-9-10,..)
        and non consecutive ids ( 0-5-6-8-10, 1-4-3-2-9,...) but no other
        configuration yet. Works with ints as well as str.

        """
        neighbors = defaultdict(set)
        seen = set()
        for each in lists:
            for item in each:
                neighbors[item].update(each)
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


    def wrap_pbc(self,delta):
        Box = np.array([self.L,self.L,self.L])
        delta = np.where(delta > 0.5 * Box, delta-Box, np.where(delta < -0.5 * Box, delta+Box , delta))
        return delta
    
