
    def convert_gaussian(self,calc_gaussian):
        '''
        Convert Gaussian simulation to Buildingblock

        Input:
           calc_gaussian (Gaussian) simulation object

        Units:
            Bohr -> angstroms

        '''
        import gaussian
        
        if isinstance(calc_gaussian, gaussian.Gaussian):
            #
            # Convert Structure 
            #
            if( calc_gaussian.strucC.n_particles > 0):
                # Add lattice to GROMACS object
                matrix_o = calc_gaussian.strucC.lat._matrix 
                matrix_i = self.lat._matrix 
                for m in range(calc_gaussian.strucC.lat.n_dim):
                    for n in range(calc_gaussian.strucC.lat.n_dim):
                        matrix_i[m][n] = units.convert_bohr_ang(matrix_o[m][n] )
                self.lat.set_matrix(matrix_i)
                #
                # Add particles and positions to LAMMPS object 
                #
                for pkey_o, particle_o  in calc_gaussian.strucC.particles.iteritems():
                    particle_i = copy.deepcopy(particle_o)
                    pos_o = calc_gaussian.strucC.positions[pkey_o]       
                    pos_i = [ units.convert_bohr_ang(v_i) for v_i in pos_o ] 
                    self.add_partpos(particle_i,pos_i)

                if( len(calc_gaussian.strucC.bonded_nblist.list) > 0 ):
                    self.bonded_nblist = copy.deepcopy(calc_gaussian.strucC.bonded_nblist)
                if( len(calc_gaussian.strucC.nonbonded_nblist.list) > 0 ):
                    self.nonbonded_nblist = copy.deepcopy(calc_gaussian.strucC.nonbonded_nblist)
                    
                for bkey_o, bond_o  in calc_gaussian.strucC.bonds.iteritems():
                    bond_i = copy.deepcopy(bond_o)
                    self.add_bond(bond_i)
                for akey_o, angle_o in calc_gaussian.strucC.angles.iteritems():
                    angle_i = copy.deepcopy(angle_o)
                    self.add_angle(angle_i)
                for dkey_o,dih_o in calc_gaussian.strucC.dihedrals.iteritems():
                    dih_i = copy.deepcopy(dih_o)
                    self.add_dihedral(dih_i)
                for ikey_o,imp_o in calc_gaussian.strucC.impropers.iteritems():
                    imp_i = copy.deepcopy(imp_o)
                    self.add_improper(imp_i)
            else:
                raise TypeError("GAUSSIAN simulation structure contains no atoms")
                    