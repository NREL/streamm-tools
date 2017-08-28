


import copy 


'''

    def guess_oplsa(self):
        """
        Guess OPLS-aa atom types based on coordination 
        """
        for pkey_i,particle_i in self.particles.iteritems():
            nb_cnt_i = self.bonded_nblist.calc_nnab(pkey_i)
            el_cnt_i = self.calc_elcnt(pkey_i,self.bonded_nblist)
            #
            # label carbons 
            #
            if particle_i.properties["number"] == 6 :
                if int(nb_cnt_i) == 4 :
                    particle_i.properties["fftype"] = 'CT' # Alkane
                if  int(nb_cnt_i) == 3 :
                    particle_i.properties["fftype"] = 'CA'  # Conjugated 
                if int(nb_cnt_i) == 2 :
                    particle_i.properties["fftype"] = 'C:'   # Allene


                if int(nb_cnt_i) == 1 :
                    particle_i.properties["fftype"] = '' # Aromatic C
                    error_line =  " WARNING!!! carbon index ",pkey_i," bonded to single atom "
                    sys.exit(error_line)
            #
            # label oxygens
            #
            if( particle_i.properties["number"] == 8 ):
                if int(nb_cnt_i) == 1 :
                    particle_i.properties["fftype"] = 'O' # double bonded
                if int(nb_cnt_i) == 2 :
                    particle_i.properties["fftype"] = 'OS' # ether

            #
            # label nitrogens 
            #
            if particle_i.properties["number"] == 7 :
                if int(nb_cnt_i) == 3 :      # amide
                    particle_i.properties["fftype"] = 'N' 

            #
            # label sulfurs
            #
            if( particle_i.properties["number"] == 16 ):
                if int(nb_cnt_i) == 2 :
                    particle_i.properties["fftype"] = 'S'   #  Thioether RSR (UA)


        #
        # label hydrogens
        #
        for pkey_i,particle_i in self.particles.iteritems():
            nb_cnt_i = self.bonded_nblist.calc_nnab(pkey_i)
            el_cnt_i = self.calc_elcnt(pkey_i,self.bonded_nblist)
            if( particle_i.properties["number"] == 1 ):
                if ( nb_cnt_i > 1  ):
                    sys.exit(' over coordinated H')
                if ( nb_cnt_i < 1  ):
                    error_line = ' unbonded H %d '%(pkey_i)
                    sys.exit(error_line)
                key_j = self.bonded_nblist.getnbs(pkey_i)[0]
                particle_j = self.particles[key_j]
                el_cnt_j = self.calc_elcnt(key_j,self.bonded_nblist)

                if ( particle_j.properties["fftype"]== 'CA' ):
                    particle_i.properties["fftype"] = 'HA' #
                if ( particle_j.properties["fftype"]== 'CT' ):
                    particle_i.properties["fftype"] = 'HC' #
        return 

    def test_guess_oplsa(self):
        self.strucC.guess_oplsa()
        self.assertEqual(self.strucC.particles[0].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[1].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[2].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[3].properties['fftype'],'CA')
        self.assertEqual(self.strucC.particles[4].properties['fftype'],'S')
        self.assertEqual(self.strucC.particles[5].properties['fftype'],'HA')
        self.assertEqual(self.strucC.particles[6].properties['fftype'],'HA')
        self.assertEqual(self.strucC.particles[7].properties['fftype'],'HA')
        self.assertEqual(self.strucC.particles[8].properties['fftype'],'HA')
        
'''


