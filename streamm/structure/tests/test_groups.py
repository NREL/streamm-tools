

class TestGroupsProps(unittest.TestCase):
    # 
    def setUp(self):
        
        self.th = containers.Container("thiophene")
            
        symbols = ['C','C','C','C','S','H','H','H','H']
        positions = [ ]
        positions.append([-1.55498576,-1.91131218,-0.00081000])
        positions.append([-0.17775976,-1.91131218,-0.00081000])
        positions.append([0.34761524,-0.57904218,-0.00081000])
        positions.append([-0.65884476,0.36101082,0.00000000])
        positions.append([-2.16948076,-0.35614618,-0.00000800])
        positions.append([-2.18966076,-2.79526518,-0.00132100])
        positions.append([0.45389024,-2.80145418,-0.00106400])
        positions.append([1.41682424,-0.35961818,-0.00138200])
        positions.append([-0.51943676,1.44024682,0.00064700])
        for i in range(len(symbols)):
            pt_i = atoms.Atom(symbols[i])
            if( symbols[i] == 'C' ):
                pt_i.charge = -0.75
            elif(  symbols[i] == 'H' ):
                pt_i.charge = 0.80
            
            pt_i.mol = 0
            #pt_i.properties = periodictable.element_symbol()
            pos_i = positions[i]
            self.th.add_partpos(pt_i,pos_i)
        
        
        self.th.lat_cubic(100.0)
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.symbol == 'C' ):
                particle_i.resname = "SCP2"
                particle_i.residue = 1
            if( particle_i.symbol == 'S' ):
                particle_i.resname = "ThS"
                particle_i.residue = 2
            if( particle_i.symbol == 'H' ):
                particle_i.resname = "HA"
                particle_i.residue = 3
        self.strucC = containers.Container()
        self.strucC.lat_cubic(100.0)
        seed = 82343
        self.strucC = self.strucC.add_struc(self.th,10,seed,verbose=False)

    def test_molnumbers(self):
        for pkey_i, particle_i  in self.strucC.particles.iteritems():
            mol_i = int(pkey_i/self.th.n_particles) 
            self.assertEqual(str(particle_i.mol),str(mol_i))

    def test_groupmol(self):
        group_tag = 'mol'
        self.strucC.group_prop('mol',group_tag)
        groupset_i = self.strucC.groupsets[group_tag]
        self.assertEqual(str(len(groupset_i.groups)),str(10))
        groupset_i.calc_cent_mass()
        groupset_i.calc_cent_mass()
        
        cm = []
        cm.append('[ 61.022463  12.212374  55.579404]')
        cm.append('[ 94.589545   0.548985  40.058567]')
        cm.append('[ 13.025619  22.458819  96.090279]')
        cm.append('[ 45.93974   30.752004  73.031331]')
        cm.append('[ 28.945124  70.792119  10.476723]')
        cm.append('[ 26.732501  56.981684  23.793239]')
        cm.append('[ 48.205917  63.191955  94.038944]')
        cm.append('[ 28.343741  95.032088  28.668735]')
        cm.append('[ 83.906182   8.100332  26.885987]')
        cm.append('[ 97.987557  38.078986  85.843074]')

        groupset_i.calc_radius()
        groupset_i.calc_radius()
        groupset_i.calc_radius()
        for gkey,group_i in groupset_i.groups.iteritems():
            self.assertEqual(str(group_i.cent_mass),str(cm[gkey]))
            self.assertEqual(str(group_i.radius),'2.57775210944')
            self.assertEqual(str(group_i.r_gy_sq),'3.6041389371')
            # print "r_gy.append(\'%s\')"%str(group_i.properties)
        for gkey in groupset_i.keys:
            self.assertEqual(str(groupset_i.cent_mass[gkey]),str(cm[gkey]))
            self.assertEqual(str(groupset_i.radius[gkey]),'2.57775210944')
            self.assertEqual(str(groupset_i.r_gy_sq[gkey]),'3.6041389371')
            
        groupset_i.group_pbcs()

        os.chdir(os.path.dirname(__file__))
        groupset_i.write_cm_xyz()
        groupset_i.write_xyzs()
        groupset_i.dump_json()

        
    def test_groupres(self):
        group_tag = 'residue'
        self.strucC.group_prop('residue',group_tag)
        groupset_i = self.strucC.groupsets[group_tag]
        self.assertEqual(str(len(groupset_i.groups)),str(30))
        
        groupset_i.calc_cent_mass()
        groupset_i.calc_radius()

        #for gkey,group_i in groupset_i.groups.iteritems():
        self.assertEqual(round(groupset_i.radius[2],6),2.587885)
        self.assertEqual(round(groupset_i.r_gy_sq[2],6),4.967159)
        self.assertEqual(round(groupset_i.Q_mn[2][0][0],6),0.005067)
        self.assertEqual(round(groupset_i.Rgy_eignval[0][0],6),1.002185)
        self.assertEqual(round(groupset_i.Rgy_eignval[0][1],6),0.410354)
        self.assertEqual(round(groupset_i.Rgy_eignval[0][2],6),0.0)
        self.assertEqual(round(groupset_i.A_sphere[0],6),0.381661)
        self.assertEqual(round(groupset_i.A_sphere_num[0],6),1.52303)
        self.assertEqual(round(groupset_i.A_sphere_dem[0],6),1.995267)
        
        groupset_i.calc_dl()

        self.assertEqual(round(groupset_i.dl_sq[0],6),5.966521)
        self.assertEqual(round(groupset_i.dl_sq[2],6),20.729214)
                
        os.chdir(os.path.dirname(__file__))
        groupset_i.write_cm_xyz()
        groupset_i.write_xyzs()
        groupset_i.dump_json()
        

                
    def tearDown(self):
        del self.th         
        del self.strucC
        
class TestGroupsHtermSp2(unittest.TestCase):
    # 
    def setUp(self):
        self.th = structure.Container("thiophene")
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.th.tag)
        self.th.read_xyz(file_i)
        self.th.lat_cubic(100.0)
        # If no bonds guess based on radii 
        self.th.bonded_nblist.guess_nblist(self.th.lat,self.th.particles,self.th.positions,"cov_radii",radii_buffer=1.25)
        # Build bonds from nblist for reference 
        self.th.bonded_bonds()
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.properties['symbol'] == 'C' ):
                particle_i.properties['resname'] = "ThSC"
                particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'S' ):
                particle_i.properties['resname'] = "ThSC"
                particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'H' ):
                particle_i.properties['resname'] = "HA"
                particle_i.properties['residue'] = 3
        group_tag = 'residue'
        self.th.group_prop('residue',group_tag)
        self.groupset_i = self.th.groupsets[group_tag]
        self.assertEqual(len(self.groupset_i.groups),2)  
        
    def test_hterm(self):

        os.chdir(os.path.dirname(__file__))
        group_i =  self.groupset_i.groups[0]
        group_i.write_xyz('Th_SC.xyz')
        hterm_i = group_i.hterm_group()
        hterm_i.write_xyz('Th_SC_hterm.xyz')
        
    def tearDown(self):
        del self.th         

class TestGroupsHtermSp3(unittest.TestCase):
    # 
    def setUp(self):
        self.struc_i = structure.Container("ethane")
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.struc_i.tag)
        self.struc_i.read_xyz(file_i)
        self.struc_i.lat_cubic(100.0)
        # If no bonds guess based on radii 
        self.struc_i.bonded_nblist.guess_nblist(self.struc_i.lat,self.struc_i.particles,self.struc_i.positions,"cov_radii",radii_buffer=1.25)
        # Build bonds from nblist for reference 
        self.struc_i.bonded_bonds()
        #
    def test_hterm1(self):

        rmHcnt = 0 
        for pkey_i, particle_i  in self.struc_i.particles.iteritems():
            particle_i.properties['resname'] = 'CRES'
            particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'H' and rmHcnt < 1 ):
                particle_i.properties['resname'] = 'HRES'
                particle_i.properties['residue'] = 3
                rmHcnt += 1
                
        group_tag = 'residue'
        self.struc_i.group_prop('residue',group_tag)
        self.groupset_i = self.struc_i.groupsets[group_tag]
        self.assertEqual(len(self.groupset_i.groups),2)  
        
        os.chdir(os.path.dirname(__file__))
        group_i =  self.groupset_i.groups[0]
        group_i.write_xyz('Eth_C.xyz')
        hterm_i = group_i.hterm_group()

        for pkey_i, particle_i  in hterm_i.particles.iteritems():
            if( particle_i.properties['symbol'] == 'C' ):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),4)  
            if( particle_i.properties['symbol'] == 'H'):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),1)  
                        
        hterm_i.write_xyz('Eth_C_hterm1.xyz')
        
    def test_hterm2(self):

        rmHcnt = 0 
        for pkey_i, particle_i  in self.struc_i.particles.iteritems():
            particle_i.properties['resname'] = 'CRES'
            particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'H' and rmHcnt < 2 ):
                particle_i.properties['resname'] = 'HRES'
                particle_i.properties['residue'] = 3
                rmHcnt += 1
            
        group_tag = 'residue'
        self.struc_i.group_prop('residue',group_tag)
        self.groupset_i = self.struc_i.groupsets[group_tag]
        self.assertEqual(len(self.groupset_i.groups),2)  
        
        os.chdir(os.path.dirname(__file__))
        group_i =  self.groupset_i.groups[0]
        group_i.write_xyz('Eth_C.xyz')
        hterm_i = group_i.hterm_group()

        for pkey_i, particle_i  in hterm_i.particles.iteritems():
            if( particle_i.properties['symbol'] == 'C' ):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),4)  
            if( particle_i.properties['symbol'] == 'H'):
                self.assertEqual(hterm_i.bonded_nblist.calc_nnab(pkey_i),1)  
                        
        hterm_i.write_xyz('Eth_C_hterm2.xyz')
        
    def tearDown(self):
        del self.struc_i         


class TestGroup_dr(unittest.TestCase):
    # 
    def setUp(self):
        self.th = structure.Container("thiophene")
        file_i = os.path.join(os.path.dirname(__file__), "%s.xyz"%self.th.tag)
        self.th.read_xyz(file_i)
        self.th.lat_cubic(100.0)
        #
        for pkey_i, particle_i  in self.th.particles.iteritems():
            if( particle_i.properties['symbol'] == 'C' ):
                particle_i.properties['resname'] = "SCP2"
                particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'S' ):
                particle_i.properties['resname'] = "ThS"
                particle_i.properties['residue'] = 1
            if( particle_i.properties['symbol'] == 'H' ):
                particle_i.properties['resname'] = "HA"
                particle_i.properties['residue'] = 1
        self.strucC = structure.Container('th_x2')
        self.strucC.lat_cubic(30.0)
        seed = 82343
        self.strucC = self.strucC.add_struc(self.th,3,seed,verbose=False)
        self.strucC.tag = 'th_x3'
        #
    def test_finddr(self):
        self.strucC.lat_cubic(300.0)
        os.chdir(os.path.dirname(__file__))
        self.strucC.write_xyz()
        #self.strucC.write_cply()
        self.list_i = []
        for pkey,par_i in self.strucC.particles.iteritems():
            # print  pkey,par_i.properties['mol'],par_i.properties['symbol'] 
            if( par_i.properties['symbol'] == 'C' or par_i.properties['symbol'] == 'S' ):
                self.list_i.append(pkey)
                print pkey ,par_i.properties['mol'] , par_i.properties['symbol'] 
                
        # self.strucC.bonded_nblist.build_nblist(self.strucC.particles,self.strucC.bonds)
        self.strucC.bonded_nblist.guess_nblist(self.strucC.lat,self.strucC.particles,self.strucC.positions,"cov_radii",radii_buffer=1.25)
        group_id = 'mol'
        self.strucC.group_prop(group_id,group_id,particles_select=self.list_i)        
        
        groupset_i = self.strucC.groupsets[group_id]
        groupset_i.calc_cent_mass()
        groupset_i.write_cm_xyz()
        groupset_i.calc_radius()     

        list_i = groupset_i.groups.keys()
        list_j = groupset_i.groups.keys()        
        
        group_i = groupset_i.groups[0]
        group_i.write_xyz()
        
        #pairbuffer = 2.5        
        pairs_ij = groupset_i.find_pairs(list_i,list_j,mol_inter=True,mol_intra=False)
        
        
        r_cut = 25.0
        bin_size = 0.10            
        close_contacts = True
        
        n_bins = int(r_cut/bin_size) + 1 
        bin_r = np.zeros(n_bins)    
        bin_r_nn = np.zeros(n_bins)    # Nearest neighbor count 
        bin_r_pp = np.zeros(n_bins)    
        probabilityperpair = 1
        volumes = []
        #
        N_i = len(list_i)
        N_j = len(list_j)
        
        
        self.strucC.calc_volume()
        volumes.append(self.strucC.volume)
        
        npos_i = groupset_i.properties['cent_mass']
        npos_j = groupset_i.properties['cent_mass']        
        npos_ij,nd_ij = self.strucC.lat.delta_npos(npos_i,npos_j)
        
        
        for ref_i in range(N_i):
            a_i_hasnieghbor = False
            r_ij_nn = r_cut   # Nearest Neighbor distance  
            g_i = list_i[ref_i]
            for ref_j in range(N_j):
                if(  pairs_ij[ref_i][ref_j] > 0.0 ):
                    dr_ij =  nd_ij[ref_i,ref_j]
                    if(  dr_ij <= r_cut ):
                            # bin distance =
                            bin_index = int( round( dr_ij / bin_size) )
                            #
                            # print " dist / bin / bin_sit", dist[ref_i,ref_j],bin_index,bin_size*float(bin_index)
                            #
                            bin_r[bin_index] += probabilityperpair
                            # Find nearest neighbor distance 
                            a_i_hasnieghbor = True
                            if( dr_ij < r_ij_nn ):
                                r_ij_nn = dr_ij
                                p_ij_nn = pairs_ij[ref_i][ref_j]
                            # 
                            if( close_contacts ):
                                g_j = list_i[ref_j]
                                dr_pi_pj = groupset_i.dr_particles(g_i,g_j,r_cut)
                                bin_pp_index = int( round( dr_pi_pj / bin_size) )
                                bin_r_pp[bin_pp_index] += probabilityperpair

            # Record nearest neighbor distance 
            if( a_i_hasnieghbor ):
                bin_nn_index = int( round( r_ij_nn /bin_size) )
                bin_r_nn[bin_nn_index] += p_ij_nn  
        # 
        # print bin_r_nn
        # print bin_r_pp
        # 
        from datetime import datetime
        import math 
        
        n_frames = len(volumes)
        n_bins = len(bin_r)
        total_cnts = np.sum( bin_r )
        total_nn_cnts = np.sum( bin_r_nn )
        box_vol_ave = np.mean(volumes)
        # 
        cnt_sum_j = 0.0 
        nn_cnt_sum_j = 0.0
        pp_cnt_sum_j = 0.0 

        rdf = dict()
        rdf['index'] = []
        for key in ['r_val','g_r_box','g_r_nn_box','g_r_pp_box','nb_cnt','cnt_sum_j','nn_nb_cnt','nn_cnt_sum_j','pp_nb_cnt','pp_cnt_sum_j']:
            rdf[key] = []
            

        for bin_index in range(n_bins):
            r_val = bin_size*float(bin_index)
            dr_sq = r_val*r_val
            r_in = r_val - bin_size*0.5
            r_out = r_val + bin_size*0.5
            dr_vol = 4.0*math.pi/3.0*( r_out**3 - r_in**3 )
            cnt_r_frame = float( bin_r[bin_index] ) /float(n_frames) 
            nn_cnt_r_frame = float( bin_r_nn[bin_index] ) /float(n_frames)
            pp_cnt_r_frame = float( bin_r_pp[bin_index] ) /float(n_frames)
            # n(r)  = 1/N_i  sum_j^{N_j} \gamma( r - r_{ij}) 
            nb_cnt = cnt_r_frame/float( N_i )
            cnt_sum_j += nb_cnt
            nn_nb_cnt = nn_cnt_r_frame/float( N_i )
            nn_cnt_sum_j += nn_nb_cnt
            pp_nb_cnt = pp_cnt_r_frame/float( N_i )
            pp_cnt_sum_j += pp_nb_cnt
            # g(r) = <V> * n(r) / dV 
            g_r_box = box_vol_ave*nb_cnt/dr_vol/float( N_j )
            g_r_nn_box = box_vol_ave*nn_nb_cnt/dr_vol/float( N_j )
            g_r_pp_box = box_vol_ave*pp_nb_cnt/dr_vol/float( N_j )
            # 
            rdf['index'].append(bin_index)
            rdf['r_val'].append(r_val)
            rdf['g_r_box'].append(g_r_box)
            rdf['nb_cnt'].append(nb_cnt)
            rdf['cnt_sum_j'].append(cnt_sum_j)
            rdf['g_r_nn_box'].append(g_r_nn_box)
            rdf['nn_nb_cnt'].append(nn_nb_cnt)
            rdf['nn_cnt_sum_j'].append(nn_cnt_sum_j)
            rdf['pp_nb_cnt'].append(pp_nb_cnt)
            rdf['pp_cnt_sum_j'].append(pp_cnt_sum_j)
            rdf['g_r_pp_box'].append(g_r_pp_box)
            
            
        for i in rdf['index']:
            r_val = rdf['r_val'][i]
            g_r_box = rdf['g_r_pp_box'][i]
            nb_cnt = rdf['pp_nb_cnt'][i]
            cnt_sum_j = rdf['pp_cnt_sum_j'][i]
            
            if( nb_cnt > 0 ):
                print "total",r_val,nb_cnt,cnt_sum_j

        # 
        # Write data file 
        # 
        rdf_tag = 'TestGroup_dr'
        import json
        with open('%s.json'%(rdf_tag), 'w') as fp:
            json.dump(rdf, fp)
            
            
        print "color Labels Bonds black"
        print "color Display Background white"
        print "mol addfile {/Users/tkemper/Development/streamm-tools/tests/th_x3.xyz} type {xyz} first 0 last -1 step 1 waitfor 1 1"
        print "mol addfile {/Users/tkemper/Development/streamm-tools/tests/mol_0.xyz} type {xyz} first 0 last -1 step 1 waitfor 1 1"
        print "mol addfile {/Users/tkemper/Development/streamm-tools/tests/mol_cm.xyz} type {xyz} first 0 last -1 step 1 waitfor 1 1"
        
    def tearDown(self):
        del self.th         
        del self.strucC         

