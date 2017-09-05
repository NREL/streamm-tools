
    def test_finddr(self):
        self.strucC.lat_cubic(300.0)
        os.chdir(os.path.dirname(__file__))
        self.strucC.write_xyz()
        #self.strucC.write_cply()
        self.list_i = []
        for pkey,par_i in self.strucC.particles.iteritems():
            # print  pkey,par_i.mol,par_i.symbol 
            if( par_i.symbol == 'C' or par_i.symbol == 'S' ):
                self.list_i.append(pkey)
                print pkey ,par_i.mol , par_i.symbol 
                
        # self.strucC.bonded_nblist.build_nblist(self.strucC.particles,self.strucC.bonds)
        self.strucC.bonded_nblist.guess_nblist(self.strucC.lat,self.strucC.particles,self.strucC.positions,"cov_radii",radii_buffer=1.25)
        group_id = 'mol'
        self.strucC.group_prop(group_id,group_id,particles_select=self.list_i)        
        
        groupContainer_i = self.strucC.groupContainers[group_id]
        groupContainer_i.calc_cent_mass()
        groupContainer_i.write_cm_xyz()
        groupContainer_i.calc_radius()     

        list_i = groupContainer_i.groups.keys()
        list_j = groupContainer_i.groups.keys()        
        
        group_i = groupContainer_i.groups[0]
        group_i.write_xyz()
        
        #pairbuffer = 2.5        
        pairs_ij = groupContainer_i.find_pairs(list_i,list_j,mol_inter=True,mol_intra=False)
        
        
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
        
        npos_i = groupContainer_i.cent_mass
        npos_j = groupContainer_i.cent_mass        
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
                                dr_pi_pj = groupContainer_i.dr_particles(g_i,g_j,r_cut)
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
        