
def checkprep(bbC_i,bbC_j,covbuffer=1.5,debug=False):
    '''
    Check 
    
    Arguments:
            bblockC_i (Container) Buildingblock container 1
            Xo_i (int) key of attachment point particle container 1
            bblockC_j (Container) Buildingblock container 2 
            Xo_j (int) key of attachment point particle container 2
    Retrun
        True - if no overlap
        False - if overlap is found 
    '''
    debug = False
    if( debug ):
        print " >checkprep debug on "
        
    # Set attachment points 
    Xo_i = bbC_i.attach_p
    Xo_j = bbC_j.attach_p
    
    npos_i = bbC_i.positions
    #npos_i.pop(Xo_i)
    npos_j = bbC_j.positions 
    #npos_j.pop(Xo_j)
    lat_i = bbC_i.lat
    npos_ij,nd_ij = lat_i.delta_npos(npos_i,npos_j)

    for pkey_i, particle_i  in bbC_i.particles.iteritems():
        if( pkey_i != Xo_i and particle_i.properties['number'] != 1 ):
            for pkey_j, particle_j  in bbC_j.particles.iteritems():
                if( pkey_j != Xo_j and particle_j.properties['number'] != 1 ):
                    dij = nd_ij[pkey_i][pkey_j]
                    radii_i = bbC_i.particles[pkey_i].properties["cov_radii"]
                    radii_j = bbC_j.particles[pkey_j].properties["cov_radii"]
                    cut_off = radii_i +radii_j
                    cut_off = cut_off*covbuffer
                    if( dij < cut_off ):
                        if( debug ):
                            log_line = "      >checkprep \n"
                            log_line += "           particle i %s %d \n"%(particle_i.properties['symbol'],pkey_i)
                            log_line += "           particle j %s %d \n"%(particle_j.properties['symbol'],pkey_j)
                            log_line += "           radii_i %f \n"%(radii_i)
                            log_line += "           radii_j %f \n"%(radii_j)
                            log_line += "           dij %f \n"%(dij)
                            print log_line
                        return False 
                    
    return True 
    
def shiftprep(bblockC_i,bblockC_j,debug=False):
    '''

    Concatenate two building block containers that had prepattach already ran on them 

        This is handy for batch attachments to not repeat the intial steps in the attach function
    
    Arguments:
            bblockC_i (Container) Buildingblock container 1
            Xo_i (int) key of attachment point particle container 1
            bblockC_j (Container) Buildingblock container 2 
            Xo_j (int) key of attachment point particle container 2
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            bbid_j   (str) Connecting atom bbid in container 2
            n_j      (int) number of connection to used in  container 2

        { bblockC_i - X_i - R_i } +  { R_j - X_j - bblockC_j }
                          =>
        { bblockC_i - X_i - X_j - bblockC_j }
    '''

    bbC_i = copy.deepcopy(bblockC_i)
    bbC_j = copy.deepcopy(bblockC_j)
    #
    Xo_i = bbC_i.attach_p
    Xo_j = bbC_j.attach_p
    # Shift  building block j to correct bond length 
    radii_i = bbC_i.particles[Xo_i].properties["cov_radii"]
    radii_j = bbC_j.particles[Xo_j].properties["cov_radii"]
    bond_vec = np.array([radii_i + radii_j,0.0,0.0])

    if( debug ):
        print " >attachprep Xo_i ",Xo_i
        print " >attachprep bbC_i.position[Xo_i] ", bbC_i.positions[Xo_i]
        print " >attachprep radii_i ",radii_i
        print " >attachprep Xo_j ",Xo_j
        print " >attachprep bbC_j.position[Xo_j] ", bbC_j.positions[Xo_j]
        print " >attachprep radii_j ",radii_j
    
    bbC_j.shift_pos(-1.0*bond_vec)

    return bbC_i,bbC_j
    
def attachprep(bbC_i,bbC_j,debug=False):
    '''

    Concatenate two building block containers that had prepattach already ran on them 

        This is handy for batch attachments to not repeat the intial steps in the attach function
    
    Arguments:
            bblockC_i (Container) Buildingblock container 1
            Xo_i (int) key of attachment point particle container 1
            bblockC_j (Container) Buildingblock container 2 
            Xo_j (int) key of attachment point particle container 2
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            bbid_j   (str) Connecting atom bbid in container 2
            n_j      (int) number of connection to used in  container 2

        { bblockC_i - X_i - R_i } +  { R_j - X_j - bblockC_j }
                          =>
        { bblockC_i - X_i - X_j - bblockC_j }
    '''
    # Set attachment points 
    Xo_i = bbC_i.attach_p
    Xo_j = bbC_j.attach_p
    #
    # Add j to i
    bbC_i += bbC_j
    #
    # Get updated atom key to form bond 
    Xp_j = bbC_i.keyupdate[Xo_j]
    #
    # Create bond between  X_i - X_j
    bond_ij = structure.Bond(Xo_i,Xp_j)
    bond_ij.properties['bondorder'] = 1
    if( debug ):
        print ">attachprep set bondorder ",Xo_i,Xp_j,bond_ij.properties['bondorder']
    bbC_i.add_bond(bond_ij)
    #
    # Remake neighbor list based on updated bonds 
    bbC_i.bonded_nblist = structure.NBlist() 
    bbC_i.bonded_nblist.build_nblist(bbC_i.particles,bbC_i.bonds )
    # Update number of attachment points
    bbC_i.calc_attachments()                

    return bbC_i
    
        
def attach(bblockC_i,bblockC_j,bbid_i="R",n_i=0,bbid_j="R",n_j=0,tag="blank"):
        '''

        Concatenate two building block containers

        Arguments:
            bblockC_i (Container) Buildingblock container 1 
            bblockC_j (Container) Buildingblock container 2 
            bbid_i   (str) Connecting atom bbid in container 1 
            n_i      (int) number of connection to used in  container 1 
            bbid_j   (str) Connecting atom bbid in container 2
            n_j      (int) number of connection to used in  container 2

        { bblockC_i - X_i - R_i } +  { R_j - X_j - bblockC_j }
                          =>
        { bblockC_i - X_i - X_j - bblockC_j }
        '''
        # Empty container checks
        if bblockC_j.n_particles == 0:             # If struc2 is empty (has no particles)
            return copy.deepcopy(bblockC_i)                 #   simply return unchanged current container
        if bblockC_i.n_particles == 0:              # If struc1 (this struc) is empty (has no particles)
            return copy.deepcopy(bblockC_j)                 #   simply return unchanged current container

        # If no tag for the new structure is provided c
        # ombine the tag of i and j
        if( tag == "blank" ):
            tag = bblockC_i.tag +  bblockC_j.tag
        
        Xn_i = 0  # Number of term atom in neighbor list of cap atom 
        Xn_j = 0  # Number of term atom in neighbor list of cap atom
        #
        # Make copies of containers to modify 
        bbC_i = copy.deepcopy(bblockC_i)
        bbC_i.tag = tag
        bbC_j = copy.deepcopy(bblockC_j)
        #
        # Record attachment 
        attachment_i = Attachment(tag,bblockC_i.tag,bblockC_j.tag,bbid_i,n_i,bbid_j,n_j)
        bbC_i.attachments.append(attachment_i)
        #
        # Find keys of attachment points 
        Rkey_i,Xkey_i = bbC_i.find_XR(bbid_i,n_i,Xn_i)
        Rkey_j,Xkey_j = bbC_j.find_XR(bbid_j,n_j,Xn_j)
        #
        # Sum charges of particles to be removed into attachment points
        bbC_i.sum_prop("charge",Xkey_i,Rkey_i)
        bbC_j.sum_prop("charge",Xkey_j,Rkey_j)
        # Align building blocks along bonds of attachment atoms 
        #
        bbC_i.align_bond(Rkey_i,Xkey_i)
        bbC_i.shift_pos(-1.0*bbC_i.positions[Xkey_i] )
        bbC_j.align_bond(Xkey_j,Rkey_j)
        #
        # Shift  building block j to correct bond length 
        radii_i = bbC_i.particles[Xkey_i].properties["cov_radii"]
        radii_j = bbC_j.particles[Xkey_j].properties["cov_radii"]
        bond_vec = np.array([radii_i + radii_j,0.0,0.0])
        bbC_j.shift_pos(-1.0*bond_vec)
        #
        # Remove atoms at R in building block i
        bbC_i.del_particle(Rkey_i)
        Xo_i = bbC_i.keyupdate[Xkey_i]
        bbC_j.del_particle(Rkey_j)
        Xo_j = bbC_j.keyupdate[Xkey_j]
        #
        # Add j to i
        bbC_i += bbC_j
        #
        # Get updated atom keys to form bond 
        Xp_i = Xo_i                   
        Xp_j = bbC_i.keyupdate[Xo_j]
        #
        # Create bond between  X_i - X_j
        bond_ij = structure.Bond(Xp_i,Xp_j)
        bbC_i.add_bond(bond_ij)
        #
        # Remake neighbor list based on updated bonds 
        bbC_i.bonded_nblist = structure.NBlist() 
        bbC_i.bonded_nblist.build_nblist(bbC_i.particles,bbC_i.bonds )
        # Update number of attachment points
        bbC_i.calc_attachments()
        # Set points of attachments
        bbC_i.Xp_i = Xp_i
        bbC_i.Xp_j = Xp_j
        #
        # Return final structure
        # 
        return bbC_i


def calc_cos(v_i,v_j):
    """
    Calculate the cos theta between two vectors 
    """
    return np.dot(v_i/np.linalg.norm(v_i),v_j/np.linalg.norm(v_j))
