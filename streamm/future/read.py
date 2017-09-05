
def parse_cplytag(self):
    '''
    Define bbid pased on <=v0.2 cplytag

    bb - X - T  - terminal site for polymerization 
    bb - X - R  - functionalizable site for functional groups 
    bb - X - S  - substitutable site for atomic substitutions 

    this should be flipped ...
    but it makes it easier to id the H since it's only bonded to 1 thing
     
    '''
    self.n_term = 0 
    self.n_func = 0
    self.n_sub = 0

    self.terms = []
    self.funcs = []
    self.subs = []
            
    for pkey_i, particle_i  in self.particles.iteritems():
       
        if( particle_i.cplytag[:8] == "termcap_" ):
            particle_i.rsite = "T"
            self.terms.append(pkey_i)
            self.n_term += 1 
        elif( particle_i.cplytag[:8] == "funccap_" ):
            particle_i.rsite = "R"
            self.funcs.append(pkey_i)
            self.n_func += 1 
        elif( particle_i.cplytag[:6] == "rgcap_" ):
            particle_i.rsite = "R"
            self.funcs.append(pkey_i)
            self.n_func += 1
        elif( particle_i.cplytag[:5] == "fcap_" ):
            particle_i.rsite = "S"
            self.subs.append(pkey_i)
            self.n_sub += 1
        elif( particle_i.cplytag[:5] == "term_" ):
            particle_i.rsite = "X"
        elif( particle_i.cplytag[:5] == "func_" ):
            particle_i.rsite = "X"
        else: 
            particle_i.rsite = particle_i.cplytag
            
        
    def read_cply(self,cply_file=''):
        """
        Read cply file
        """
        if( len(cply_file) == 0 ):
            cply_file = "%s.cply"%(self.tag)
        read_lattice = False
        read_attachments = False

        # Set tag of structure based on file
        if( self.tag == "blank" ):
            tag_s1 = cply_file.split("/")   # split string based on dir 
            tag_s2 = tag_s1[-1]             # take last string 
            self.tag = tag_s2[:-5]          # remove .cply
            self.common_tag = tag_s2[:-5]         # remove .cply
        #
        with open(cply_file) as f:
            for line in f:
                col = line.split()
                if( read_lattice ):
                    matrix = self.lat._matrix
                    matrix[lv_cnt][0] = float( col[0] )
                    matrix[lv_cnt][1] = float( col[1] )
                    matrix[lv_cnt][2] = float( col[2] )
                    if( lv_cnt == 2 ):
                        self.lat.set_matrix(matrix)
                        read_lattice = False 
                    lv_cnt += 1
                    
                elif( len(col) >= 4  and col[0] != "name" and col[0] != "bond"  and col[0] != "attachment"  and col[0] != "replication"  and col[0] != "#" and col[0] != "type"  and col[0] != "func" ):

                    BBatom_i = BBatom(str(col[0]))
                    
                    pos_i = [ float(col[1]),float(col[2]),float(col[3]) ]

                    if (len(col) >= 14 ):
                        BBatom_i.label = str(col[4])
                        BBatom_i.fftype = str(col[5])
                        BBatom_i.mass = float(col[6])                        
                        BBatom_i.charge  = float(col[7])
                        BBatom_i.qgroup = int(col[8])                        
                        BBatom_i.ring = int(col[9])                        
                        BBatom_i.residue = int(col[10])
                        BBatom_i.resname = str(col[11])                       
                        BBatom_i.mol = int(col[12])                        
                        BBatom_i.cplytag = str(col[13])

                    elif (len(col) == 13 ):
                        BBatom_i.label = str(col[4])
                        BBatom_i.fftype = str(col[5])
                        BBatom_i.mass = float(col[6])                        
                        BBatom_i.charge = float(col[7])
                        BBatom_i.qgroup = int(col[8])                        
                        BBatom_i.ring = int(col[9])                        
                        BBatom_i.residue = int(col[10])
                        BBatom_i.resname = str(col[11])                       
                        BBatom_i.mol = int(col[12])

                    elif (len(col) == 8 ):
                        BBatom_i.residue = int(col[5])
                        BBatom_i.resname = str(col[6])                        
                        BBatom_i.cplytag = str(col[7])
                        
                    elif (len(col) == 7 ):
                        BBatom_i.charge = float(col[4])
                        BBatom_i.residue = int(col[5])
                        BBatom_i.resname = str(col[6])

                    elif (len(col) == 5 ):
                        BBatom_i.cplytag = str(col[4])

                    self.add_partpos(BBatom_i,pos_i, deepcopy = True)
                   
                elif(len(col) >= 3 ):
                    if( col[0] == "bond"):
                        pkey1 = int(col[1]) - 1
                        pkey2 = int(col[2]) - 1
                        bond_i = Bond(pkey1,pkey2)
                        if( len(col) >= 4 ):
                            bond_i.border =  int(col[3]) 
                        #print "process_line bond line ",col
                        self.add_bond(bond_i)
                    if( col[0] == "name"):
                        self.tag = str(col[1])
                        self.common_tag = str(col[1])
                        self.deptag = str(col[2])
                        self.ctag = str(col[2])
                        self.name = str(col[3])
                        self.IUPAC = ''
                        for s in col[4:]:
                            self.IUPAC += str(s)
                    if( col[0] == "type"):
                        self.moltype = str(col[1])
                        self.backbone = str(col[2])
                        self.bblist = str(col[3:])

                # Key word search 
                elif( len(col) > 0 ):
                    if( str(col[0]) == 'lattice' ):
                        read_lattice = True
                        lv_cnt = 0                        
          
        # Convert old tags to new bbid 
        self.parse_cplytag()
        # Build neighbor list of bonded atoms 
        if( self.n_bonds == 0 ):
            # If no bonds guess based on radii
            self.bonded_nblist = self.guess_nblist(0,radii_buffer=1.25)
            # Build bonds from nblist for reference 
            self.bonded_bonds()
        else:
            self.bonded_nblist.build_nblist(self.particles,self.bonds )

    def write_cply(self,cply_file='', write_ff=True,write_bonds=True,write_latvec=True,write_attachments=True):
        """
        Write cply file
        """
        if( len(cply_file) == 0 ):
            cply_file = "%s.cply"%(self.tag)

        F = open(cply_file,'w')
        cply_line = "name %s %s %s %s \n"%(self.tag,self.deptag,self.name,self.IUPAC)
        cply_line += "type %s %s %s \n"%(self.moltype,self.backbone,self.bblist)
        
        if(write_ff ):
            cply_line += "# atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,ptclObj.mass,charge,qgroup,ring,residue,resname, mol, cplytag \n"
        else:
            cply_line += "# atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]), cplytag \n"
            
        F.write(cply_line)
        for pkey_i, particle_i  in self.particles.iteritems():
            pos_i = self.positions[pkey_i]
            atomic_symb = particle_i.element.symbol
            bbid = particle_i.rsite
            cplytag = particle_i.cplytag
            if(write_ff ):
                mass = particle_i.mass
                charge = particle_i.charge
                residue = particle_i.residue
                resname = particle_i.resname
                label = particle_i.label
                fftype = particle_i.fftype
                qgroup = particle_i.qgroup
                mol = particle_i.mol
                ring = particle_i.ring
                cply_line =  " %5s %16.8f %16.8f %16.8f %s %s %12.8f %12.8f  %d %d %d %s  %d %s \n"%(atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,mass,charge,qgroup,ring,residue,resname, mol, bbid )
                #cply_line =  "%d %5s %16.8f %16.8f %16.8f %s %s %12.8f %12.8f  %d %d %d %s  %d %s \n"%(pkey_i,atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]),label,fftype,mass,charge,qgroup,ring,residue,resname, mol, cplytag )
            else:
                cply_line =  " %5s %16.8f %16.8f %16.8f %s \n"%(atomic_symb ,float(pos_i[0]), float(pos_i[1]),float(pos_i[2]), cplytag )
            F.write(cply_line)

        if( self.lat._matrix[0][0] == 0.0 ):
            write_latvec = False 
        if( write_latvec ):
            F.write(" lattice \n")
            for d in range(3):
                cply_line = " %16.8f %16.8f %16.8f \n"%( self.lat._matrix[d][0],self.lat._matrix[d][1],self.lat._matrix[d][2])
                F.write(cply_line)
            

        if( write_bonds ):
            if( self.n_bonds > 0 ):
                for bkey_i, bond_i  in self.bonds.iteritems():
                    b_i = bond_i.pkey1 + 1 
                    b_j = bond_i.pkey2 + 1
                    # F.write("  bond %d %d  %d \n"%(b_i,b_j,bond_i.order']))
                    F.write("  bond %d %d  \n"%(b_i,b_j))

        if( write_attachments ):
            for attachment_i in self.attachments :
                cply_line = " attachment "
                cply_line += str(attachment_i) +"\n"
                F.write(cply_line)
            for replication_i in self.replications :
                cply_line = " replication "
                cply_line += str(replication_i) +"\n"
                F.write(cply_line)                
        F.close()

        return            