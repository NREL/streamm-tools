import os,sys,fileinput,random
from math import *

from basics import *
from fragments import *

class Options:
    def __init__(self, do_6_6 = True, do_5_6 = False, do_inversion = True, do_mono=False, do_bis=True, do_tris=False, 
                 with_equivalent=False, accuracy=1, nstates=12):
        self.do_6_6 = do_6_6
        self.do_5_6 = do_5_6
        self.do_inversion = do_inversion
        self.do_bis = do_bis
        self.do_tris = do_tris
        self.do_mono = do_mono
        self.with_equivalent = with_equivalent
        self.accuracy = accuracy
        self.nstates = nstates

class PosDesc:
    def __init__(self, bond="6-6", inverted=False):
        self.bond = bond
        self.inverted = inverted
    def to_string(self):
        if self.inverted: 
            c="B" 
        else: 
            c="A"
        s = "%s%s" % (self.bond, c)
        return s

class C60(Fragment):
    def __init__(self, fname=None, frag=None):
        Fragment.__init__(self, fname, frag)
        self.oxpos = []
        self.b56 = []
        self.b66 = []
        self.p5 = []
        self.p6 = []
        self.oxpos = []
        self.indene_pos = []
        self.indene_pos_idx = []
        self.com = []

        if (frag != None):
            self.b56 = frag.b56
            self.b66 = frag.b66
            self.p5 = frag.p5
            self.p6 = frag.p6
            self.oxpos = frag.oxpos
            self.indene_pos = frag.indene_pos
            self.indene_pos_idx = frag.indene_pos_idx
            self.com = frag.com

    def not_c60_part(self,test_idx):
        for i in test_idx:
            if i >59:
                return True
        return False


    def filter_add(self,test_idx):
        # given a plane (indices of points in test_idx), add it to the list of planes if it is not already in the list,
        # or add new points from test_idx to existing entry.
        new = True
        if (not self.connected(test_idx)):
            return False
        for i in range(0,len(self.planes)):
            if (is_subset(test_idx, self.planes[i])):
                new = False
            else:
                if (is_subset(self.planes[i], test_idx)):
                    new = False
                    for t in test_idx:
                        if (t not in self.planes[i]):
                            self.planes[i].append(t)
        if (self.not_c60_part(test_idx)):
            new = False
        if (new):
            self.planes.append(test_idx) # is not repeat

    def find_planes(self):
        # find all the planes in the structure
        self.planes = []
#        for i1 in range(0,len(self.atoms)):
        for i1 in range(0,60):
            a1 = self.atoms[i1]
            for j in range(0,len(a1.nb)):
                i2 = a1.nb[j]
                a2 = self.atoms[i2]
                for k in range(0,len(a2.nb)):
                    i3 = a2.nb[k]
                    if (i3 != i1):
                        a3 = self.atoms[i3]
                        plane_idx, plane_atm = self.add_coplanar([i1,i2,i3], [a1.pos,a2.pos,a3.pos])
#                        print plane_idx
                        self.filter_add(plane_idx)
#        print self.planes

    def sort_planes(self):
        # given the plane list, divide it into pentagons and hexagons
        self.p5 = []
        self.p6 = []
        for p in self.planes:
            if (len(p) == 5):
                self.p5.append(p)
            elif (len(p) == 6):
                self.p6.append(p)
            else:
                print "WARNING: plane ", p, " dangling."
        print len(self.p5), " pentagons"
        print len(self.p6), " hexagons"
        print "Plane indices"
        for p in self.p5:
            print p
        for p in self.p6:
            print p

    def find_shapes(self):
        # do the steps to find the pentagons and hexagons
        self.find_planes()
        self.sort_planes()

    def find_5_6_bonds(self):
        # given pentagon and hexagon lists, find all 5-6 bonds
        self.b56 = []
        for ip5 in range(0,len(self.p5)):
            p5 = self.p5[ip5]
            # search all _bonds_ in p5 (not just all pairs of points)
            for i1 in range(0,len(p5)):
                for i2 in range(i1+1, len(p5)):
                    if (self.bonded(p5[i1], p5[i2])):
                        # for each, search all hexagons to see if edge is shared
                        for ip6 in range(0,len(self.p6)):
                            p6 = self.p6[ip6]
                            if (is_subset([p5[i1], p5[i2]], p6)):
                                # found 5-6 bond
                                # but final check that it's not the bond below the "spinach"
#                                if (not not_under_spinach or (self.atoms[p5[i1]].no_non_c60_neighbors() and
                                ## COMMENTING THIS OUT FOR MORE GENERIC VERSION
#                                                              self.atoms[p5[i2]].no_non_c60_neighbors())):
                                idx_inter = get_subset([p5[i1], p5[i2]], p6)
                                self.b56.append([ip5,ip6,i1,i2,idx_inter[0], idx_inter[1]])
        print "5-6 bonds"
        for b in self.b56:
            print b, self.p5[b[0]], self.p6[b[1]]


    def find_6_6_bonds(self):
        # given hexagon list, find all 6-6 bonds
        self.b66 = []
           # codewise, stealing from above 
        for ip61 in range(0,len(self.p6)):
            p61 = self.p6[ip61]
            # search all _bonds_ in p5 (not just all pairs of points)
            for i1 in range(0,len(p61)):
                for i2 in range(i1+1, len(p61)):
                    if (self.bonded(p61[i1], p61[i2])):
                        # for each, search all hexagons to see if edge is shared
                        for ip62 in range(ip61+1,len(self.p6)):
                            p62 = self.p6[ip62]
                            if (is_subset([p61[i1], p61[i2]], p62)):
                                # found 6-6 bond
                                # but final check that it's not the bond below the "spinach"
#                                if (not not_under_spinach or (self.atoms[p5[i1]].no_non_c60_neighbors() and
                                ## COMMENTING THIS OUT FOR MORE GENERIC VERSION
#                                                              self.atoms[p5[i2]].no_non_c60_neighbors())):
                                idx_inter = get_subset([p61[i1], p61[i2]], p62)
                                self.b66.append([ip61,ip62,i1,i2,idx_inter[0], idx_inter[1]])
        print "6-6 bonds"
        for b in self.b66:
            print b, self.p6[b[0]], self.p6[b[1]]

    def find_oxy_pos(self):
        # from 5-6 bond list, find oxygen position
        # from above, we have b = [idx into p5, idx into p6, idx1 in particular p5, idx2 in particular p5]
        # note points in p5 accessed by idx1 and idx2 are ALSO in p6; that's the bond 
        # And entries of p5,p6 are themselves indices into atom list
        self.oxpos = []
        print "Ox positions:"
        for b in self.b56:
            p5 = self.p5[b[0]]
            atom1idx = p5[b[2]]
            atom2idx = p5[b[3]]
            mid = vec_ave(self.atoms[atom1idx].pos, self.atoms[atom2idx].pos)
            v1 = vec_sub(mid, self.com)
            oxpos = vec_axpy(ox_height, v1, self.com)
            print oxpos
            self.oxpos.append(oxpos)

    def find_indene_pos(self, do_5_6=False, do_6_6=True, do_inversion=True):
        # from 5-6 bond list, find oxygen position
        # from above, we have b = [idx into p5, idx into p6, idx1 in particular p5, idx2 in particular p5]
        # And entries of p5,p6 are themselves indices into atom list        
        self.indene_pos = []
        self.indene_pos_idx = [] # indices in global atom list of atoms under each indene pos
        self.indene_pos_type = []  # further description of this position (which bond type, fwd vs. inverse indene orientation)
        print "Indene anchor positions:"
        if (do_5_6):
            for b in self.b56:
                p5 = self.p5[b[0]]
                atom1idx = p5[b[2]]
                atom2idx = p5[b[3]]
                v1 = vec_sub(self.atoms[atom1idx].pos, self.com)
                pos1 = vec_axpy(indene_height, v1, self.com)
                v2 = vec_sub(self.atoms[atom2idx].pos, self.com)
                pos2 = vec_axpy(indene_height, v2, self.com)
                print pos1, pos2
                self.indene_pos.append([pos1, pos2])
                self.indene_pos_idx.append([atom1idx, atom2idx])
                self.indene_pos_type.append(PosDesc(bond="5-6", inverted=False))
                if (do_inversion):
                    self.indene_pos.append([pos2, pos1])
                    self.indene_pos_idx.append([atom2idx, atom1idx])
                    self.indene_pos_type.append(PosDesc(bond="5-6", inverted=True))
        if (do_6_6):
            for b in self.b66:
                p6 = self.p6[b[0]]
                atom1idx = p6[b[2]]
                atom2idx = p6[b[3]]
                v1 = vec_sub(self.atoms[atom1idx].pos, self.com)
                pos1 = vec_axpy(indene_height, v1, self.com)
                v2 = vec_sub(self.atoms[atom2idx].pos, self.com)
                pos2 = vec_axpy(indene_height, v2, self.com)
                print pos1, pos2
                self.indene_pos.append([pos1, pos2])
                self.indene_pos_idx.append([atom1idx, atom2idx])
                self.indene_pos_type.append(PosDesc(bond="6-6", inverted=False))
                if (do_inversion):
                    self.indene_pos.append([pos2, pos1])
                    self.indene_pos_idx.append([atom2idx, atom1idx])
                    self.indene_pos_type.append(PosDesc(bond="6-6", inverted=True))

    def add_indene(self, bond_idx, indene):
        pos = self.indene_pos[bond_idx]
        # for debuging: attach extra atoms at "pin" (connector atom) locations
        self.atoms.append(Atom("Al", pos[0]))
        self.atoms.append(Atom("Al", pos[1]))
        frag = self.prepare_fragment(bond_idx, indene)
        print "Frag: len ", len(frag.atoms)
        for a in frag.atoms:
            self.atoms.append(a)
            print a.elmnt, a.pos
    

    def prepare_indene_set(self, bond_idx_set, indene):
        rot_frags = []
        min_allowed_interfrag_dist= 1.0
#        min_allowed_interfrag_dist= 0.0  # turns this off except for trying to attach two indenes at root with opposite orientation
        for i in range(0,len(bond_idx_set)):
            # for each bond location: 
            pos = self.indene_pos[bond_idx_set[i]]
            # for debuging: attach extra atoms at "pin" (connector atom) locations
#            self.atoms.append(Atom("Al", pos[0]))
#            self.atoms.append(Atom("Al", pos[1]))
            # make a transformed fragment
            frag = self.prepare_fragment(bond_idx_set[i], indene)
            # check whether it's compatible with previously prepared fragments
            for j in range(0,i):
                d = frag.dist_to_other_structure(rot_frags[j])
                print "min interfrag distance %d to %d is %f" % (bond_idx_set[i], bond_idx_set[j], d)
                if (d < min_allowed_interfrag_dist):
                    # report that we failed
                    return None
            # if so, add to list
            rot_frags.append(frag)

#        print "rot frags", rot_frags

        # report that we succeeded via list of "placed" fragments
        return rot_frags

    def prepare_fragment(self,bond_idx, indene):
        ## This function handles fairly generically the connection of 2 atoms on one fragment to 2 atoms on another.
        ## this leaves a degree of freedom, which we handle by scanning.  in this case the criteria is keeping
        ## the placed fragment away from the center of mass.
        # real job is transform indene "connector atoms" to those identified here, and all other indene atoms
        pos = self.indene_pos[bond_idx]
        # get transform for given connector points
        TR = get_2pt_transform(pos[0], pos[1], indene.atoms[indene.pin_idx[0]].pos, indene.atoms[indene.pin_idx[1]].pos) 
        # but there is still a degree of freedom:
        free_axis = vec_normalize(vec_sub(pos[0],pos[1]))
        nj = 10
        maxmindist = 0
        mintheta = None
        # pick by scan of different angles, pick one by _maximizimng_ minimal distance of atoms in 
        # rotated indene from center of mass of c60
        for j in range(0,nj):
            mindist = 1000
            theta = 2.0 * pi * float(j) / float(nj)
            R2 = get_axis_angle_rotation(free_axis, theta)
            for i in range(0,len(indene.atoms)):
                a = indene.atoms[i]
                ta = apply_transforms(TR, R2, a, pos[0])
                d = vec_norm(vec_sub(ta, self.com))
                mindist = min(d, mindist)
            if (mintheta == None or mindist > maxmindist):
                mintheta = theta
                maxmindist = mindist

#        print "mintheta, maxmindist ", mintheta, maxmindist
        R2 = get_axis_angle_rotation(free_axis, mintheta)

        # now, with optimal orientation of indene, copy 
        # to a new fragment
        frag = Fragment()
        for i in range(0,len(indene.atoms)):
            a = indene.atoms[i]
            ta = apply_transforms(TR, R2, a, pos[0])
            e = a.elmnt
## debug: highlight attachment point atoms in VMD.            
#            if (i in indene.pin_idx):
#                e = "O"
            frag.atoms.append(Atom(e, ta))
        
        return frag

    def shift_indene_carbons(self, idx=None):
        # shift in place, so no structure re-use for different indenizations
        if (idx == None):
            idx = range(0,len(self.indene_pos_idx))
        for i in idx:
            for j in range(0,2):
                atomidx = self.indene_pos_idx[i][j]
                v = vec_sub(self.atoms[atomidx].pos, self.com)
                p = vec_axpy(indene_carbon_height, v, self.com)
#                print atomidx, self.atoms[atomidx].pos, self.com, v, p
                self.atoms[atomidx].pos = p
            

    def shift_carbons(self, idx=None):
        ## for output, shift the carbons in the bond in question
        if (idx == None):
            idx = range(0,len(self.oxpos))
        self.shifted_atoms = []
        for a in self.atoms:
            self.shifted_atoms.append(Atom(a.elmnt, a.pos))

        # go back through 5-6 bonds and shift carbons a little
        # from above, we have that oxpos order is same as b56 order, so idx can refer to either.
#        for b in self.b56:
        for i in idx:
            #old code
            b = self.b56[i]
            p5 = self.p5[b[0]]
            atom1idx = p5[b[2]]
            atom2idx = p5[b[3]]

            # new logic: bonds may be 5-6 list or 6-6 list, either way 2nd of pair is in self.p6.
            # indices of points in particular p6 plane are in b[4] and b[5]
#            b = bond_list[i]
#            p5 = self.p5[b[0]]
#            atom1idx = p5[b[4]]
#            atom2idx = p5[b[5]]

            mid = vec_ave(self.atoms[atom1idx].pos, self.atoms[atom2idx].pos)
            # now bump atom1 and atom2 outward a little
            v1 = vec_sub(self.atoms[atom1idx].pos, mid)
            v2 = vec_sub(self.atoms[atom2idx].pos, mid)
            self.shifted_atoms[atom1idx].pos = vec_axpy(carbon_expansion_factor, v1, mid)
            self.shifted_atoms[atom2idx].pos = vec_axpy(carbon_expansion_factor, v2, mid)

    def get_idx_set_num_5_6(self, bond_idx_set):
        cnt66=0
        cnt56=0
        for b in bond_idx_set:
            if (self.indene_pos_type[b].bond == "6-6"):
                cnt66 += 1
            if (self.indene_pos_type[b].bond == "5-6"):
                cnt56 += 1
        return cnt56

    def write_oxpos_subset(self, idx=None, fname="oxpos.xyz"):
        # write the results
        f = file(fname, "w")
        if (idx == None):
            idx = range(0,len(self.oxpos))
        f.write("%d\noxpos\n" % len(idx))
        for i in idx:
            ox = self.oxpos[i]
            f.write("O %f %f %f\n" % (ox[0],ox[1],ox[2]))
        f.close()

    def write_gaussian_input(self, idx, i, M, N, dirname):
        ### this is for the oxidizer case, circa late 2009 (?) ###
        self.shift_carbons(idx)

        fname = "subsets/%s/%dO_%dof%d_PCBM.inp" % (dirname, M, i+1, N)
        f = file(fname, "w")
        """
        %chk=3O_7of10_PCBM.chk
        # opt b3lyp/6-31g*
        
        Name of the job (can be similar to the chk file title)
        
        0 1
        """
        f.write("%%chk=%dO_%dof%d_PCBM.chk\n" % (M,i+1,N))
        f.write("# opt b3lyp/6-31g*\n\n")
        f.write("%dO_%dof%d_PCBM\n\n0 1\n" % (M, i+1, N))
        
        # pcbm part:
        for a in self.shifted_atoms:
            f.write("%s %f %f %f\n" % (a.elmnt, a.pos[0], a.pos[1], a.pos[2]))

        #oxygen part:
        for i in idx:
            ox = self.oxpos[i]
            f.write("O %f %f %f\n" % (ox[0],ox[1],ox[2]))

        f.close()

        # side effect/debug: write xyz file
        f = file("subsets/%s/full_structure.xyz" % (dirname), "w")
        f.write("%d\ntest\n" % (len(self.shifted_atoms) + len(idx)))
        for a in self.shifted_atoms:
            f.write("%s %f %f %f\n" % (a.elmnt, a.pos[0], a.pos[1], a.pos[2]))
        #oxygen part:
        for i in idx:
            ox = self.oxpos[i]
            f.write("O %f %f %f\n" % (ox[0],ox[1],ox[2]))
        f.close()


def one_dist_set(mids, indenes, i0, i1):
    ctheta = vec_dot(indenes[i0],indenes[i1]) / (vec_norm(indenes[i0]) * vec_norm(indenes[i1])) 
    tip_dist = vec_norm(vec_sub(vec_add(mids[i0],indenes[i0]),vec_add(mids[i1],indenes[i1])))
    bond_dist = vec_norm(vec_sub(mids[i0], mids[i1]))
    return ctheta, tip_dist, bond_dist

def build_meta_dict(c60, bond_idx_set, job_name):
    # details unavoidable.
    # alpha = angle between first two indenes
    # beta = angle betwee their plane and 3rd indene. YOW, good fun!
    ## need midpoints between all pairs in c60.indene_pos[bond_idx_set]
    ## first two define alpha, 3rd beta.
    # we also calculate "gamma"s, pairwise angles
    # might as well store the bond index set also
    mids = []
    bonds = []
    indenes = []
    dall = {}
    dm = {}
    dg = {}
    dall['metadata'] = dm
    dall['geometry'] = dg

    dm['name'] = job_name
    dm["td_method"] = "B3LYP"
    dm["opt_method"] = "B3LYP"
    dm["opt_basis"] = "6-31G(d)"
    dm["td_basis"] = "6-31G(d)"
    dm["nstates"] = "12"

    s = ""
    cnt66 = 0
    cnt56 = 0
    for b in bond_idx_set:
        p1 = vec_sub(c60.indene_pos[b][0],c60.com)
        p2 = vec_sub(c60.indene_pos[b][1],c60.com)
        mid = vec_ave( p1 , p2)
        mids.append(mid)
        s +=  "%d:%s " % (b, c60.indene_pos_type[b].to_string())
        bond = vec_sub(p2,p1)
        bonds.append(bond)
        ind = vec_cross(mid, bond)
        indenes.append(ind)
        if (c60.indene_pos_type[b].bond == "6-6"):
            cnt66 += 1
        if (c60.indene_pos_type[b].bond == "5-6"):
            cnt56 += 1
    # try to get at symmetry via various vectors, distances, and angle.  bis case first:      
    cosa = vec_dot(mids[0],mids[1]) / (vec_norm(mids[0]) * vec_norm(mids[1])) 
    a = acos(cosa)
    dg['alpha'] = a
    dg['num_6-6'] = cnt66
    dg['num_5-6'] = cnt56
    
    ctheta, tip_dist, bond_dist = one_dist_set(mids, indenes, 0,1)
    dg['cos_theta'] = ctheta         
    dg['tip_dist'] = tip_dist
    dg['theta'] = acos(min(1.0,max(-1.0,ctheta)))
    dg['bis_dist'] = bond_dist

    # alpha is for bis- case, just angle between center of c60 and the two indene attachment points
    if (len(bond_idx_set) > 2):
        ### first section is original stuff Nikos et al suggested
        #  for tris- case, we store the separate angles between each pair
        cosg02 = vec_dot(mids[0],mids[2]) / (vec_norm(mids[0]) * vec_norm(mids[2])) 
        g02 = acos(cosg02)
        cosg12 = vec_dot(mids[1],mids[2]) / (vec_norm(mids[1]) * vec_norm(mids[2])) 
        g12 = acos(cosg12)
        dg['gamma01'] = a
        dg['gamma02'] = g02
        dg['gamma12'] = g12
        # and store the angle between the 3rd and the plane made by the first two
        n = vec_cross(mids[0], mids[1])
        cosb = vec_dot(n,mids[2]) / (vec_norm(n) * vec_norm(mids[2])) 
        b = acos(cosb)
        # now b is angle between 3rd and _normal_ to plane, so angle we want is 90-b
        dg['beta'] = pi/2.0 - b

        ### this stuff if analog of what in bis case proved necessary for symmetry checking
        idxs = [[0,1],[0,2],[1,2]]
        # theta = angle between oriented indenes
        # tip_dist = distance between oriented indene tips
        # tris_dist = distance between centers of c60 bonds that indenes are attached to
        for iset in idxs:
            ctheta, tip_dist, bond_dist = one_dist_set(mids, indenes, iset[0], iset[1])
            dg['cos_theta%d%d' % (iset[0], iset[1])] = ctheta         
            dg['tip_dist%d%d' % (iset[0], iset[1])] = tip_dist
            dg['theta%d%d' % (iset[0], iset[1])] = acos(min(1.0,max(-1.0,ctheta)))
            dg['tris_dist%d%d' % (iset[0], iset[1])] = bond_dist

    dg['bond indices'] = s
    dg['bond_idx_set'] = bond_idx_set
    return dall
        

########################
# hack / globals for now to handle removal of equivalent structures
gPosDict = []
gBondIdxSets = []

def add_to_global_pos_dict(dall):
    d = dall['geometry']

    nindene = d['num_6-6'] + d['num_5-6']
    if (nindene == 2):
        # bis case
        gPosDict.append([ d['bis_dist'], d['cos_theta'], d['tip_dist']])
    elif (nindene == 3):
        gPosDict.append([ d['tris_dist01'], d['tris_dist02'], d['tris_dist12'], d['theta01'], d['theta02'], d['theta12'],
                          d['tip_dist01'], d['tip_dist02'], d['tip_dist12']])        

    gBondIdxSets.append(d['bond_idx_set'])

def clear_pos_dict():
    del gPosDict[:]
    del gBondIdxSets[:]

def permute_pd_entry(pde, perm):
    # reorder by sets of three
    res = [pde[i] for i in range(0,len(pde))]
    for i in range(0,len(perm)):
        pi = perm[i]
        for j in range(0,3):
            res[pi*3 + j] = pde[3*i + j]
#    print "permuting ", pde, "  to  ", res 
    return res

def cluster_pd_perm(pd, bidxs, perm):
    # a list of cosines and distances between indene attachment points and between indene tips
    # for tris case, there are three sets of three corresponding to different pairs. 
    pos_dict_dist_thresh = 0.01   
    m = []
    wt = [1.0 for i in range(0,len(pd[0]))]
    mm = [[100,-100] for i in range(0,len(wt))]
    for i in range(0,len(pd)):
        for j in range(0,len(wt)):
            mm[j][0] = min(mm[j][0], pd[i][j])
            mm[j][1] = max(mm[j][1], pd[i][j])
    for j in range(0,len(wt)):
        wt[j] = abs(1/(mm[j][1] - mm[j][0]))

#    print "using weight ", wt
    for i in range(0,len(pd)):
        row = []
        for j in range(0,len(pd)):
            pde = permute_pd_entry(pd[j], perm)
            row.append(vec_wtnorm(wt, vec_sub(pd[i],pde)))
        m.append(row)
#    print "pos dict distances:"
#    print m
    clusters = []
    for i in range(0,len(pd)):
        k = 0
        found = False
        while (k < len(clusters) and not found):
            if (m[i][clusters[k][0]] < pos_dict_dist_thresh):    # ith pos close to kth cluster, exemplified by 0th elmnt 
                clusters[k].append(i)
                #print "appending ", i, "  to cluster ", k
                found = True
            else:
                k += 1
        if (not found):
            #print "not found, append new cluster starts with ", i
            clusters.append([i])
    print "For permutation ", perm, "  ", len(clusters), " inequivalent clusters found: ", clusters
    return clusters

def cluster_pos_dict(pd, bidxs, base_name, full_set_dir, set_dir):
    import glob
    from string import replace
    print "pos_dict"
    for i in range(0,len(pd)):
        print i,pd[i], bidxs[i]
    print "Clustering equivalent structures"
    clusters = None
    if (len(pd[0]) == 3):
        perms = [[0]]
    else:  ## this assumes tris case!
        perms = [[0,1,2], [0,2,1], [1,0,2], [1,2,0], [2,0,1], [2,1,0]]
    for perm in perms:
        new_clusters = cluster_pd_perm(pd, bidxs, perm)
        if (clusters == None or len(new_clusters) < len(clusters)):
            clusters = new_clusters

    print "Deleting redundant structures."
    bond_idx_list = []
    for cl in clusters:
#        print cl
        print "cluster ", cl
        for i in range(0,len(cl)):
            print pd[cl[i]], bidxs[cl[i]]
        for i in range(0,len(cl)):  # keep 0th, delete others
            xyz_name, job_name = get_names(bidxs[cl[i]], base_name, full_set_dir, set_dir)
            if (i==0):
                #print "keeping: ", xyz_name, job_name
                bond_idx_list.append(bidxs[cl[i]])
            else:
                root_name = replace(xyz_name, ".xyz", "")
                star_name = "%s.*"  % root_name
                for filename in glob.glob(star_name) :
#                    print "removing: ", filename
                    os.remove( filename ) 
    return bond_idx_list



###############################################

def get_names(bond_idx_set, base_name, full_set_dir, set_dir):
    s = ""
    for p in bond_idx_set:
        s += ".%d" % p
    xyz_name = "%s/%s.indene%s.xyz" % (full_set_dir, base_name, s)
    job_name = "%s.indene%s" % (base_name, s)
    return xyz_name, job_name
            

def finish_one_frag_set(c60, bond_idx_set, indene, opts, full_set_dir =".", set_dir=""):
    placed_indenes = c60.prepare_indene_set(bond_idx_set, indene)    
    if placed_indenes == None:
        print "indene attachment FAILED for bond index set", bond_idx_set
        return 0
    else:
        print "indene attachment SUCCEEDED for bond index set", bond_idx_set
        decor = C60(frag=c60)
        decor.append_fragments(placed_indenes)
        decor.shift_indene_carbons(bond_idx_set)
        xyz_name, job_name = get_names(bond_idx_set, decor.fname, full_set_dir, set_dir)
        decor.write_xyz(xyz_name)
        decor.write_com("run.com.template", xyz_name, job_name, get_basis_str(opts.accuracy), opts.nstates)
#        frag.write_com("donoracceptor.com.template",  xyz_name, job_name, get_basis_str(options.accuracy), options.nstates)
        decor.write_slurm("run.slurm.template", xyz_name, job_name)
#        frag.write_slurm("donoracceptor.slurm.template",  xyz_name, job_name)
        d = build_meta_dict(c60, bond_idx_set, job_name)
        decor.write_meta(d, xyz_name)
        add_to_global_pos_dict(d)
        return 1

def prepare_structs_and_dirs(c60, indene, opts, mols_dir):
    # will prepare one of 6-6, 5-6, or mixed dirs depending on opts, subsequent
    # assumed contents of c60.indene_pos
    import os.path, os
    if (not os.path.isdir(mols_dir)):
        os.mkdir(mols_dir)
    if (opts.do_6_6 and opts.do_5_6):
        set_dir = "mixed"
    elif (opts.do_6_6):
        set_dir = "6-6"
    elif (opts.do_5_6):
        set_dir = "5-6"
    else:
        print "NOTHING TO PREPARE"
        return
    full_set_dir = "%s/%s" % (mols_dir, set_dir)
    if (not os.path.isdir(full_set_dir)):
        os.mkdir(full_set_dir)

    use_reverse_as_base_too = True  ## this is necessary in 5-6 case especially
    totalgood = 0
    totaltry = 0
    if (opts.do_mono):
        bond_idx_set = [0]
        totalgood += finish_one_frag_set(c60, bond_idx_set, indene, opts, full_set_dir, set_dir)
        totaltry += 1

    bis_bond_idx_list = None
    if (opts.do_bis):
        for i in range(1, len(c60.indene_pos)):
            bond_idx_set = [0, i]
            if (set_dir != "mixed" or c60.get_idx_set_num_5_6(bond_idx_set) == 1):
                totalgood += finish_one_frag_set(c60, bond_idx_set, indene, opts, full_set_dir, set_dir)
                totaltry += 1
                if (opts.do_inversion and use_reverse_as_base_too and i > 1):
                    bond_idx_set = [1, i]
                    totalgood += finish_one_frag_set(c60, bond_idx_set, indene, opts, full_set_dir, set_dir)
                    totaltry += 1
        print "There were %d successful bis indene set placements out of %d tries" % (totalgood, totaltry)
        if (not opts.with_equivalent):
            bis_bond_idx_list = cluster_pos_dict(gPosDict, gBondIdxSets, c60.fname, full_set_dir, set_dir)
    
    clear_pos_dict()
    print "GLOBAL POS DICT NOW HAS LENGTH ", len(gPosDict)
    totaltry = 0
    totalgood = 0
    if (opts.do_tris):
        if (bis_bond_idx_list != None):
            for i in range(0, 10):
#            for i in range(0, len(bis_bond_idx_list)):
                bi = bis_bond_idx_list[i]
                for j in range(0, 20):
#                for j in range(0, len(c60.indene_pos)):
                    bond_idx_set = [bi[0], bi[1], j]
                    totalgood += finish_one_frag_set(c60, bond_idx_set, indene, opts, full_set_dir, set_dir)
                    totaltry += 1
        else:
            for i in range(1, len(c60.indene_pos)):
                for j in range(i, len(c60.indene_pos)):
                    bond_idx_set = [0, i, j]
                    totalgood += finish_one_frag_set(c60, bond_idx_set, indene, opts, full_set_dir, set_dir)
                    totaltry += 1

        print "There were %d successful tris indene set placements out of %d tries" % (totalgood, totaltry)
        if (not opts.with_equivalent):
            tris_bond_idx_list = cluster_pos_dict(gPosDict, gBondIdxSets, c60.fname, full_set_dir, set_dir)


def enumeration_tests(c60, indene, opts, mols_dir):
    # test 1: can we just place a SET once
#    bond_idx_set = [0, 13, 2]
    # for just 0, this is the mono- case
    if (True):
        bond_idx_set = [0]
        finish_one_frag_set(c60, bond_idx_set, indene, 0)


    # test 2: run a loop and discover possible locations
    # this is the di-, bis- case
    if (False):
        write_all = True
        for i in range(1, len(c60.indene_pos)):
            bond_idx_set = [0, i]
            placed_indenes = c60.prepare_indene_set(bond_idx_set, indene)
            if placed_indenes == None:
                print "indene attachment FAILED for bond index set", bond_idx_set
            else:
                print "indene attachment SUCCEEDED for bond index set", bond_idx_set
                if (write_all):
                    decor = Structure(c60)
                    decor.append_fragments(placed_indenes)
                    decor.write_xyz("%s/%s.indene.%d.xyz" % (mols_dir, decor.fname, i))

    # test 3: tri-, tris-  case
    if (False):
        write_all = True
        for i in range(1, len(c60.indene_pos)):
            for j in range(i, len(c60.indene_pos)):
                bond_idx_set = [0, i, j]
                placed_indenes = c60.prepare_indene_set(bond_idx_set, indene)
                if placed_indenes == None:
                    print "indene attachment FAILED for bond index set", bond_idx_set
                else:
                    print "indene attachment SUCCEEDED for bond index set", bond_idx_set
                    if (write_all):
                        decor = Structure(c60)
                        decor.append_fragments(placed_indenes)
                        decor.write_xyz("%s/%s.indene.%d.%d.xyz" % (mols_dir, decor.fname, i, j))


    # test 4: tetra-,  kis-  case, getting silly!
    if (False):
        write_all = True
        for i in range(1, len(c60.indene_pos)):
            if (c60.prepare_indene_set([0,i], indene) != None):
                for j in range(i+1, len(c60.indene_pos)):
                    if (c60.prepare_indene_set([0,i,j], indene) != None):
                        for k in range(j+1, len(c60.indene_pos)):
                            bond_idx_set = [0, i, j, k]
                            placed_indenes = c60.prepare_indene_set(bond_idx_set, indene)                    
                            if placed_indenes == None:
                                print "indene attachment FAILED for bond index set", bond_idx_set
                            else:
                                print "indene attachment SUCCEEDED for bond index set", bond_idx_set
                                if (write_all):
                                    decor = Structure(c60)
                                    decor.append_fragments(placed_indenes)
                                    decor.write_xyz("%s/%s.indene.%d.%d.%d.xyz" % (mols_dir, decor.fname, i, j, k))

                                

def oxidizer_main():
    if (len(sys.argv) < 3):
        print "python atomer.py <structure xyz><number_of_oxygens> <number_of_realizations> "
#        print "python atomer.py <structure xyz><number_of_oxygens> <number_of_realizations> [nus]"
#        print "attach oxygen to PCBM. nus(==not under spinach) will prevent putting"
        print "attach oxygen to PCBM. "
#        print "oxygen between carbons where one of the carbons attaches to the spinach"
        sys.exit()

    M = int(sys.argv[2])
    N = int(sys.argv[3])
#    not_under_spinach = False
#    if (len(sys.argv) == 5):
#        if (sys.argv[4] == "nus"):
#            not_under_spinach = True

    # load structure
    pcbm = C60(fname=sys.argv[1])
    pcbm.find_neighbors()
    show_raw = True
    if (show_raw):
        print "Raw structure"
        for i in  range(0,len(pcbm.atoms)):
            atm = pcbm.atoms[i]
            print i, atm.elmnt, atm.pos, atm.nb
    # find pentgons and hexagons
    pcbm.find_shapes()
    # find 5-6 bonds
    pcbm.find_5_6_bonds()
    # find oxygen positions
    pcbm.find_oxy_pos()

    # write results
    pcbm.write_oxpos_subset()
    pcbm.write_xyz("pcbm_shifted.xyz")

    # potential use:  setup N randomly selected subsets of size M
    os.system("mkdir subsets 2>/dev/null")
    for i in range(0,N):
        idx = random.sample(range(0,len(pcbm.oxpos)), M)
        dirname = "pcbm_%dO_%dof%d" % (M, i+1, N)
        os.system("mkdir subsets/%s 2>/dev/null" % dirname)
        pcbm.write_oxpos_subset(idx, "subsets/%s/%s_ox_only.xyz" % (dirname,dirname))
        pcbm.write_gaussian_input(idx, i, M, N, dirname)


def indenizer_main():
    from optparse import OptionParser
    import os.path
    usage = "usage: %prog <c60_atoms_file> <indene_file> [options] "
    parser = OptionParser(usage=usage)    
    parser.add_option("-6", "--without_6-6", dest="no_6_6", help="do NOT connect at 6-6 bonds (False)", action="store_true", default=False)
    parser.add_option("-5", "--with_5-6", dest="do_5_6", help="connect at 5-6 bonds (False)", action="store_true", default=False)
    parser.add_option("-i", "--without_inversion", dest="no_inversion", help="do NOT consider both \'left\' and \'right\' indene orientation (False)", action="store_true", default=False)
    parser.add_option("-f", "--files", dest="mols_dir",  type="string", default="mols",
                                    help="where to dump output files")
    parser.add_option("-b", "--without_bis", dest="no_bis", help="do NOT do bis- case enumeration (False)", action="store_true", default=False)
    parser.add_option("-t", "--with_tri", dest="do_tris", help="do tris- case enumeration (False)", action="store_true", default=False)
    parser.add_option("-m", "--with_mono", dest="do_mono", help="do mone- [just attach one, good for debugging] (False)", action="store_true", default=False)
    parser.add_option("-q", "--with_equivalent", dest="with_equivalent", help="allow equivalent structures", action="store_true", default=False)
    parser.add_option("-a", "--accuracy", dest="accuracy", help="0-3, selects Gaussian basis set", type="int", default=1)
    parser.add_option("-n", "--nstates", dest="nstates", help="number of states Gaussion calcs.", type="int", default=12)

    
    (options, args) = parser.parse_args()
    # internal options versus command line options, slight difference
    opts = Options(do_6_6=not options.no_6_6, do_5_6=options.do_5_6, do_inversion=not options.no_inversion, 
                   do_bis=not options.no_bis, do_tris=options.do_tris, do_mono=options.do_mono, 
                   with_equivalent=options.with_equivalent, accuracy=options.accuracy, nstates=options.nstates)

    print options.no_6_6, options.do_5_6, options.no_inversion
    print opts.do_6_6, opts.do_5_6, opts.do_inversion
    print "required args: ", args
    if (len(args) != 2):
        print "not all require args supplied, see --help"
        sys.exit()
    mols_dir = options.mols_dir
    if (not os.path.isdir(mols_dir)):
        os.mkdir(mols_dir)
    #############

    c60 = C60(fname=sys.argv[1])
    c60.find_neighbors()

    # find pentgons and hexagons
    c60.find_shapes()
    # find bonds between 5 and 6 and between 6 and 6 member rings
    c60.find_5_6_bonds()
    c60.find_6_6_bonds()
    # find indene positions
    c60.find_indene_pos(opts.do_5_6, opts.do_6_6, opts.do_inversion)

    indene = Fragment(fname=sys.argv[2])
    indene.find_neighbors()
    indene.pin_idx = [6,8] # indices in indene.atoms of atoms that bond to c60
    fragments = [c60,indene]
    for f in fragments:
        print f.fname
#        f.write_xyz

    print "there are %d possible indene locations" % ( len(c60.indene_pos))

    # test 0: can we just add it once, no checking for compatibility
#    c60.add_indene(0, indene)
#    c60.add_indene(1, indene)

    if (False):
        enumeration_tests(c60, indene, opts, mols_dir)

    prepare_structs_and_dirs(c60, indene, opts, mols_dir)

def main():
    do_indene = True
    do_oxygen = False
    
    if (do_oxygen):
        oxidizer_main()
    if (do_indene):
        indenizer_main()


if __name__ == "__main__":
    main()

