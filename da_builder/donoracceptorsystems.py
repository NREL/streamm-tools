#!/usr/bin/env python

import os, sys, fileinput, random
from copy import deepcopy
from basics import *
from fragments import *


def get_connecting_elmnts(don, acc, dt, at):
    """
    Args:
        don (..type?..): donor element
        acc (..type?..): acceptor element
        dt             : index for donor element
        at   
    """

    don_elmnt = don.atoms[dt].elmnt
    acc_elmnt = acc.atoms[at].elmnt
    return (don_elmnt, acc_elmnt)


def is_pair(a,b, pair):
    """
    Args:
        a 
        b
        pair

    Returns:
        bool-->?
    """
    if (a == pair[0] and b == pair[1]) or (a == pair[1] and b == pair[0]):
        return True
    else:
        return False


def illegal_connection(don, acc, dt, at):
    """
    enforce "no N-O bonds", etc.
    Args: ?

    Returns: (bool)
    """
    don_elmnt, acc_elmnt = get_connecting_elmnts(don, acc, dt, at)
    if ((don_elmnt == "N" and acc_elmnt=="O") or (acc_elmnt == "N" and don_elmnt=="O")):
        print "ERROR: no N-O bonds allowed"
        return True
    elif (don_elmnt != "C" and acc_elmnt != "C"):
        print "ERROR: no hetero-hetero bonds allowed (all bonds must involve carbon)"
        return True
    else:
        return False

def get_bond_len(don, acc, dt, at):
    """
    Args: ?

    Returns: bond length
    """

    don_elmnt, acc_elmnt = get_connecting_elmnts(don, acc, dt, at)
    if (is_pair(don_elmnt, acc_elmnt, ["C","F"])):
        return 1.35
    elif (is_pair(don_elmnt, acc_elmnt, ["C","Si"])):
        return 1.94
    elif (is_pair(don_elmnt, acc_elmnt, ["C","H"])):
        return 1.09
    else:
        return 1.4


def connect_parts(don, acc, dt, dc, at, ac, rot_ctrl):
    """
    'don' and 'acc' are now just labels, each can be any fragment
    d,a for donor, acceptor, t,c for term, cap
    """
    #  print "dt, dc, at, ac", dt, dc, at, ac
    if illegal_connection(don, acc, dt, at):
        print "ERROR: illegal connection attempted"
        return None

    bond_len = get_bond_len(don, acc, dt, at)
    newfrag = Fragment(frag=don)
    # temp shift of h-bond to C-C bond length, then attachment
    #e.g. maybe 9 is donor H idx, 0 is acceptor C idx, 7 is donor C idx, 8 is acceptor H idx.
    hpos = newfrag.atoms[dc].pos
    newpos = vec_axpy(bond_len, vec_normalize(vec_sub(newfrag.atoms[dc].pos, newfrag.atoms[dt].pos)), newfrag.atoms[dt].pos)
    newfrag.atoms[dc].pos = newpos
    pacc = newfrag.place_fragment(acc, [[dc, at],[dt, ac]], rot_ctrl)  # connect newfragor term to acceptor cap, and vice versa
    newfrag.atoms[dc].pos = hpos
    # delete hydrogens at bond
    # if (newfrag.count_funccaps_for(dt) > 1 or newfrag.count_fcaps_for(dt) > 1):
    if (newfrag.count_funccaps_for(dt) > 1):
        # print "saveing newfrags link (funcs, funccaps)", dt, newfrag.count_funccaps_for(dt)
        pass
    else:
        # print "deleting newfrag link ", dt
        newfrag.delete_link(dt)

    pacc.delete_link(at)
    newfrag.delete_atom(dc)
    pacc.delete_atom(ac)
    newfrag.append_fragments([pacc])
    return newfrag


def connect_frags(don, acc, term_don, term_acc, rot_ctrl):
    """
    'don' and 'acc' are now just labels, each can be any fragment
     d,a for donor, acceptor, t,c for term, cap
     """
    # print "get term for donor ", term_don 
    dt,dc = don.get_term(term_don)
    # print "get term for acc ", term_acc
    at,ac = acc.get_term(term_acc)
    # whole = connect_parts(don, acc, dt, dc, at, ac, bond_len, "match_normal")
    # whole = connect_parts(don, acc, dt, dc, at, ac, bond_len, "stick_up")
    whole = connect_parts(don, acc, dt, dc, at, ac, rot_ctrl)
    return whole


def connect_rgroup(don, res, func_don, func_res, rot_ctrl):
    """
    'don' and 'res' are now just labels, each can be any fragment
     d,r for donor, residue, f,c for func, cap
     """
    #    print "get func for donor ", func_don 
    df,dc = don.get_func(func_don)
    #    print "get func for res ", func_res
    rf,rc = res.get_rgroup(func_res)
    #    whole = connect_parts(don, res, df, dc, rf, rc, bond_len, "match_normal")
    whole = connect_parts(don, res, df, dc, rf, rc, rot_ctrl)
    return whole


def connect_fgroup(don, res, func_don, func_res, rot_ctrl):
    """
    'don' and 'res' are now just labels, each can be any fragment
    d,r for donor, residue, f,c for func, cap
    """
    df,dc = don.get_f(func_don)
    rf,rc = res.get_rgroup(func_res)
    whole = connect_parts(don, res, df, dc, rf, rc, rot_ctrl)
    return whole


def tokenize(s):
    import re
    #  print "tokenizing ", s
    p1 = re.split("\(([^)]*)\)", s)
    # splits "D (R1 R2 ) A (R4)" into D, R1 R2, A, R4
    #    print "p1 = ", p1
    allp = []
    even = False
    for p in p1:
        if even:
            allp.append ('(')
        p2 = p.split()
        #        print "p = ", p
        #        print "p2 = ", p2
        for g in p2:
            allp.append(g)
        if even:
            allp.append(')')
        if even:
            even=False
        else:
            even=True
    return allp

def make_token_sets(s):
    """
    ['R1', 'R2', 'R1'] -> [[0,2], [1]]
    """
    found = []
    idx_sets = []
    for i in range(0,len(s)):
        targ = s[i]
        if (targ not in found):
            idx_sets.append([j for j in range(0,len(s)) if s[j] == targ])
        found.append(targ)
    return idx_sets


def illegal_rgroup_spec(subfragnames, curparent):
    """
    two things to check: 1) no "None" rgroup fragments 2) rgroups matching template of parent
    Note: order of these comparisons is important.  Want to check that we have enough subfrags before
    allowing NULL groups [TODO: rethink, bad logic].
    """
    if ("ERROR" in subfragnames):
        print "ERROR: no null rgroups allowed"
        return True
    elif (len(subfragnames) < len(curparent.rgroup_tokens)):
        print "ERROR: too few rgroups given"
        return True
    elif (len(subfragnames) > len(curparent.rgroup_tokens)):
        print "ERROR: too many rgroups given"
        return True
    
    #    print "checking that ", subfragnames
    #    print "matches ", curparent.rgroup_spec_str, " ie ", curparent.rgroup_tokens 
    #e.g., checking that  [1, 1]
    #matches  D(R1,R1)  ie  ['R1', 'R1']
    idx_sets = make_token_sets(curparent.rgroup_tokens)
    for i in range(0, len(idx_sets)):
        idx = idx_sets[i]
        if (idx != []):
            val = subfragnames[idx[0]]
            for j in range(0, len(idx)):
                # print "checking ", val , "   against ", subfragnames[idx[j]]
                 if (subfragnames[idx[j]] != val or ### basically we are just trying to match a pattern, not specific groups, so "R1 R1" matches "R2 R2"
                    (val.lower() != "null" and val != "0" and val[0].lower() != curparent.rgroup_tokens[idx[j]][0].lower())):  ### but we also don't let "f<X>" fill place meant for "R<X>"
                    print "ERROR: rgroup pattern %s does not match template of parent fragment %s" % (subfragnames[idx[j]], curparent.rgroup_tokens[idx[j]])
                    print "demanded by template %s in position %d" % (curparent.rgroup_spec_str, idx[j])
                    return True

#    print "checking for NULL in illegal place"
    for i in range(len(subfragnames)):
        if (subfragnames[i]=="NULL" and curparent.rgroup_tokens[i][0].lower() == "r"):
            print "no null rgroups allowed"
            return True

#    print "legal rgroup specification"

#    elif ("NULL" in subfragnames):
#        print "null rgroups allowed"


    return False

    

def make_frags(bblocks, parts, rot_ctrl):
    """
    this routine makes all the top level fragments, so it takes
    care of attaching rgroups.  It returns a  list of the top level fragments, ie
    a list of decorated (by R groups) donors and acceptors, and spacers
    """
    frags = []
    i = 0

    curparent = None
    curresdict = None
    while (i < len(parts)):
        c = parts[i]
        c = c.strip()
        reverse_attach = False
        
        if (c[0] == "r"):
            # reverse attachment order of this block
            c = c[1:]
            reverse_attach = True
#        print "outer: ",c
        if (c in bblocks.alldondict):
            frag_template = bblocks.alldondict[c]
            newfrag = Fragment(frag=frag_template)
            newfrag.type = "donor"
            curparent = frag_template
            curresdict = bblocks.alldonresdict if len(bblocks.alldonres) > 0 else bblocks.allresdict
            i += 1
            
            # make_ff TK
            newfrag.setall_resnumb(i)
            newfrag.setall_restype(c)

        elif (c in bblocks.allaccdict):
            frag_template = bblocks.allaccdict[c]
            newfrag = Fragment(frag=frag_template)
            newfrag.type = "acceptor"
            curparent = frag_template
            curresdict = bblocks.allaccresdict if len(bblocks.allaccres) > 0 else bblocks.allresdict
            i += 1
            
            # make_ff TK
            newfrag.setall_resnumb(i)
            newfrag.setall_restype(c)

        elif (c in bblocks.allspacerdict):
            frag_template = bblocks.allspacerdict[c]
            newfrag = Fragment(frag=frag_template)
            newfrag.type = "spacer"
#            curparent = None # not allowed to attach rgroups to spacer
            curparent = frag_template # now we allow to attach rgroups to spacer
            curresdict = bblocks.allspacerresdict if len(bblocks.allspacerres) > 0 else bblocks.allresdict
            i += 1
            
            # make_ff TK
            newfrag.setall_resnumb(i)
            newfrag.setall_restype(c)

        elif (c in bblocks.alltermdict):
            frag_template = bblocks.alltermdict[c]
            newfrag = Fragment(frag=frag_template)
            newfrag.type = "terminal"
            curresdict = bblocks.alltermresdict if len(bblocks.alltermres) > 0 else bblocks.allresdict
            curparent = frag_template
            i += 1
            
            # make_ff TK
            newfrag.setall_resnumb(i)
            newfrag.setall_restype(c)

        elif c=='(':
            i = i+1
            c = parts[i]
            subfrags = []
            subfragnames = []
            while (c != ')'):
                c=parts[i]
#                print "inner: ", c
                if (c == ")"):
                    pass
#                elif (c in bblocks.allresdict):
                elif (c != "0" and c.lower() != "null" and c in curresdict):
#                    frag_template = bblocks.allresdict[c]                    
                    frag_template = curresdict[c]
                    afrag = Fragment(frag=frag_template)
                    afrag.type = "substituent"
                    if (c[0].lower() == "f"):
                        afrag.subtype = "fgroup"
                    else:
                        afrag.subtype = "rgroup"
                    subfrags.append(afrag)
                    subfragnames.append(c)
                    
            
                    # make_ff TK
                    afrag.setall_resnumb(i)
                    afrag.setall_restype(c)

                elif (c == "0" or c.lower() == "null"):
                    subfrags.append(None)
                    subfragnames.append("NULL")
                else:
                    print "ERROR: rgroup %s does not exist in resdict:" % c, curresdict
                    return None
#                else:
#                    subfrags.append(None)
#                    subfragnames.append("ERROR")  # this will now be illegal
                i = i+1
            rtarget_idx = 0
            ftarget_idx = 0
            # target_idx on donor/acceptor gets "used up" by connecting an rgroup,
            # so we only increment it if user put "null" (anything but R) as an "argument"
            # to the relevant D or A tag
            if (len(subfrags) > newfrag.count_funccaps()):
                print "ERROR: more rgroups (%d) than functional locations on donor acceptor based on %s" % (len(subfrags), newfrag.fname)
                return None

            # rgroups need to stick up if possible
            local_rot_ctrl = deepcopy(rot_ctrl)
            local_rot_ctrl.free_rot_criteria = "stick_up"

            # wasteful (but all this is fast) but easier to check in one spot whether specified r-group
            # arrangement was legal
#            print "checking rgroup:", subfragnames, curparent
            if (illegal_rgroup_spec(subfragnames, curparent)):
                print "ERROR: illegal rgroup specification"
                return None


            for j in range(0,len(subfrags)):
                if subfrags[j] != None:
                    if (subfrags[j].subtype == "fgroup"):
                        afrag = connect_fgroup(newfrag, subfrags[j], ftarget_idx, 0, local_rot_ctrl)
                    else:
                        afrag = connect_rgroup(newfrag, subfrags[j], rtarget_idx, 0, local_rot_ctrl)
                    if (afrag == None):
                        print "ERROR: rgroup attachment failure"
                        return None
                    afrag.type = newfrag.type
                    newfrag = afrag
                else:
                    ftarget_idx += 1  ### this assumes all NULLs are in f-group locations, ie 
                                      ### no NULL rgroups.


        elif (c==')'):
            newfrag = None
            i += 1
        else:
            print "ERROR: fragment %s does not exist" % c
            return None


        if (newfrag != None and (i >= len(parts)-1 or parts[i] != "(")):
#            print "about to append a new frag, c,i=", c, i
            newfrag.reverse_attach = reverse_attach
            frags.append(newfrag)
#    print frags
    return frags


def build_from_str(bblocks, s, options):
    """
    eg s = "D1 D3 A4(R1 R2) D D(R4 R8) A"
    """
    parts = tokenize(s)

    if (options.alternation_pattern != None):
        criteria = "follow_pattern"
    elif (options.make_trans):
        criteria = "alternate"
    else:
        criteria = "max_dist_to_frag"
    rot_ctrl = RotationControl(criteria, options.alternation_pattern)

    # build list of decorated (rgroups applied) top level fragments
    frags = make_frags(bblocks, parts, rot_ctrl)
    if (frags == None):
        print "Structure generation FAILED for string ", s
        return None
    for newfrag in frags:
        if (newfrag == None):
            print "Structure generation FAILED for string ", s
            return None

    # now join them
    fragidx = 0
    rot_ctrl.counter = 0
    for newfrag in frags:
        
        # make_ff TK
        #newfrag.setall_resnumb(fragidx+1)
        #newfrag.setall_restype(s)
        
        print "processing fragment of type : ", newfrag.type," fragidx ",fragidx," s ",s
        #  print "processing fragment of type : ", newfrag.type, " rot ctrl crit= ", rot_ctrl.free_rot_criteria
        if (fragidx == 0):
            frag = Fragment(frag=newfrag)
            frag.type = newfrag.type
            frag.reverse_attach = newfrag.reverse_attach
            firstfrag = frag
            #  print "frag idx", fragidx, "frag.type ", frag.type, "reverse_attach ", frag.reverse_attach 
        elif (newfrag != None):
            # by default, 1st (idx 0) link point on new frag gets plugge in to 
            # 2nd link point on existing (old) fragment, ie, "end-to-end" 
            oldfrag_conn_idx = 1
            newfrag_conn_idx = 0
            # then we adjust for special cases:
            if (newfrag.reverse_attach):
                newfrag_conn_idx = 1  # connect to second instead of first
            if (firstfrag.reverse_attach and fragidx == 1):
                oldfrag_conn_idx = 0  # connect to first instead of second, where initial fragment was reversed.
            if (frag.type == "terminal"):
                oldfrag_conn_idx = 0  # only has one, will crash if not 0
            if (newfrag.type=="terminal"):
                newfrag_conn_idx = 0  # only has the one.

#            print "frag idx", fragidx, "frag.type newfrag.type ", frag.type, newfrag.type, "reverse_attach ", newfrag.reverse_attach, "attachment indices: ", oldfrag_conn_idx, newfrag_conn_idx
            bigfrag = connect_frags(frag, newfrag, oldfrag_conn_idx, newfrag_conn_idx, rot_ctrl)   # connect newfrag to frag in first available way (0,0)
            bigfrag.type = frag.type
            frag = bigfrag
            # this was changed from connect_frags(frag, newfrag, 1, 0) 10/11/11, PG, installing "terminals" (only have 1 term, so idx=1 crashed)
            rot_ctrl.counter += 1  # rotation counter increments for each bond, which only start at second fragment 
        fragidx += 1

    if (frag != None):
        if (frag.atoms_too_close()):
            print "WARNING: legal structure, but atoms too close for string ", s

        print "Structure generation SUCCEEDED for string ", s
    return frag


def collapse(str):
    """
    """
    p1 = str.split()
    p2 = []
    for p in p1:
        p2.append(p.strip())
    
    s = ""
    for p in p2:
        for c in p:
            if (c == "("):
                c = "_"
            elif (c == ")"):
                c = "_"
            s+=c
    return s



def first_tests(don, acc):
    """
    Driver for write_xyz ?
    """
    new1 = connect_frags(don, acc, 1, 1)
    print "donor + acceptor"
    new1.dump()
    new1.write_xyz("DA.xyz")

    new2 = connect_frags(new1, don, 1, 0)
    print "donor + acceptor + donor"
    new2.dump()
    new2.write_xyz("DAD.xyz")

    new3 = connect_frags(new2, acc, 1, 1)
    new3.write_xyz("DADA.xyz")


def load_one_dir(dondir, has_hdr, key, subsets, allownokey):
    """
    """
    import re,os
    files = os.listdir(dondir)
    frags = []
    names = []
    asdict = {}
    for f in files:
        test2 = re.match("([^.].*)\.cply$",f)
        if (test2 != None):
            fbase = test2.group(1)
#            print "testing ", fbase, key, subsets, allownokey
            if ((key not in subsets and allownokey) or (key in subsets and fbase in subsets[key])):  # if subsets for this bb type specified, make sure it's one we want to load
                fname = "%s/%s" % (dondir,f)
#                print "loading ", fname,  "  fbase= ", fbase, "for key ", key
                if (has_hdr):
                    frag = FragmentWithHeader(fname)
                else:
                    frag = Fragment(fname)
                frags.append(frag)
                names.append(fbase)
                asdict[fbase] = frag

    print "key is:", key

    # need special case for null spacer:
    if (key == "spacers"):
        if ("spacers" not in subsets or "0" in subsets['spacers']):
            frags.append(None)
            names.append("")
            print "added null spacer"
            asdict["0"] = None
        
    # more generally, null substituents now allowed
    if (key in ['donor_functional_groups', 'acceptor_functional_groups','terminal_functional_groups','spacer_functional_groups']):
        if (key not in subsets or "0" in subsets[key]):
            frags.append(None)
            names.append("0")
            print "added null %s" % key
            asdict["0"] = None

    
    return frags, names, asdict


def load_cfg_file(fname):
    """
    """
    import fileinput
    d = {}
    if (fname == None):
        return d
    for ln in fileinput.input(fname):
        ln = ln.strip()
        if (ln == "" or ln[0] == "#"):
            pass
        else:
            ln = ln.split('=')
            if (len(ln) < 2):
                print "syntax error in config file at",ln
            key = ln[0].strip()
            val = ln[1].strip()
            val = val.split("#")[0].strip()
            d[key] = val
    return d


def load_subsets_dict(fname):
    """
    """
    d = load_cfg_file(fname)
    dd = {}
    if ("donors" in d):
        dd['donors'] = d['donors'].split()
    if ("acceptors" in d):
        dd['acceptors'] = d['acceptors'].split()
    if ("donors2" in d):
        dd['donors2'] = d['donors2'].split()
    if ("acceptors2" in d):
        dd['acceptors2'] = d['acceptors2'].split()
    if ("fdonors" in d):
        dd['fdonors'] = d['fdonors'].split()
    if ("facceptors" in d):
        dd['facceptors'] = d['facceptors'].split()
    if ("functional_groups" in d):
        dd['functional_groups'] = d['functional_groups'].split()
    if ("spacers" in d):
        dd['spacers'] = d['spacers'].split()
    if ("terminals" in d):
        dd['terminals'] = d['terminals'].split()
    if ("donor_functional_groups" in d):
        dd['donor_functional_groups'] = d['donor_functional_groups'].split()
    if ("acceptor_functional_groups" in d):
        dd['acceptor_functional_groups'] = d['acceptor_functional_groups'].split()
    if ("terminal_functional_groups" in d):
        dd['terminal_functional_groups'] = d['terminal_functional_groups'].split()
    if ("spacer_functional_groups" in d):
        dd['spacer_functional_groups'] = d['spacer_functional_groups'].split()
    return dd



class BuildingBlocks:
    """
    
    """
    def __init__(self, bblocks_dir, subsets_file):
        """
        Constructor
        """

        self.alldon = None
        self.allacc = None
        self.allres = None
        self.allterm = None
        self.allspacer = None
        self.bblocks_dir = bblocks_dir
        subsets_dict = load_subsets_dict(subsets_file) # this is outside class?
        self.load_building_blocks(self.bblocks_dir, subsets_dict)

    def load_building_blocks(self, bblocks_dir, subsets_dict):
        """
        """   
        import os,re
        bbdir = bblocks_dir
        dondir = "%s/donors" % (bbdir)
        accdir = "%s/acceptors" % (bbdir)
        termdir = "%s/terminals" % (bbdir)
        resdir = "%s/functional_groups" % (bbdir)
        spacerdir = "%s/spacers" % (bbdir)
        self.alldon, self.alldonname, self.alldondict = load_one_dir(dondir, True, "donors", subsets_dict, True)
        self.allacc, self.allaccname, self.allaccdict = load_one_dir(accdir, True, "acceptors", subsets_dict, True)
        self.alldon2, self.alldonname2, self.alldondict2 = load_one_dir(dondir, True, "donors2", subsets_dict, False)
        self.allacc2, self.allaccname2, self.allaccdict2 = load_one_dir(accdir, True, "acceptors2", subsets_dict, False)
        self.alldonf, self.alldonnamef, self.alldondictf = load_one_dir(dondir, True, "fdonors", subsets_dict, False)
        self.allaccf, self.allaccnamef, self.allaccdictf = load_one_dir(accdir, True, "facceptors", subsets_dict, False)
        self.allres, self.allresname, self.allresdict = load_one_dir(resdir, False, "functional_groups", subsets_dict, True)
        self.alldonres, self.alldonresname, self.alldonresdict = load_one_dir(resdir, False, "donor_functional_groups", subsets_dict, True)
        self.allaccres, self.allaccresname, self.allaccresdict = load_one_dir(resdir, False, "acceptor_functional_groups", subsets_dict, True)
        self.alltermres, self.alltermresname, self.alltermresdict = load_one_dir(resdir, False, "terminal_functional_groups", subsets_dict, True)
        self.allspacerres, self.allspacerresname, self.allspacerresdict = load_one_dir(resdir, False, "spacer_functional_groups", subsets_dict, True)
        self.allterm, self.alltermname, self.alltermdict = load_one_dir(termdir, True, "terminals", subsets_dict, True)
        self.allspacer, self.allspacername, self.allspacerdict = load_one_dir(spacerdir, True, "spacers", subsets_dict, True)
        #
        # Add molecules to donors
        #    this is done so molecule exist in a different directory so they will not be called by enumerate.py
        #    but can be found when running donoracceptorsystems.py
        # 
        # TWK TRAVIS
        # 
        moldir = "%s/fullerene" % (bbdir)
        self.allmol, self.allmolname, self.allmoldict = load_one_dir(moldir, True, "mols", subsets_dict, True)
        for key,val in self.allmoldict.iteritems():
            if key not in self.alldondict:
                self.alldondict[key] = val
        
        for key,val in self.alldondict2.iteritems():
            if key not in self.alldondict:
                self.alldondict[key] = val
        for key,val in self.allaccdict2.iteritems():
            if key not in self.allaccdict:
                self.allaccdict[key] = val
        for key,val in self.alldondictf.iteritems():
            if key not in self.alldondict:
                self.alldondict[key] = val
        for key,val in self.allaccdictf.iteritems():
            if key not in self.allaccdict:
                self.allaccdict[key] = val
        for key,val in self.alldonresdict.iteritems():
            if key not in self.allresdict:
                self.allresdict[key] = val
        for key,val in self.allaccresdict.iteritems():
            if key not in self.allresdict:
                self.allresdict[key] = val
        for key,val in self.alltermresdict.iteritems():
            if key not in self.allresdict:
                self.allresdict[key] = val

        print "Loaded building blocks:"
        print "donors", self.alldonname
        print "acceptors", self.allaccname
        print "donors2", self.alldonname2
        print "acceptors2", self.allaccname2
        print "fdonors", self.alldonnamef
        print "facceptors", self.allaccnamef
        print "terminals", self.alltermname
        print "spacers", self.allspacername
        print "functional_groups", self.allresname
        print "donor_functional_groups", self.alldonresname
        print "acceptor_functional_groups", self.allaccresname
        print "terminal_functional_groups", self.alltermresname
        print "spacer_functional_groups", self.allspacerresname








def build_usages(input_string, bblocks):
    """
    extract names of which donors, acceptors, etc were used
    could be done at construction time, but bookkeeping sort of a pain right now, 
    and easily partially recreated here

    Args:
        input_string (?):?
        bblocks      (?):?
    Returns:  don, acc, dres, ares, tres, sres, spacer, term
    """
    don = set([])
    acc = set([])
    dres = set([])
    ares = set([])
    tres = set([])
    sres = set([])
    term = set([])
    spacer = set([])
    parts = tokenize(input_string)
    for c in parts:
        if c[0] == "r":
            c = c[1:]
        if (c in bblocks.alldondict):
            don.add(c)
            inpart = "don"
        elif (c in bblocks.allaccdict):
            acc.add(c)
            inpart = "acc"
        elif (c in bblocks.alltermdict):
            term.add(c)
            inpart = "term"
        elif (c in bblocks.allspacerdict):
            spacer.add(c)
            inpart = "spacer"
        elif (inpart == None):
            print "ERROR: first part is not D, A, or T!, expect crash"
        elif (c in bblocks.allresdict):
            if (inpart=="don"):
                dres.add(c)
            elif (inpart=="acc"):
                ares.add(c)
            elif (inpart=="term"):
                tres.add(c)
            elif (inpart=="spacer"):
                sres.add(c)
    return don, acc, dres, ares, tres, sres, spacer, term


def build_meta(frag, short_name, input_string, bblocks, options):
    """
    Args:?

    Returns: dd ?
    """

    import time, datetime
    d = {}
    dd = {}
    dd["metadata"] = d 
    d['tag'] = short_name
    d['basis'] = get_basis_str(options.accuracy)
    dt = datetime.datetime.fromtimestamp(time.time())
    d['date_time'] = dt.isoformat()
    don, acc, dres, ares, tres, sres, spacer, term = build_usages(input_string, bblocks)
    d['donors'] = list(don)
    d['acceptors'] = list(acc)
    d['donor_substituents'] = list(dres)
    d['acceptor_substituents'] = list(ares)
    d['terminal_substituents'] = list(tres)
    d['spacer_substituents'] = list(sres)
    d['terminals'] = list(term)
    d['spacers'] = list(spacer)
    d['nstates'] = options.nstates
    d['accuracy'] = options.accuracy
    return dd



def gen_struct(base_input_str, bblocks, options, number, write_files = True):
    """
    Args:

    Returns: bool ?
    """
    from time import strftime, gmtime

    input_str = ""
    
    for i in range (0,number):
        input_str += base_input_str 
        input_str += " "

        print "concatonating  ",base_input_str," to input_str string "

    frag = build_from_str(bblocks, input_str, options)


    # Set REPOLOCATION in pbs template script
    repoPath=options.repoPath

    # If repoPath not set from options set to default /projects/opv/PUBLISH-'user'-datestamen
    if (repoPath == None):
        userName = os.getenv("USER")
        timeStamp = strftime("%Y_%b_%d_%a_%H_%M", gmtime())
        repoPath="/projects/opv/PUBLISH-" + userName + "-" + timeStamp

    # Remote template file directories checks
    opvPath=options.opvPath
    if ( not os.path.exists(opvPath) ):
        print "OPV pbs/com template files not found. Check opv-project repo location"
        sys.exit(0)


    if (frag != None and write_files):
        short_name = collapse(base_input_str)
        struct_dir = "%s/%s" % (options.mols_dir, short_name)
        if (not os.path.isdir(struct_dir)):
            os.mkdir(struct_dir)
        xyz_name = "%s/acc%d_%s_n%d.xyz" % (struct_dir, options.accuracy, short_name, number)
        job_name = "acc%d_%s_n%d" % (options.accuracy, short_name, number)
        frag.write_xyz(xyz_name)

        print "------------------------------------------------------------------------------------------------------"

        fullFilePath=os.path.join(opvPath,"donoracceptor.com.template")
        if ( not os.path.exists(fullFilePath) ):
            print "Template file", fullFilePath, " not found. Check --opvPath option"
            #sys.exit(0)
        else:
            frag.write_com(fullFilePath, xyz_name, job_name, get_basis_str(options.accuracy), options.nstates)

        fullFilePath=os.path.join(opvPath,"donoracceptor.pbs.template")
        if ( not os.path.exists(fullFilePath) ):
            print "Template file", fullFilePath, " not found. Check --opvPath option"
            #sys.exit(0)
        else:
            frag.write_pbs(fullFilePath, xyz_name, job_name, repoPath)


        #
        # SWS: adding new file templates explicitly for restart com files
        #
        fullFilePath=os.path.join(opvPath,"donoracceptor.com.template.r1")
        if ( not os.path.exists(fullFilePath) ):
            print "Template file", fullFilePath, " not found. Check --opvPath option"
            # sys.exit(0)
        else:
            frag.write_com_restart(fullFilePath,  xyz_name, job_name, get_basis_str(options.accuracy), options.nstates)

        fullFilePath=os.path.join(opvPath,"donoracceptor.com.template.r2")
        if ( not os.path.exists(fullFilePath) ):
            print "Template file", fullFilePath, " not found. Check --opvPath option"
            #sys.exit(0)
        else:
            frag.write_com_restart(fullFilePath,  xyz_name, job_name, get_basis_str(options.accuracy), options.nstates)

        print "------------------------------------------------------------------------------------------------------"

        # Build meta data 
        json_data = build_meta(frag, short_name, base_input_str, bblocks, options)
        json_data['metadata']['number'] = number
        json_data['metadata']['n'] = number

        
        frag.write_meta(json_data, xyz_name)

    if (frag != None):
        return True
    else:
        return False


def gen_monomer(input_str, bblocks, options, write_files = True): 
    """
    Driver gen_struct
    """
    return gen_struct(input_str, bblocks, options, 1, write_files)


def gen_dimer(input_str, bblocks, options):
    """
    Driver gen_struct
    """
    return (gen_struct(input_str,  bblocks, options, 2))


def get_options():
    """
    Input options (mirrored in enumerate?)
    """
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog <structure defining string (e.g. \'D1 (R2 R2) A2 (R1)\')> [options] "
    parser = OptionParser(usage=usage)    
    parser.add_option("-f", "--files", dest="mols_dir",  type="string", default="mols",
                                    help="where to dump output files")
    parser.add_option("-b", "--bblocks", dest="bblocks_dir",  type="string", default="BuildingBlocks",
                                    help="where to look for the building blocks")
    parser.add_option("-a", "--accuracy", dest="accuracy", help="0=low, 1=medium, 2=high, 3=super-high, selects Gaussian basis set", type="int", default=1)
    parser.add_option("-n", "--nstates", dest="nstates", help="number of electronic states to calculate", type="int", default=12)
    parser.add_option("-d", "--dimer", dest="do_dimer", help="generate input for dimers", action="store_true", default=False)
    parser.add_option("-t", "--trimer", dest="do_trimer", help="generate input for trimer", action="store_true", default=False)
    parser.add_option("-q", "--tetramer", dest="do_tetramer", help="generate input for tetramers", action="store_true", default=False)
    parser.add_option("-s", "--subsets", dest="subsets_file",  type="string", default=None,
                                    help="file where subsets of building blocks to use are specified (e.g. \"donors = D1 D2\")")
    parser.add_option("-r", "--repeat", dest="repeat_units", help="generate up to this number of repeat units", type="int", default=1)
    parser.add_option("-z", "--trans", dest="make_trans", help="generate \'trans\' (alternating) structure", action="store_true", default=False)
    parser.add_option("-p", "--pattern", dest="alternation_pattern",  type="string", default=None, help="fragment rotation pattern (list of angles)")

    # Force field options 
    parser.add_option("--make_ff", dest="make_ff", help="Generate force field input files from a finished calculation  ", action="store_true", default=False)
    parser.add_option("--norm_dihparam", dest="norm_dihparam",default=False,action="store_true",help="Normalize dihedral potential terms if single dihedral is specified in itp file  ")
    parser.add_option("--ff_charges", dest="ff_charges", help="Use force field charges rather than charges in cply file  ", action="store_true", default=False)

    # these are really for enumerator, but they share options
    parser.add_option("-c", "--class", dest="struct_class",  type="string", default="DA", help="space separated string of classes of donor acceptor systems to enumerate, among DA ADA DAD DD AA DAD DDA, e.g. \"DA DAD\"")
    parser.add_option("-e", "--equal_blocks_ok", dest="equal_blocks_ok", help="generate DiADi cases, etc", action="store_true", default=False)

    # this path is where results will be stored (eg /projects/opv/...)
    parser.add_option("--repoPath", dest="repoPath",  type="string", default=None,
                                    help="location of the REPOLOCATION to set in pbs template file. Default is /projects/opv/PUBLISH-$USER-'datestamp'")

    # this path points to the top-level OPV-repo to get pbs/com template files
    parser.add_option("--opvPath", dest="opvPath",  type="string", default=".",
                                        help="location of opv-project repo where com/pbs templates are located. Default is current directory")

    (options, args) = parser.parse_args()
    print "accuracy = ", options.accuracy
    print "nstates = ", options.nstates
    print "molecules base output directory = ", options.mols_dir
    print "using building blocks in ", options.bblocks_dir
    # options.subsets_file = "subset.txt"
    if (options.subsets_file != None):
        print "Using subsets of building blocks specified in ", options.subsets_file
    else:
        print "Not using subsets because options.subsets_file = ", options.subsets_file

    if (not os.path.isdir(options.mols_dir)):
        os.mkdir(options.mols_dir)


    # Track dates when mols output directories created
    mols_fulldir = os.path.abspath(options.mols_dir)
    molsTrackFile=os.path.join(options.opvPath, ".molsTrackFile")

    # Get time info for repo version tracking
    import datetime as dt
    tobj = dt.datetime.today()

    year   = tobj.year
    month  = tobj.month
    day    = tobj.day
    hour   = tobj.hour
    second = tobj.second
    dateStamp = str(year) + "-" + str(month) + "-" + str(day)

    if ( os.path.exists(molsTrackFile) ):
        fobj = open(molsTrackFile, 'a')
    else:
        fobj = open(molsTrackFile, 'w')

    repoTrackString = "Tracking " + mols_fulldir + "\n"
    fobj.write(repoTrackString)
    fobj.close()

    molsTimeStampFile = os.path.join(mols_fulldir, 'timestamp.dat')

    print "molsTrackFile     = ", molsTrackFile
    print "mols_fulldir      = ", mols_fulldir
    print "molsTimeStampFile = ", molsTimeStampFile

    if not os.path.exists(molsTimeStampFile):
        fobj = open(molsTimeStampFile, 'w')
        repoTrackString = "This directory created_on " + dateStamp + "\n"
        fobj.write(repoTrackString)
        fobj.close()

    return options, args


    
def main():
    options, args = get_options()
    input_str = args[0]
    print "donoracceptorsystems: setting up files for input string = ", input_str
    
    # load the building blocks
    bblocks = BuildingBlocks(options.bblocks_dir, options.subsets_file)

    # generate the structures and input files
    success = gen_monomer(input_str, bblocks, options)
    if (success == True and options.do_dimer):
        gen_dimer(input_str, bblocks, options)
    if (success == True and options.do_trimer):
        gen_struct(input_str,  bblocks, options, 3)
    if (success == True and options.do_tetramer):
        gen_struct(input_str,  bblocks, options, 4)

    if (success == True and options.repeat_units != 1):
        for irpt in range(1,options.repeat_units+1):
            gen_struct(input_str,  bblocks, options, irpt)

if __name__=="__main__":
    main()

