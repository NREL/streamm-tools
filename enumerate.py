from basics import *
from fragments import *
from donoracceptorsystems import *

import string

class MultiIndex(object):
    # multi-index class, indexing goes left to right (as opposed to counting in 
    # base "maxval", which would go right to left).
    def __init__(self,length, maxvals):
        self.maxvals = maxvals
        self.length = length
        self.midx = [0 for i in range(length)]
        
    def incr(self):
        done = True
        for i in range(self.length):
            if self.midx[i] < self.maxvals[i] - 1:
                self.midx[i] += 1
                done = False
                break;
            else:
                self.midx[i] = 0
        return done

        
class BuildingBlockEnumerator(BuildingBlocks):
    def __init__(self, bblocks_dir, subsets_file):
        BuildingBlocks.__init__(self, bblocks_dir, subsets_file)
        self.reset_counters()
    
    def reset_counters(self):
        self.tries = 0
        self.successes = 0
        self.good_strs = []

    def try_str(self, input_str, options, real_deal):
        # generate the structures and input files
        success = gen_monomer(input_str, self, options, write_files = real_deal)
        if (success == True and options.do_dimer and real_deal):
            gen_dimer(input_str, self, options)
        if (success == True and options.do_trimer and real_deal):
            gen_struct(input_str,  self, options, 3)
        if (success == True and options.do_tetramer and real_deal):
            gen_struct(input_str,  self, options, 4)
        if (success):
            self.successes += 1
            self.good_strs.append(input_str)
        self.tries += 1
        return success


    def get_resname(self, typename):
        from copy import copy
        if (typename=="donor"):
            specialres = self.alldonresname
        elif (typename=="acceptor"):
            specialres = self.allaccresname
        elif (typename=="terminal"):
            specialres = self.alltermresname
        elif (typename=="spacer"):
            specialres =  self.allspacerresname
        else:
            specialres = []

        if len(specialres) > 0:
            resname = copy(specialres)
        else:
            resname = copy(self.allresname)
            resname.append("0")  # we allow NULL substituent; ie no subst. on given location.
#        resname = [r for r in resname]  # w/out this, python is adding to the original list; this just makes it a copy
        return resname

    def try_set(self, frags, fragnames, typename, options, write_inputs, fluorinate):
        good_frag = []
        resname = self.get_resname(typename)

        print "Trying set ", fragnames, " of type ", typename, "implying residues " , resname

        for ifrag in range(0,len(frags)):
            frag = frags[ifrag]
            if (frag == None):  ## this is for null spacer, which is in "frags" if we're calling this fn. for the spacers
                continue
            print "fragment %d (name: %s) signature : " % (ifrag, fragnames[ifrag]), frag.rgroup_tokens
            print frag.unique, frag.template            
            midx = MultiIndex(len(frag.unique), [len(resname) for im in range(len(frag.unique))])
            done = False
            while not done:
#                print midx.midx
                # midx enumerates the possible substituents for this frag, taking into account
                # the signature, ie the fact that some subst's need to occur in appropriate pairs, etc.
                rvec = [resname[midx.midx[i]] for i in range(midx.length)]
                # check lead characters of rvec against those in frag.unique; they must match, so we don't mix, e.g.
                # fluorine substs. and regular substs.  Whether we allow NULL depends on whether its a "regular" rgroup
                # or a fluorination.  TODO: Ask Ross, is this the right idea?
#                print "trying rvec ", rvec
                goodtry = True
                for ilead in range(len(rvec)):                                                               ### bad ==   
                    rval = rvec[ilead]
                    if ((rval != "0" and (rval[0].lower() != frag.unique[ilead][0].lower())) or  ## not null and not a match to template
     #                   (rval == "0" and frag.unique[ilead][0].lower() == "r") or                       ## no nulls allowed as rgroups
                        (not fluorinate and (rval[0].lower() == "f" or 
                                             (len(rval) > 1 and rval[1].lower() == "f")))):## trying to fluorinate, but we're not doing that now.
                        goodtry = False
#                        print "fail: ", rvec, fluorinate
                        break
#                print ["f" in r.lower() for r in rvec]
                if (fluorinate and sum(["f" in r.lower() for r in rvec]) == 0):
                    goodtry = False  # delete case of all null f's when fluorinating; this case covered when we don't fluorinate
                    print "deleting allnullF", rvec
                if (goodtry):
                    s = "%s ( " % (fragnames[ifrag])
                    for i in range(len(frag.template)):
                        s += "%s " % (rvec[frag.template[i]])
                    s += ")"

                    print "try: ", s
                    if (self.try_str(s, options, write_inputs)):
                        good_frag.append(s)
#                print "curmidx:", midx.midx
                done = midx.incr()
            
        return good_frag


    def try_class(self, stuff, idxmap, spacers, options):
        # a bit tricky!
        ## eg stuff = [donors, acceptors], idxmap = [0,1,0] means we're making DAD, (where 1st and 2nd D's have to be the same)
        # stuff = [don, acc, don2], idxmap = [0,1,2] means we're making DAD2 (the D and D2 representatives need not be the same)
        midx = MultiIndex(len(stuff), [len(stuff[i]) for i in range(len(stuff))])
        empty_stuff = sum([len(bblist)==0 for bblist in stuff]) > 0 # true if one of the building block lists in stuff is empty.
        done = empty_stuff
        while (not done):
            for ispc in range(len(spacers)):
                s = ""
                for i in range(len(idxmap)):
                    # print i, midx.midx, idxmap, stuff
                    s += "%s " % stuff[idxmap[i]][midx.midx[idxmap[i]]]  # midx[i] is where we are in the enumeration of idxmap[i]th element of stuff.
                    s += "%s " % spacers[ispc]

                print "TRYING %s" % s
                skipBool=False

                #####################################
                # Selecting out strings (takes s)

                # skipBool = self.exclude_string(s)
                # s = self.dropSpacer(s)
                #####################################
  
                if (not skipBool):
                    self.try_str(s, options, True)

            done = midx.incr()


    #########################################################################
    #
    # Selecting out strings (takes s) with specific rules (hard-wired)
    # Arg 1 -- s is the string that try_str passes to donoracceptorsystems
    # Pass back flag to skip/execute
    #
    # NOTE:  Needs D61 at beginning and not copied in 2nd donor spot
    #########################################################################
    def exclude_string(self, eStr):

        splitStr = eStr.split(' ')
        strDlist = []
        for frag in splitStr:
            if 'D' in frag:
                strDlist.append(frag)
                        
        # Needs D61 at beginning and not copied in 2nd donor spot
        if ( (strDlist[0] != "D61") or (strDlist[0] == strDlist[1]) ):
            print "frag = ", strDlist, " is copied... skipping generation"
            skipFlag = True
        else:
            skipFlag = False

        return skipFlag


    #########################################################################
    #
    # Takes a donoracceptorsystems string and drops an element 'by hand'
    # NOTE:  Changing element to drop needs editing manually
    #########################################################################
    def dropSpacer(self, eStr):

        dropString = " Sp2 ( 0 0 ) "
        slen = -1*len(dropString)
        newString = eStr[:slen]
        print "dropString = ", dropString
        print "slen = ", slen
        print "newString = ", newString
        return newString


    def enum_smarter(self, options):
        only_unequal_bb = not options.equal_blocks_ok 
        self.reset_counters()

        classes = options.struct_class.split()        
        save_decorated_bb = False
        if ("BB" in classes):
            save_decorated_bb = True
            print "Generating input for lone donor and acceptor calcs"

        print "Before self.allspacername = ", self.allspacername

        good_don = self.try_set(self.alldon, self.alldonname, "donor", options, write_inputs=save_decorated_bb, fluorinate = False)
        good_acc = self.try_set(self.allacc, self.allaccname, "acceptor", options, write_inputs=save_decorated_bb, fluorinate = False)
        good_term = self.try_set(self.allterm, self.alltermname, "terminal", options, write_inputs=save_decorated_bb, fluorinate = False)
        good_spacer = self.try_set(self.allspacer, self.allspacername, "spacer", options, write_inputs=save_decorated_bb, fluorinate = False)

        ## fluorinated building blocks
        if (self.alldonf != []):
            good_donf = self.try_set(self.alldonf, self.alldonnamef, "donor", options, write_inputs=False, fluorinate = True)
        else:
            good_donf = []
        if (self.allaccf != []):
            good_accf = self.try_set(self.allaccf, self.allaccnamef, "acceptor", options, write_inputs=False, fluorinate = True)
        else:
            good_accf = []
        good_spacerf = self.try_set(self.allspacer, self.allspacername, "spacer", options, write_inputs=False, fluorinate = True)

        use_decorated_spacers = True
        print "good_spacer        = ", good_spacer
        print "good_spacerf       = ", good_spacerf
        print "self.allspacername = ", self.allspacername
#       Commented line from old-version with if statement on use_decorated_spacers
#       allthespacers = good_spacer + good_spacerf + self.allspacername

        # When decorating this line prevents redundant spacer listings
        allthespacers = good_spacer + good_spacerf
        print "length allthespacers = ", len(allthespacers)

        # If no spacers listed this sets a default of allthespacers = ['']
        # so entire enumerate doenst fail
        if (len(allthespacers) == 0):
            allthespacers = self.allspacername

        print "All spacers: ", allthespacers

        ## exit here if we're JUST doing building blocks
        if ("BB" in classes and len(classes)==1):
            print "Generating JUST input for lone donor and acceptor calcs, returning"
            return

        good_don2 = self.try_set(self.alldon2, self.alldonname2, "donor", options, write_inputs=False, fluorinate = False)
        good_acc2 = self.try_set(self.allacc2, self.allaccname2, "acceptor", options, write_inputs=False, fluorinate = False)
        if (good_don2 == []):
            good_don2 = good_don
        if (good_acc2 == []):
            good_acc2 = good_acc
        print "Building block decoration done; SETUP: Tried %d, successfully built %d" % (self.tries, self.successes)
        for d in good_don:
            print d
        for a in good_acc:
            print a
        if (good_donf == []):
            print "no fluorinated donors"
        else:
            for d in good_donf:
                print d
        if (good_accf == []):
            print "no fluorinated acceptors"
        else:
            for a in good_accf:
                print a
        for t in good_term:
            print t
        self.reset_counters()

        if "D" in classes:
            self.try_class([good_don], [0], allthespacers, options)
        if "A" in classes:
            self.try_class([good_acc], [0], allthespacers, options)
        if "DA" in classes:
            self.try_class([good_don, good_acc], [0,1], allthespacers, options)
        if "ADA" in classes:
            self.try_class([good_acc, good_don], [0,1,0], allthespacers, options)
        if "AAD" in classes:
            self.try_class([good_acc, good_don], [1,1,0], allthespacers, options)
        if "DAD" in classes:
            self.try_class([good_don, good_acc], [0,1,0], allthespacers, options)
        if "DDA" in classes:
            self.try_class([good_don, good_acc], [0,0,1], allthespacers, options)
        if "ADA2" in classes:
            self.try_class([good_acc, good_don, good_acc2], [0,1,2], allthespacers, options)
        if "DAD2" in classes:
            self.try_class([good_don, good_acc, good_don2], [0,1,2], allthespacers, options)
        if "DD" in classes:
            self.try_class([good_don, good_don2], [0,1], allthespacers, options)
        if "AA" in classes:
            self.try_class([good_acc, good_acc2], [0,1], allthespacers, options)

## classes for fluorination

        print "good_donf = ", good_donf

        if "DfD" in classes:
            self.try_class([good_donf, good_don], [0,1], allthespacers, options)
        if "DfDf" in classes:
            self.try_class([good_donf, good_donf], [0,1], allthespacers, options)
        if "DfA" in classes:
            self.try_class([good_donf, good_acc], [0,1], allthespacers, options)

        if "AfD" in classes:
            self.try_class([good_accf, good_don], [0,1], allthespacers, options)
        if "AfA" in classes:
            self.try_class([good_accf, good_acc], [0,1], allthespacers, options)

        if "DfADf" in classes:
            self.try_class([good_donf, good_acc, good_donf], [0,1,0], allthespacers, options)
        if "DAfD" in classes:
            self.try_class([good_don, good_accf, good_don], [0,1,0], allthespacers, options)
        if "ADfA" in classes:
            self.try_class([good_acc, good_donf, good_acc], [0,1,0], allthespacers, options)
        if "AfDfAf" in classes:
            self.try_class([good_accf, good_donf, good_accf], [0,1,0], allthespacers, options)
##

        # A new feature: "terminal" groups for small molecule exploration.  So far they are on either end:
        # TADAT, TDADT, TAA2AT, TDD2DT
        if "TDADT" in classes:            
            for iterm1 in range(0,len(good_term)):
                for idon1 in range(0,len(good_don)):
                    for iacc in range(0,len(good_acc)):                                    
                        for ispc in range(0,len(allthespacers)):
                            sp = "%s" % allthespacers[ispc]  ## one of the spacer names can be "", for no spacer
                            str = "%s %s %s %s %s %s %s" % (good_term[iterm1], good_don[idon1], sp, good_acc[iacc], sp, good_don[idon1], good_term[iterm1])
                            print "TRYING %s" % str
                            self.try_str(str, options, True)
                                    

        if "TADAT" in classes:            
            for iterm1 in range(0,len(good_term)):
                for iacc1 in range(0,len(good_acc)):
                    for idon in range(0,len(good_don)):                                    
                        for ispc in range(0,len(allthespacers)):
                            sp = "%s" % allthespacers[ispc]  ## one of the spacer names can be "", for no spacer
                            str = "%s %s %s %s %s %s %s" % (good_term[iterm1], good_acc[iacc1], sp, good_don[idon], sp, good_acc[iacc1], good_term[iterm1])
                            print "TRYING %s" % str
                            self.try_str(str, options, True)

        if "TDD2DT" in classes:            
            for iterm1 in range(0,len(good_term)):
                for idon1 in range(0,len(good_don)):
                    for idon2 in range(0,len(good_don)):
                        if (idon1 != idon2):
                            for ispc in range(0,len(allthespacers)):
                                sp = "%s" % allthespacers[ispc]  ## one of the spacer names can be "", for no spacer
                                str = "%s %s %s %s %s %s %s" % (good_term[iterm1], good_don[idon1], sp, good_don[idon2], sp, good_don[idon1], good_term[iterm1])
                                print "TRYING %s" % str
                                self.try_str(str, options, True)

        if "TDfD2fDfT" in classes:
            for iterm1 in range(0,len(good_term)):
                for idon1 in range(0,len(good_donf)):
                    for idon2 in range(0,len(good_donf)):
                        if (idon1 != idon2):
                            for ispc in range(0,len(allthespacers)):
                                sp = "%s" % allthespacers[ispc]  ## one of the spacer names can be "", for no spacer
                                str = "%s %s %s %s %s %s %s" % (good_term[iterm1], good_don[idon1], sp, good_don[idon2], sp, good_don[idon1], good_term[iterm1])
                                print "TRYING %s" % str
                                self.try_str(str, options, True)

        if "TAA2AT" in classes:            
            for iterm1 in range(0,len(good_term)):
                for iacc1 in range(0,len(good_acc)):
                    for iacc2 in range(0,len(good_acc)):
                        if (iacc1 != iacc2):
                            for ispc in range(0,len(allthespacers)):
                                sp = "%s" % allthespacers[ispc]  ## one of the spacer names can be "", for no spacer
                                str = "%s %s %s %s %s %s %s" % (good_term[iterm1], good_acc[iacc1], sp, good_acc[iacc2], sp, good_acc[iacc1], good_term[iterm1])
                                print "TRYING %s" % str
                                self.try_str(str, options, True)

        print "FINAL_ASSEMBLY: Tried %d, successfully built %d" % (self.tries, self.successes)


def main():
    import copy
    options, args = get_options()
    
    # load the building blocks
    enumerator = BuildingBlockEnumerator(options.bblocks_dir, options.subsets_file)

    ## would be cool:
    #    enumerator.enum_template("D* (R0 R0) A*", options)
    ## but we can do better:
    enumerator.enum_smarter(options)
    good2 = copy.deepcopy(enumerator.good_strs)
    print "Found %d valid structures" % (len(good2))
    for g in good2:
        print g

if __name__=="__main__":
    main()

#    print "Found %d good ones" % (len(good1))
#    for g in good1:
#        print g
#    for g in good1:
#        if g not in good2:
#            print g , " not found by smart enum"
