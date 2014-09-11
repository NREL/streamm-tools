import os,sys,fileinput,string
import json
from math import *
from os import system
from string import center

# Arguments: 
# argv[1] = File name
# This script is designed to read through a gaussian output file (*.out)
# and extract information about the molecule, including:
# Total energy
# occupied and unoccupied levels, in eV
# Excitation energies/wavelengths (ev/nm) and oscillator strengths

# Convert a.u. to eV
auenev = float(27.21138386)

# This information is written to files, and used to generate plots of
# absorption spectra with gnuplot

def ns_get(fline):
    """ Parses a line containing 'td=nstates=ns/' embedded in other
        text to find and return the number of excited states from a
        Gaussian TDDFT calculation """
    # line has format * td=nstates=ns/* where ns = number of states
    l = fline[fline.index('nstates='):] # Remove everything before 'nstates='
    ll = l[l.index('='):] # Remove everything before =ns*
    lll = ll[:ll.index('/')] # Remove / and everything after
    ns = int(lll.strip('=')) # Remove the = and assign to integer ns
    return ns

def spec_get(fline):
    """ Parses a line containing "Excited State" to return the energy gap
        (eV), the absorption wavelength (nm) and the oscillator strength """
    st_list = line.split()
    s = st_list[8] # The 9th element is f=XXX
    st_list[8] = s[2:] # Remove the f= part from this element

    (st_list[4],st_list[6],st_list[8])
    egap = float(st_list[4])
    lam = float(st_list[6])
    fosc= float(st_list[8])
    return (egap,lam,fosc)

def check_get(fline):
    """Parses a line containing "%chk=" to get the root part of the checkpoint
       file name.
       For example, "%chk=file.chk" returns the string "file"
    """
    l = fline[fline.index('=')+1:]
    ll = l[:l.index('.chk')] 
    return ll

##############################################################################
##############################################################################
# Begin main part of the program
##############################################################################
##############################################################################
usage = "python tdfilt.py filename"
if (len(sys.argv) != 3):
    print "usage: python tdfilt.py filename.log filename.json"
    sys.exit("Incorrect argument list")

fname = sys.argv[1]
jname = sys.argv[2]
print "jname = ",jname

str_list = ["\nOpening file ","\"",fname,"\"\n"]
print ''.join(str_list)
ftddft = open(fname,"r")

# Read through a gaussian output file (*.out) to get information about
# the energy levels and generate an absorption spectrum

stinfo = []
a_occvals = []
b_occvals = []
a_virtvals = []
b_virtvals = []
converged = False
betas = False
ns=0
timedep = True

for line in ftddft.readlines():

    if 'td=nstates=' in line.lower(): # Find the number of excited states
        timedep = True
        llow = line.lower()
        ns = ns_get(llow)
        print "ns = ", ns
        # Note that this will find two lines with nstates in normal output

    if 'Excited State' in line: # Find information about excited states
        ex_info = (egap,lam,fosc) = spec_get(line)
#        print "%.4f %.2f %.4f" % ex_info # ev lam osc.strength
        stinfo.append(ex_info)

    if '%chk=' in line: # Extract the name of the checkpoint file (w/o .chk)
        cfile = check_get(line.lower())
        print cfile,"\n"

    if 'Optimized Parameters' in line or timedep==True: converged=True
    if timedep==True: converged=True

    if 'Tot=' in line: # Get the dipole vector and length (in Debye)
        dip_list = line.split()
        x = center("X",10)
        y = center("Y",10)
        z = center("Z",10)
        tot = center("Total",10)
        dx = center(dip_list[1],10)
        dy = center(dip_list[3],10)
        dz = center(dip_list[5],10)
        dtot = center(dip_list[7],10)
        
# Get the occupied and unoccupied levels
# (after convergence has been assured)
    if 'occ.' in line and converged == True:
        lsplit = line.split()
        if 'Alpha' in line:
            for v in lsplit[4:]: a_occvals.append(auenev*float(v))
        if 'Beta' in line:
            betas = True
            for v in lsplit[4:]: b_occvals.append(auenev*float(v))
    if 'virt.' in line and converged == True:
        lsplit = line.split()
        if 'Alpha' in line:
            for v in lsplit[4:]: a_virtvals.append(auenev*float(v))
        if 'Beta' in line:
            for v in lsplit[4:]: b_virtvals.append(auenev*float(v))

# Get the total energy, after convergence has been assured
    if 'B3LYP)' in line and converged==True:
#        print >> eout , line
        esplit = line.split()
        etot = float(esplit[4])
        
ftddft.close()

# Get working directory to provide absolute paths in the gnuplot stuff below
# This was added during debugging and does not seem to be strictly necessary 
#workdir = os.getcwd()
# A working directory "." seems to work fine
workdir = "."

# if nstates= found, then do the spectral output 
#print "timedep = ",timedep
if timedep:

# Create file to have the HOMO, LUMO, Gap, Optical LUMO, total energy
    fall = open("".join((cfile,"_all_ens.out")),"w")
    thomo = a_occvals[-1]
    tlumo = a_virtvals[0]
    (tgap, tlam, tosc) = stinfo[0]
    print >> fall, "HOMO   LUMO   Gap    Opt.LUMO    Tot.En"
    print >> fall, "%.4f %.4f %.4f %.4f %.6f" % (thomo, tlumo, tgap, thomo+tgap,etot)
#    print "%.4f %.4f %.4f %.4f %.6f" % (thomo, tlumo, tgap, thomo+tgap,etot)
    fall.close()

# Get information about HOMO and LUMO contours that contain 80% of the 
# electron density (generated by ~/OPENCUBMAN/fetocv *MO.cube 0.8 0 > *data)
    fhdata = open("hdata","r")
    for line in fhdata.readlines():
        hsplit = line.split()
        hval = float(hsplit[3])
    fhdata.close()

    fldata = open("ldata","r")
    for line in fldata.readlines():
        lsplit = line.split()
        lval = float(lsplit[3])
    fldata.close()
    
    hval = sqrt(hval*hval)
    lval = sqrt(lval*lval)
#    print "hval and lval = ",hval,lval
    
# Write out all desired information to the metadata file
    properties = {
        'homo':thomo,
        'lumo':tlumo,
        'gap':tgap,
        'olumo':thomo+tgap,
        'etot':etot,
        'homo_contour':hval,
        'lumo_contour':lval,
        'mu_tot':float(dtot),
        'mu_x':float(dx),
        'mu_y':float(dy),
        'mu_z':float(dz),
    }


#   add results to the metadata properties json file
    f = open(jname, 'r')
    doc = json.load(f)
    f.close()
    doc['results'] = properties
    f = open(jname, 'w')
    json.dump(doc, f, indent=2)
    f.close()
