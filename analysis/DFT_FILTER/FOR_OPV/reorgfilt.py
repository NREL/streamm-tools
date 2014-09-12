import os,sys,fileinput,string
import time, datetime
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
usage = "python tdfilt.py filename doc.json"
#print len(sys.argv)
#for i in sys.argv:
#    print i
#print 
if (len(sys.argv) != 3):
    print "usage: python reorgfilt.py fname.log jname.json"
    sys.exit("Incorrect argument list")

fname = sys.argv[1]
jname = sys.argv[2]
print "jname = ",jname

str_list = ["\nOpening file ","\"",fname,"\"\n"]
print ''.join(str_list)
ftddft = open(fname,"r")

# Read through a gaussian output file (*.log) to get information about
# the energy of the cation for the initial geometry, for the cation geometry
# and of the neutral for the cation geometry

stinfo = []
first = 0

for line in ftddft.readlines():

    if '%chk=' in line: # Extract the name of the checkpoint file (w/o .chk)
        cfile = check_get(line.lower())
        #print cfile,"\n"

    if 'E(UB3LYP)' in line:
        esplit = line.split()
        cat_en = float(esplit[4])
        first = first + 1
        if (first==1):
            cat_en_first = cat_en

    if 'Optimized Parameters' in line:
        converged=True
        cat_en_final = cat_en

# Get the total energy, after convergence has been assured
    if 'E(RB3LYP)' in line:
        esplit = line.split()
        etot = float(esplit[4])

ftddft.close()

# Write out all desired information to the metadata file
cat_props = {
        'vert_cat_en':cat_en_first,
        'relaxed_cat_en':cat_en_final,
        'relaxed_neut_en':etot,
}

#   add cation results and data to the metadata properties json file
f = open(jname, 'r')
doc = json.load(f)
print 'Initial json file'
print doc
f.close()

now = datetime.datetime.now()
doc['metadata']['date_cation_calc_finished'] = now.isoformat()

doc['cation_results'] = cat_props

#jjname = ''.join(['NEW_',jname])
jjname = jname
f = open(jjname, 'w')
json.dump(doc, f, indent=2)
f.close()

