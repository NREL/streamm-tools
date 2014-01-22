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
    kns = int(lll.strip('=')) # Remove the = and assign to integer ns
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

def seconv_data(lmin,lmax,dl,ewidth,stinfo,spectrum):
    """Convolve the stick spectrum (lam,fosc) with a gaussian defined by lwidth
       exp(-(E-Estick)**2/2*ewidth**2)
       spec(l) = sum_lam fosc*exp(-(1/l-1/lam)**2/(1242/ewidth)**2/2)*scale
       scale = 2.64/ewidth if ewidth is in units of eV
       this gives output as x 10^4 L / mol cm
       Print (l,spec(l)) to a file called spect_lam.out for printing by gnuplot
    """
    nlen = int((lmax-lmin)/dl)
    l_it = range(nlen+1)
    for i in l_it:
        l_it[i] = dl*float(i) + lmin

    for l in l_it:
        spect = 0.0
        for (egap,lam,fosc) in stinfo:
            dlam = (1/l - 1/lam)*(1242/ewidth)
            gaussl = exp(-dlam**2 / 2.)
            spect = spect + fosc*gaussl*2.64/ewidth
        callout = (l,spect)
        spectrum.append(callout)

##############################################################################
##############################################################################
# Begin main part of the program
##############################################################################
##############################################################################
usage = "python tdfilt.py filename"
#print len(sys.argv)
#for i in sys.argv:
#    print i
#print 
if (len(sys.argv) != 7):
    print "usage: python tdfilt.py fname.log jname.json energy width lmin lmax dl"
    sys.exit("Incorrect argument list")

fname = sys.argv[1]
jname = sys.argv[2]
print "jname = ",jname

ewidth = float(sys.argv[3])
lmin = float(sys.argv[4])
lmax = float(sys.argv[5])
dl = float(sys.argv[6])
print ewidth,lmin,lmax,dl

#lwidth = fwhm/(2.*log(2.))
#print fwhm,lwidth

str_list = ["\nOpening file ","\"",fname,"\"\n"]
print ''.join(str_list)
ftddft = open(fname,"r")

# Read through a gaussian output file (*.out) to get information about
# the energy levels and generate an absorption spectrum

stinfo = []
gaplist = []
lamlist = []
fosclist = []
a_occvals = []
b_occvals = []
a_virtvals = []
b_virtvals = []
converged = False
betas = False
ns=0
timedep = False

for line in ftddft.readlines():

    if 'td=nstates=' in line.lower(): # Find the number of excited states
        timedep = True
        llow = line.lower()
        ns = ns_get(llow)
        #print 'System has %d excited states ' % ns
        # Note that this will find two lines with nstates in normal output

    if 'Excited State' in line: # Find information about excited states
        ex_info = (egap,lam,fosc) = spec_get(line)
        #print "%.4f %.2f %.4f" % ex_info # ev lam osc.strength
        stinfo.append(ex_info)
        gaplist.append(egap)
        lamlist.append(lam)
        fosclist.append(fosc)

    if '%chk=' in line: # Extract the name of the checkpoint file (w/o .chk)
        cfile = check_get(line.lower())
        #print cfile,"\n"

    if 'Optimized Parameters' in line or timedep==True: converged=True

    if timedep==True: converged=True

    if ('Tot=' in line and 'X=' in line and converged==True): # Get the dipole vector for the ground state and length (in Debye)
        dip_list = line.split()
        x = center("X",10)
        y = center("Y",10)
        z = center("Z",10)
        tot = center("Total",10)
        dx = center(dip_list[1],10)
        dy = center(dip_list[3],10)
        dz = center(dip_list[5],10)
        dtot = center(dip_list[7],10)
        print dx, dy, dz, dtot

    if ('Tot=' in line and 'X=' in line and converged==True and timedep==True): # Get the dipole vector for the first excited state and length (in Debye)
        dip_list = line.split()
        x = center("X",10)
        y = center("Y",10)
        z = center("Z",10)
        tot = center("Total",10)
        d2x = center(dip_list[1],10)
        d2y = center(dip_list[3],10)
        d2z = center(dip_list[5],10)
        d2tot = center(dip_list[7],10)
        #print d2x, d2y, d2z, d2tot
        
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
    if 'B3LYP)' in line:
#        print >> eout , line
        esplit = line.split()
        etot = float(esplit[4])

    if 'wPBE)' in line:
#        print >> eout , line
        esplit = line.split()
        etot = float(esplit[4])

ftddft.close()

# Define the HOMO, LUMO, Gap, Optical LUMO
thomo = a_occvals[-1]
tlumo = a_virtvals[0]
(tgap, tlam, tosc) = stinfo[0]
opt_lumo = thomo + tgap

#sticks = {'sticklabels':['Energy gap','Wavelength','Oscillator Strength'],'sticklist':stinfo}


sticks = {'Energy gap (eV)':gaplist,'Wavelength (nm)':lamlist,'Oscillator strength':fosclist}

#print sticks

#print "HOMO LUMO GAP OPT_LUMO"
#print thomo, tlumo, tgap, opt_lumo
#print
print 'Total energy = ', etot

# Call the function to create the convolved spectrum "cfile_spect_lam.out"
spectrum = []
seconv_data(lmin,lmax,dl,ewidth,stinfo,spectrum)

#print spectrum


spectra = {'broadening':ewidth,'xlabel':'Wavelength (nm)','ylabel':'Absorption (x 10^4 L/mol cm)','spectdata':spectrum}

#print spectra

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
print "hval and lval = ",hval,lval

# Write out all desired information to the metadata file
properties = {
        'homo':thomo,
        'lumo':tlumo,
        'gap':tgap,
        'optical_lumo':thomo+tgap,
        'homo_contour':hval,
        'lumo_contour':lval,
        #'mu_tot':float(dtot),
        #'mu_x':float(dx),
        #'mu_y':float(dy),
        #'mu_z':float(dz),
        #'mu2_tot':float(d2tot),
        #'mu2_x':float(d2x),
        #'mu2_y':float(d2y),
        #'mu2_z':float(d2z),
        'total energy':etot,
}

#   add results, sticks, and spectrum to the metadata properties json file
f = open(jname, 'r')
doc = json.load(f)
print 'Initial json file'
print doc
f.close()

now = datetime.datetime.now()
doc['metadata']['date_calculation_finished'] = now.isoformat()

doc['results'] = properties
doc['sticks'] = sticks
doc['spectrum'] = spectra

#jjname = ''.join(['NEW_',jname])
jjname = jname
f = open(jjname, 'w')
json.dump(doc, f, indent=2)
f.close()

