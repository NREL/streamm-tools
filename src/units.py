#! /usr/bin/env python
"""
Unit conversions

File will be read in and all units will be converted to the default units

   distance - Angstroms 
   mass - AMU
   energy - kcal/mol
   charge - electron
   density - AMU Angstroms^-3

Upon output default units will be converted to the appropriate units,
which will be displayed in the header of the output file
or given for a certian type of input file 

Conversions are taken from the

CRC Handbook of Chemistry and Physics 95th see CRC_HCP_energy.pdf in doc director
"""

# Distance 
nm_angstroms = 10.0
bohr_angstrom = 0.5291772086

# Energy 
kJmol_kcalmol = 0.239006
kJmol_eV = 0.01036427
kJmol_H = 0.0003808798
kcalmol_eV = 0.04336411 
kcalmol_H = 0.001593601 
eV_H = 0.03674931 
#kJmol_incm = 83.5935
#eV_kcalmol = 96.4853


# Constants 
# http://physics.nist.gov/cgi-bin/cuu/Value?na
const_avo = 6.02214129 # x10^23 mol^-1 


# grosig_kcalmol = 5.61230943


def convert_bohr_ang(dist_bohr):
    """
    Convert bohr to angstroms
    """
    return dist_bohr*bohr_angstrom


def convert_ang_bohr(dist_ang):
    """
    Convert angstroms to bohr 
    """
    return dist_ang/bohr_angstrom
        

def convert_gcm3_AMUA3(den_gmc3):
    """
    Convert density from g/cm^3 to AMU Angstroms^-3
    """
    return den_gmc3*const_avo/10.0
    

def convert_AMUA3_gcm3(den_AMUA3):
    """
    Convert density from AMU Angstroms^-3 to  g/cm^3
    """
    return den_AMUA3/const_avo*10.0
    
def convert_nm_angstroms(d_nm):
    """
    convert nm to angstroms
    """
    return d_nm*nm_angstroms
    
def convert_angstroms_nm(d_angstroms):
    """
    convert angstroms to nm 
    """
    return d_angstroms/nm_angstroms

    
def convert_kJmol_kcalmol(en_kJmol):
    """
    convert kJmol to kcalmol
    """
    return en_kJmol*kJmol_kcalmol
    
def convert_kJmol_eV(en_kJmol):
    """
    convert kJmol to eV
    """
    return en_kJmol*kJmol_eV
    

def convert_kJmol_H(en_kJmol):
    """
    convert kJmol to Hartree
    """
    return en_kJmol*kJmol_H

def convert_kcalmol_kJmol(en_kcalmol):
    """
    convert kcalmol  to kJmol
    """
    return en_kcalmol/kJmol_kcalmol
    
def convert_kcalmol_eV(en_kcalmol):
    """
    convert kcalmol to eV
    """
    return en_kcalmol*kcalmol_eV

def convert_kcalmol_H(en_kcalmol):
    """
    convert kcalmol to Hartree
    """
    return en_kcalmol*kcalmol_H

def convert_eV_kJmol(en_eV):
    """
    convert eV to kJmol
    """
    return en_eV/kJmol_eV
    
def convert_eV_kcalmol(en_eV):
    """
    convert eV to kcalmol
    """
    return en_eV/kcalmol_eV
    
def convert_eV_H(en_eV):
    """
    convert eV to Hartree
    """
    return en_eV*eV_H
    

def convert_H_kJmol(en_H):
    """
    convert Hartree to kJmol
    """
    return en_H/kJmol_H


def convert_H_kcalmol(en_H):
    """
    convert Hartree to kcalmol
    """
    return en_H/kcalmol_H


def convert_H_eV(en_H):
    """
    convert Hartree to kcalmol
    """
    return en_H/eV_H


    
#def convert_kJmol_incm(en_kJmol):
#    """
#    convert kJmol to cm^-1
#    """
#    return en_kJmol*kJmol_incm
    
#def convert_gro_sig(sig_gro):
#    """
#    convert epsilon from gromacs to kcal/mol 
#    """
#    return sig_gro*grosig_kcalmol
    
    

def convert_g_bond_kb(g_kb):
    """
    convert gromacs harmonic bond parameter in
       2.0 kJ/mol /nm /nm
    to general bond parameter defined by lammps
       kcal/mol /angstrom/angstrom
    """
    return g_kb*kJmol_kcalmol/nm_angstroms/nm_angstroms/2.0


def convert_kb_g_bond(g_kb):
    """
    convert to gromacs harmonic bond parameter in
       2.0 kJ/mol /nm /nm
    from general bond parameter defined by lammps
       kcal/mol /angstrom/angstrom
    """
    return g_kb/kJmol_kcalmol*nm_angstroms*nm_angstroms*2.0



def convert_g_angle_kb(g_ka):
    """
    convert gromacs harmonic angle parameter in
       1/2.0 K kJ/mol/radina^2
    to general angle parameter defined by lammps
       kcal/mol /radina^2
    """
    return g_ka*kJmol_kcalmol/2.0

    
def convert_kb_g_angle(g_ka):
    """
    convert to gromacs harmonic angle parameter in
        1/2.0 K kJ/mol/radina^2
    from general angle parameter defined by lammps
        kcal/mol /radina^2
    """
    return g_ka/kJmol_kcalmol*2.0
