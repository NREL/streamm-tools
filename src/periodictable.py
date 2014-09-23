"""
Class of elements with associated atomic properties 
"""

# Dr. Travis Kemper
# NREL
# Initial Date 7/01/2014
# travis.kemper@nrel.gov

import numpy 

class element():
    """
    Gives the basic properties of the elements in the periodic table
    """


    def __init__(self,symbol,number,mass,cov_radii,vdw_radii):
        """
        Constructor for a element
        
        Args:
        symbol,number,ionic_radii
        """

        if isinstance(symbol, str):
            self.symbol = symbol
        else:
            print "1st arg should be str"
            raise TypeError


        if isinstance(number, int):
            self.number = number
        else:
            print "2nd arg should be int"
            raise TypeError

        if isinstance(mass, float):
            self.mass = mass
        else:
            print "3rd arg should be float"
            raise TypeError


        if isinstance(cov_radii, float):
            self.cov_radii = cov_radii
        else:
            print "4th arg should be float"
            raise TypeError


        if isinstance(vdw_radii, float):
            self.vdw_radii = vdw_radii
        else:
            print "5th arg should be float"
            raise TypeError

    def __del__(self):
        """
        Destructor, clears object memory
        """
        
        del self.symbol 
        del self.number 
        del self.mass 
        del self.cov_radii
        del self.vdw_radii

    def symbol(self):
        return self.symbol

    def mass(self):
        return self.mass

    def cov_radii(self):
        return self.cov_radii

    def vdw_radii(self):
        return self.vdw_radii
    

class periodictable:
    """
    Elements of the periodic table 
    """

    def __init__(self):
        """
        Set elements of the periodic table 
        """
        
        self.elements=dict()
        self.maxgid = 110

        elsymbol = ['??']*(self.maxgid+1)

        elsymbol[1] =   'H'  
        elsymbol[2] =   'He'
        elsymbol[3] =   'Li'  
        elsymbol[4] =   'Be'  
        elsymbol[5] =   'B'  
        elsymbol[6] =   'C'  
        elsymbol[7] =   'N'  
        elsymbol[8] =   'O'  
        elsymbol[9] =   'F'  
        elsymbol[10] =   'Ne'  
        elsymbol[11] =   'Na'  
        elsymbol[12] =   'Mg'  
        elsymbol[13] =   'Al'  
        elsymbol[14] =   'Si'  
        elsymbol[15] =   'P'  
        elsymbol[16] =   'S'  
        elsymbol[17] =   'Cl'  
        elsymbol[18] =   'Ar'  
        elsymbol[19] =   'K'  
        elsymbol[20] =   'Ca'  
        elsymbol[21] =   'Sc'  
        elsymbol[22] =   'Ti'  
        elsymbol[23] =   'V'  
        elsymbol[24] =   'Cr'  
        elsymbol[25] =   'Mn'  
        elsymbol[26] =   'Fe'  
        elsymbol[27] =   'Co'  
        elsymbol[28] =   'Ni'  
        elsymbol[29] =   'Cu'  
        elsymbol[30] =   'Zn'  
        elsymbol[31] =   'Ga'  
        elsymbol[32] =   'Ge'  
        elsymbol[33] =   'As'  
        elsymbol[34] =   'Se'  
        elsymbol[35] =   'Br'  
        elsymbol[36] =   'Kr'  
        elsymbol[37] =   'Rb'  
        elsymbol[38] =   'Sr'  
        elsymbol[39] =   'Y'  
        elsymbol[40] =   'Zr'  
        elsymbol[41] =   'Nb'  
        elsymbol[42] =   'Mo'  
        elsymbol[43] =   'Tc'  
        elsymbol[44] =   'Ru'  
        elsymbol[45] =   'Rh'  
        elsymbol[46] =   'Pd'  
        elsymbol[47] =   'Ag'  
        elsymbol[48] =   'Cd'  
        elsymbol[49] =   'In'  
        elsymbol[50] =   'Sn'  
        elsymbol[51] =   'Sb'  
        elsymbol[52] =   'Te'  
        elsymbol[53] =   'I'  
        elsymbol[54] =   'Xe'  
        elsymbol[55] =   'Cs'  
        elsymbol[56] =   'Ba'  
        elsymbol[71] =   'Lu'  
        elsymbol[72] =   'Hf'  
        elsymbol[73] =   'Ta'  
        elsymbol[74] =   'W'  
        elsymbol[75] =   'Re'  
        elsymbol[76] =   'Os'  
        elsymbol[77] =   'Ir'  
        elsymbol[78] =   'Pt'  
        elsymbol[79] =   'Au'  
        elsymbol[80] =   'Hg'  
        elsymbol[81] =   'Tl'  
        elsymbol[82] =   'Pb'  
        elsymbol[83] =   'Bi'  
        elsymbol[84] =   'Po'  
        elsymbol[85] =   'At'  
        elsymbol[86] =   'Rn'  
        elsymbol[58] =   'Ce'  
        elsymbol[66] =   'Dy'  
        elsymbol[68] =   'Er'  
        elsymbol[63] =   'Eu'  
        elsymbol[64] =   'Gd'  
        elsymbol[67] =   'Ho'  
        elsymbol[57] =   'La'  
        elsymbol[60] =   'Nd'  
        elsymbol[61] =   'Pm'  
        elsymbol[59] =   'Pr'  
        elsymbol[62] =   'Sm'  
        elsymbol[65] =   'Tb'  
        elsymbol[69] =   'Tm'  
        elsymbol[70] =   'Yb'  
        elsymbol[87] =   'Fr'  
        elsymbol[88] =   'Ra'  
        elsymbol[104] =   'Rf'  
        elsymbol[105] =   'Db'  
        elsymbol[106] =   'Sg'  
        elsymbol[107] =   'Bh'  
        elsymbol[108] =   'Hs'  
        elsymbol[109] =   'Mt'  
        elsymbol[110] =   'Ds'  
        elsymbol[89] =   'Ac'  
        elsymbol[95] =   'Am'  
        elsymbol[97] =   'Bk'  
        elsymbol[98] =   'Cf'  
        elsymbol[96] =   'Cm'  
        elsymbol[99] =   'Es'  
        elsymbol[100] =   'Fm'  
        elsymbol[101] =   'Md'  
        elsymbol[102] =   'No'  
        elsymbol[93] =   'Np'  
        elsymbol[91] =   'Pa'  
        elsymbol[94] =   'Pu'  
        elsymbol[90] =   'Th'  
        elsymbol[92] =   ' U'

        ELMASS = numpy.zeros(self.maxgid+1)

        ELMASS[1] =  1.0080 #  H
        ELMASS[2] =  4.0030 #  He
        ELMASS[3] =  6.9410 #  Li
        ELMASS[4] =  9.0120 #  Be
        ELMASS[5] =  10.8110 #  B
        ELMASS[6] =  12.0110 #  C
        ELMASS[7] =  14.0070 #  N
        ELMASS[8] =  15.9990 #  O
        ELMASS[9] =  18.9980 #  F
        ELMASS[10] =  20.1800 #  Ne
        ELMASS[11] =  22.9910 #  Na
        ELMASS[12] =  24.3050 #  Mg
        ELMASS[13] =  26.9820 #  Al
        ELMASS[14] =  28.0860 #  Si
        ELMASS[15] =  30.9740 #  P
        ELMASS[16] =  32.0660 #  S
        ELMASS[17] =  35.4530 #  Cl
        ELMASS[18] =  39.9480 #  Ar
        ELMASS[19] =  39.0980 #  K
        ELMASS[20] =  40.0780 #  Ca
        ELMASS[21] =  44.9560 #  Sc
        ELMASS[22] =  47.8670 #  Ti
        ELMASS[23] =  50.9420 #  V
        ELMASS[24] =  51.9960 #  Cr
        ELMASS[25] =  54.9380 #  Mn
        ELMASS[26] =  55.8450 #  Fe
        ELMASS[27] =  58.9330 #  Co
        ELMASS[28] =  58.6930 #  Ni
        ELMASS[29] =  63.5460 #  Cu
        ELMASS[30] =  65.3900 #  Zn
        ELMASS[31] =  69.7230 #  Ga
        ELMASS[32] =  72.6100 #  Ge
        ELMASS[33] =  74.9220 #  As
        ELMASS[34] =  78.9600 #  Se
        ELMASS[35] =  79.9040 #  Br
        ELMASS[36] =  83.8000 #  Kr
        ELMASS[37] =  85.4680 #  Rb
        ELMASS[38] =  87.6200 #  Sr
        ELMASS[39] =  88.9060 #  Y
        ELMASS[40] =  91.2240 #  Zr
        ELMASS[41] =  92.9060 #  Nb
        ELMASS[42] =  95.9400 #  Mo
        ELMASS[43] =  98.00 #  Tc
        ELMASS[44] =  101.0700 #  Ru
        ELMASS[45] =  102.9060 #  Rh
        ELMASS[46] =  106.4200 #  Pd
        ELMASS[47] =  107.8680 #  Ag
        ELMASS[48] =  112.4110 #  Cd
        ELMASS[49] =  114.8180 #  In
        ELMASS[50] =  118.710 #  Sn
        ELMASS[51] =  121.7600 #  Sb
        ELMASS[52] =  127.6000 #  Te
        ELMASS[53] =  126.9040 #  I
        ELMASS[54] =  131.2900 #  Xe
        ELMASS[55] =  132.9050 #  Cs
        ELMASS[56] =  137.3270 #  Ba
        ELMASS[71] =  174.9670 #  Lu
        ELMASS[72] =  178.4900 #  Hf
        ELMASS[73] =  180.9480 #  Ta
        ELMASS[74] =  183.8400 #  W
        ELMASS[75] =  186.2070 #  Re
        ELMASS[76] =  190.2300 #  Os
        ELMASS[77] =  192.2170 #  Ir
        ELMASS[78] =  195.0780 #  Pt
        ELMASS[79] =  196.9670 #  Au
        ELMASS[80] =  200.5900 #  Hg
        ELMASS[81] =  204.3830 #  Tl
        ELMASS[82] =  207.2000 #  Pb
        ELMASS[83] =  208.9800 #  Bi
        ELMASS[84] =  210.00 #  Po
        ELMASS[85] =  210.00 #  At
        ELMASS[86] =  222.00 #  Rn
        ELMASS[58] =  140.1160 #  Ce
        ELMASS[66] =  162.5000 #  Dy
        ELMASS[68] =  167.2600 #  Er
        ELMASS[63] =  151.9640 #  Eu
        ELMASS[64] =  157.2500 #  Gd
        ELMASS[67] =  164.9300 #  Ho
        ELMASS[57] =  138.9060 #  La
        ELMASS[60] =  144.2400 #  Nd
        ELMASS[61] =  145.00 #  Pm
        ELMASS[59] =  140.9080 #  Pr
        ELMASS[62] =  150.3600 #  Sm
        ELMASS[65] =  158.9250 #  Tb
        ELMASS[69] =  168.9340 #  Tm
        ELMASS[70] =  173.0400 #  Yb
        ELMASS[87] =  223.00 #  Fr
        ELMASS[88] =  226.00 #  Ra
        ELMASS[104] =  261.00 #  Rf
        ELMASS[105] =  262.00 #  Db
        ELMASS[106] =  266.00 #  Sg
        ELMASS[107] =  264.00 #  Bh
        ELMASS[108] =  269.00 #  Hs
        ELMASS[109] =  268.00 #  Mt
        ELMASS[110] =  271.00 #  Ds
        ELMASS[89] =  227.00 #  Ac
        ELMASS[95] =  243.00 #  Am
        ELMASS[97] =  247.00 #  Bk
        ELMASS[98] =  251.00 #  Cf
        ELMASS[96] =  247.00 #  Cm
        ELMASS[99] =  252.00 #  Es
        ELMASS[100] =  257.00 #  Fm
        ELMASS[101] =  258.00 #  Md
        ELMASS[102] =  259.00 #  No
        ELMASS[93] =  237.00 #  Np
        ELMASS[91] =  231.0360 #  Pa
        ELMASS[94] =  244.00 #  Pu
        ELMASS[90] =  232.0380 #  Th
        ELMASS[92] =  238.0290 #  U

        radi_cov = numpy.zeros(self.maxgid+1)

        radi_cov[1] =  0.230 #  H
        radi_cov[2] =  1.500 #  He
        radi_cov[3] =  1.280 #  Li
        radi_cov[4] =  0.960 #  Be
        radi_cov[5] =  0.830 #  B
        radi_cov[6] =  0.680 #  C
        radi_cov[7] =  0.680 #  N
        radi_cov[8] =  0.680 #  O
        radi_cov[9] =  0.640 #  F
        radi_cov[10] =  1.500 #  Ne
        radi_cov[11] =  1.660 #  Na
        radi_cov[12] =  1.410 #  Mg
        radi_cov[13] =  1.210 #  Al
        radi_cov[14] =  1.200 #  Si
        radi_cov[15] =  1.050 #  P
        radi_cov[16] =  1.020 #  S
        radi_cov[17] =  0.990 #  Cl
        radi_cov[18] =  1.510 #  Ar
        radi_cov[19] =  2.030 #  K
        radi_cov[20] =  1.760 #  Ca
        radi_cov[21] =  1.700 #  Sc
        radi_cov[22] =  1.600 #  Ti
        radi_cov[23] =  1.530 #  V
        radi_cov[24] =  1.390 #  Cr
        radi_cov[25] =  1.610 #  Mn
        radi_cov[26] =  1.520 #  Fe
        radi_cov[27] =  1.260 #  Co
        radi_cov[28] =  1.240 #  Ni
        radi_cov[29] =  1.320 #  Cu
        radi_cov[30] =  1.220 #  Zn
        radi_cov[31] =  1.220 #  Ga
        radi_cov[32] =  1.170 #  Ge
        radi_cov[33] =  1.210 #  As
        radi_cov[34] =  1.220 #  Se
        radi_cov[35] =  1.210 #  Br
        radi_cov[36] =  1.500 #  Kr
        radi_cov[37] =  2.200 #  Rb
        radi_cov[38] =  1.950 #  Sr
        radi_cov[39] =  1.900 #  Y
        radi_cov[40] =  1.750 #  Zr
        radi_cov[41] =  1.640 #  Nb
        radi_cov[42] =  1.540 #  Mo
        radi_cov[43] =  1.470 #  Tc
        radi_cov[44] =  1.460 #  Ru
        radi_cov[45] =  1.450 #  Rh
        radi_cov[46] =  1.390 #  Pd
        radi_cov[47] =  1.450 #  Ag
        radi_cov[48] =  1.440 #  Cd
        radi_cov[49] =  1.420 #  In
        radi_cov[50] =  1.390 #  Sn
        radi_cov[51] =  1.390 #  Sb
        radi_cov[52] =  1.470 #  Te
        radi_cov[53] =  1.400 #  I
        radi_cov[54] =  1.500 #  Xe
        radi_cov[55] =  2.440 #  Cs
        radi_cov[56] =  2.150 #  Ba
        radi_cov[71] =  1.870 #  Lu
        radi_cov[72] =  1.750 #  Hf
        radi_cov[73] =  1.700 #  Ta
        radi_cov[74] =  1.620 #  W
        radi_cov[75] =  1.510 #  Re
        radi_cov[76] =  1.440 #  Os
        radi_cov[77] =  1.410 #  Ir
        radi_cov[78] =  1.360 #  Pt
        radi_cov[79] =  1.500 #  Au
        radi_cov[80] =  1.320 #  Hg
        radi_cov[81] =  1.450 #  Tl
        radi_cov[82] =  1.460 #  Pb
        radi_cov[83] =  1.480 #  Bi
        radi_cov[84] =  1.400 #  Po
        radi_cov[85] =  1.210 #  At
        radi_cov[86] =  1.500 #  Rn
        radi_cov[58] =  2.040 #  Ce
        radi_cov[66] =  1.920 #  Dy
        radi_cov[68] =  1.890 #  Er
        radi_cov[63] =  1.980 #  Eu
        radi_cov[64] =  1.960 #  Gd
        radi_cov[67] =  1.920 #  Ho
        radi_cov[57] =  2.070 #  La
        radi_cov[60] =  2.010 #  Nd
        radi_cov[61] =  1.990 #  Pm
        radi_cov[59] =  2.030 #  Pr
        radi_cov[62] =  1.980 #  Sm
        radi_cov[65] =  1.940 #  Tb
        radi_cov[69] =  1.900 #  Tm
        radi_cov[70] =  1.870 #  Yb
        radi_cov[87] =  2.600 #  Fr
        radi_cov[88] =  2.210 #  Ra
        radi_cov[89] =  2.150 #  Ac
        radi_cov[95] =  1.800 #  Am
        radi_cov[97] =  1.540 #  Bk
        radi_cov[98] =  1.830 #  Cf
        radi_cov[96] =  1.690 #  Cm
        radi_cov[99] =  1.500 #  Es
        radi_cov[93] =  1.900 #  Np
        radi_cov[91] =  2.000 #  Pa
        radi_cov[94] =  1.870 #  Pu
        radi_cov[90] =  2.060 #  Th
        radi_cov[92] =  1.960 #  U
        radi_cov[ 100:self.maxgid ] =  1.500

        radi_vdw = numpy.zeros(self.maxgid+1)

        radi_vdw[1] =  1.090 #  H
        radi_vdw[2] =  1.400 #  He
        radi_vdw[3] =  1.820 #  Li
        radi_vdw[4] =  2.000 #  Be
        radi_vdw[5] =  2.000 #  B
        radi_vdw[6] =  1.700 #  C
        radi_vdw[7] =  1.550 #  N
        radi_vdw[8] =  1.520 #  O
        radi_vdw[9] =  1.470 #  F
        radi_vdw[10] =  1.540 #  Ne
        radi_vdw[11] =  2.270 #  Na
        radi_vdw[12] =  1.730 #  Mg
        radi_vdw[13] =  2.000 #  Al
        radi_vdw[14] =  2.100 #  Si
        radi_vdw[15] =  1.800 #  P
        radi_vdw[16] =  1.800 #  S
        radi_vdw[17] =  1.750 #  Cl
        radi_vdw[18] =  1.880 #  Ar
        radi_vdw[19] =  2.750 #  K
        radi_vdw[20] =  2.000 #  Ca
        radi_vdw[21] =  2.000 #  Sc
        radi_vdw[22] =  2.000 #  Ti
        radi_vdw[23] =  2.000 #  V
        radi_vdw[24] =  2.000 #  Cr
        radi_vdw[25] =  2.000 #  Mn
        radi_vdw[26] =  2.000 #  Fe
        radi_vdw[27] =  2.000 #  Co
        radi_vdw[28] =  1.630 #  Ni
        radi_vdw[29] =  1.400 #  Cu
        radi_vdw[30] =  1.390 #  Zn
        radi_vdw[31] =  1.870 #  Ga
        radi_vdw[32] =  2.000 #  Ge
        radi_vdw[33] =  1.850 #  As
        radi_vdw[34] =  1.900 #  Se
        radi_vdw[35] =  1.850 #  Br
        radi_vdw[36] =  2.020 #  Kr
        radi_vdw[37] =  2.000 #  Rb
        radi_vdw[38] =  2.000 #  Sr
        radi_vdw[39] =  2.000 #  Y
        radi_vdw[40] =  2.000 #  Zr
        radi_vdw[41] =  2.000 #  Nb
        radi_vdw[42] =  2.000 #  Mo
        radi_vdw[43] =  2.000 #  Tc
        radi_vdw[44] =  2.000 #  Ru
        radi_vdw[45] =  2.000 #  Rh
        radi_vdw[46] =  1.630 #  Pd
        radi_vdw[47] =  1.720 #  Ag
        radi_vdw[48] =  1.580 #  Cd
        radi_vdw[49] =  1.930 #  In
        radi_vdw[50] =  2.170 #  Sn
        radi_vdw[51] =  2.000 #  Sb
        radi_vdw[52] =  2.060 #  Te
        radi_vdw[53] =  1.980 #  I
        radi_vdw[54] =  2.160 #  Xe
        radi_vdw[55] =  2.000 #  Cs
        radi_vdw[56] =  2.000 #  Ba
        radi_vdw[71] =  2.000 #  Lu
        radi_vdw[72] =  2.000 #  Hf
        radi_vdw[73] =  2.000 #  Ta
        radi_vdw[74] =  2.000 #  W
        radi_vdw[75] =  2.000 #  Re
        radi_vdw[76] =  2.000 #  Os
        radi_vdw[77] =  2.000 #  Ir
        radi_vdw[78] =  1.720 #  Pt
        radi_vdw[79] =  1.660 #  Au
        radi_vdw[80] =  1.550 #  Hg
        radi_vdw[81] =  1.960 #  Tl
        radi_vdw[82] =  2.020 #  Pb
        radi_vdw[83] =  2.000 #  Bi
        radi_vdw[84] =  2.000 #  Po
        radi_vdw[85] =  2.000 #  At
        radi_vdw[86] =  2.000 #  Rn
        radi_vdw[58] =  2.000 #  Ce
        radi_vdw[66] =  2.000 #  Dy
        radi_vdw[68] =  2.000 #  Er
        radi_vdw[63] =  2.000 #  Eu
        radi_vdw[64] =  2.000 #  Gd
        radi_vdw[67] =  2.000 #  Ho
        radi_vdw[57] =  2.000 #  La
        radi_vdw[60] =  2.000 #  Nd
        radi_vdw[61] =  2.000 #  Pm
        radi_vdw[59] =  2.000 #  Pr
        radi_vdw[62] =  2.000 #  Sm
        radi_vdw[65] =  2.000 #  Tb
        radi_vdw[69] =  2.000 #  Tm
        radi_vdw[70] =  2.000 #  Yb
        radi_vdw[87] =  2.000 #  Fr
        radi_vdw[88] =  2.000 #  Ra
        radi_vdw[89] =  2.000 #  Ac
        radi_vdw[95] =  2.000 #  Am
        radi_vdw[97] =  2.000 #  Bk
        radi_vdw[98] =  2.000 #  Cf
        radi_vdw[96] =  2.000 #  Cm
        radi_vdw[99] =  2.000 #  Es
        radi_vdw[93] =  2.000 #  Np
        radi_vdw[91] =  2.000 #  Pa
        radi_vdw[94] =  2.000 #  Pu
        radi_vdw[90] =  2.000 #  Th
        radi_vdw[92] =  1.860 #  U
        radi_vdw[100:self.maxgid] =  2.000

        for el_indx in range(1,self.maxgid):
            symbol_i  = elsymbol[el_indx]
            mass_i  = ELMASS[el_indx]
            cov_radii_i  = radi_cov[el_indx]
            vdw_radii_i  = radi_vdw[el_indx]
            el_i = element(symbol_i,el_indx,mass_i,cov_radii_i,vdw_radii_i)
            self.elements[symbol_i] = el_i

    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.elements
        del self.maxgid

    def getelementWithMass(self,mass_i):
        """
        Find element based on mass
        """
        mass_i_int = int(mass_i)
        for el_symb in self.elements:
            el_mass_int = int( self.elements[el_symb].mass )
            if( mass_i_int == el_mass_int ):
                el =  self.elements[el_symb]

        return el

    
    def getelementWithSymbol(self,symbol_i):
        """
        Find element based on mass
        """
        for el_symb in self.elements:
            if( el_symb == symbol_i ):
                el =  self.elements[el_symb]
                
        return el

    
