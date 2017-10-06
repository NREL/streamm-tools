%chk=<tag>.chk
#P opt=Cartesian iop(1/7=12000) AM1 scf=qc

<tag>_AM1OPT

<charge> <spin_mult>
<coord>


--Link1--
%chk=<tag>.chk
#P opt <method>/<basis>  geom=(check,NewDefinition) guess=check  pop=full density=current

<tag> ground state

0 1


--Link1--
%chk=<tag>.chk
#P <method> td=nstates=<nstates>/<basis> geom=check guess=check  pop=full density=current

<tag> TD-DFT

0 1

