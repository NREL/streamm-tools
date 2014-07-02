%chk=acc1_D100_n1.chk
#P opt=Cartesian iop(1/7=12000) AM1 scf=qc

acc1_D100_n1 semi-empirical pre-optimiz

0 1
C -1.554986 -1.911312 -0.000810
C -0.177760 -1.911312 -0.000810
C 0.347615 -0.579042 -0.000810
C -0.658845 0.361011 0.000000
S -2.169481 -0.356146 -0.000008
H -2.189661 -2.795265 -0.001321
H 0.453890 -2.801454 -0.001064
H 1.416824 -0.359618 -0.001382
H -0.519437 1.440247 0.000647


--Link1--
%chk=acc1_D100_n1.chk
#P opt b3lyp/6-31g(d)  geom=(check,NewDefinition) guess=check  pop=full density=current

acc1_D100_n1 ground state

0 1


--Link1--
%chk=acc1_D100_n1.chk
#P b3lyp td=nstates=12/6-31g(d) geom=check guess=check  pop=full density=current

acc1_D100_n1 TD-DFT

0 1

