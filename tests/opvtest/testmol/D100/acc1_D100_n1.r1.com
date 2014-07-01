%chk=acc1_D100_n1.chk
#P opt b3lyp/6-31g(d) geom=check guess=check pop=full density=current

acc1_D100_n1 ground state

0 1

--Link1--
%chk=acc1_D100_n1.chk
#P b3lyp td=nstates=12/6-31g(d) geom=check guess=check  pop=full density=current

acc1_D100_n1 TD-DFT

0 1

