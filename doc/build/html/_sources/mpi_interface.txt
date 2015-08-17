.. _mpi_interface:

.. index:: MPI, Boost, mpi4py, communication


MPI Interface Class
============================

The MPI module :mod:`mpiBase` provides a common interface for typical
communication operations. The implementation of these operations is
handled by derived classes that use a variety of external
libraries. Currently these include:

- boost
- mpi4py
- serial

where 'serial' is explicitly listed to emphasize that this implementation
uses the same interface as the actual MPI libraries. Therefore, a user
can write an algorithm using the MPI interface that does not have to
be changed if their code is moved to a system that does not have any
of the supported MPI libraries installed. Instructions and a sample
code for testing the MPI routines follows.

1.  Paste code below to a file --> e.g. 'commTst.py'
2.  While logged onto Peregrine load an MPI python module e.g. 'module load boost/1.54.0.a/impi-intel'
3.  Run 'mpirun -n 2 python commTst.py'


Example code: commTst.py::

    import mpiBase
    p = mpiBase.getMPIObject()

    # MPI setup
    rank = p.getRank()
    size = p.getCommSize()

    x = 0              # Sets x on all processors
    p.barrier()
    if rank == 0:      # Resets x on proc 0
        x = 12.3456
    p.barrier()

    # Check x setting on each processor
    for proc in range(size):
        if proc == rank:
            print "Before broadcast x = ", x, " on processor ", proc
        p.barrier()
    p.barrier()

    # Broadcast x to all processors (default root proc is '0')
    xAll = p.bcast(x)
    p.barrier()

    # Check x setting on each processor
    for proc in range(size):
        if proc == rank:
            print "After broadcast x = ", xAll, " on processor ", proc
        p.barrier()
    p.barrier()

Note that this file can be run if one loads another python MPI library
or none at all.
