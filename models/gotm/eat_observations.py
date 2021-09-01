from mpi4py import MPI
import numpy
import time

MPI_COMM_obs = MPI.COMM_WORLD.Split(color=2)
MPI_COMM_model = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)
MPI.COMM_WORLD.Barrier()

with open('obs_times.dat') as f:
    for i, l in enumerate(f):
        obs_time = l.rstrip('\n').strip('\'')
        MPI_COMM_obs.Send([obs_time.encode('ascii'), MPI.CHARACTER], dest=0, tag=1)
        nobs = 10000*(i + 1)
        MPI_COMM_obs.Send(numpy.array(nobs, dtype='i4'), dest=0, tag=1)
        dat = numpy.random.random(nobs)
        MPI_COMM_obs.Isend(dat, dest=0, tag=1).wait()
        time.sleep(1)

MPI_COMM_obs.Send([b'0000-00-00 00:00:00', MPI.CHARACTER], dest=0, tag=1)