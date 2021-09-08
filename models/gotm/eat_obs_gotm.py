from mpi4py import MPI
import numpy
import time

# For MPI_split
color_obs=1
color_model=2
color_filter=4

# Tags
tag_timestr=1
tag_nobs=2
tag_obs=3

# Setup communicators
MPI_COMM_obs = MPI.COMM_WORLD.Split(color=color_obs)
rank = MPI_COMM_obs.Get_rank()
size = MPI_COMM_obs.Get_size()
MPI_COMM_model = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)
MPI_COMM_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)

MPI_COMM_obs_model = MPI.COMM_WORLD.Split(color=color_obs+color_model,key=-1)
rank_obs_model = MPI_COMM_obs_model.Get_rank()
size_obs_model = MPI_COMM_obs_model.Get_size()
have_model = size_obs_model > 1
if not have_model:
   print(' obs(no model program present)')
   rank_obs_model = -1
   size_obs_model = -1
else:
   nmodel=size_obs_model-size

MPI_COMM_obs_filter = MPI.COMM_WORLD.Split(color=color_obs+color_filter,key=-1)
rank_obs_filter = MPI_COMM_obs_filter.Get_rank()
size_obs_filter = MPI_COMM_obs_filter.Get_size()
have_filter = size_obs_filter > 1
if not have_filter:
   print(' obs(no filter program present)')
   rank_obs_filter = -1
   size_obs_filter = -1

MPI_COMM_model_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)

# Communicators are created

print(' obs(ranks and sizes: O-OM-OF): %4d %4d %4d %4d %4d %4d' % (rank,rank_obs_model,rank_obs_filter,size,size_obs_model,size_obs_filter) )

MPI.COMM_WORLD.Barrier()

# Everything is initialized

with open('obs_times.dat') as f:
    for i, l in enumerate(f):
        obs_time = l.rstrip('\n').strip('\'')
        print(obs_time)
        if have_model:
           for dest in range(1, nmodel+1):
              MPI_COMM_obs_model.Send([obs_time.encode('ascii'), MPI.CHARACTER], dest=dest, tag=tag_timestr)
        if have_filter:
           nobs = 10000*(i + 1)
           MPI_COMM_obs_filter.Send(numpy.array(nobs, dtype='i4'), dest=1, tag=tag_nobs)
           dat = numpy.random.random(nobs)
           MPI_COMM_obs_filter.Isend(dat, dest=1, tag=tag_obs).wait()

if have_model:
   for dest in range(1, nmodel+1):
      MPI_COMM_obs_model.Send([b'0000-00-00 00:00:00', MPI.CHARACTER], dest=dest, tag=tag_timestr)

if have_filter:
   nobs=-1
   MPI_COMM_obs_filter.Send(numpy.array(nobs, dtype='i4'), dest=1, tag=tag_nobs)
   
