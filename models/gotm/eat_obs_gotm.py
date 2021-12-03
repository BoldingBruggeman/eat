#!/usr/bin/env python

from mpi4py import MPI
import numpy
import io
import re
import datetime
import argparse

# For MPI_split - MUST match colors defined in eat_config.F90
color_obs=1
color_model=2
color_filter=4

# Tags - MUST match tags defined in eat_config.F90
tag_timestr=1
tag_nobs=1
tag_iobs=2
tag_obs=3
tag_analysis=1
tag_forecast=2

# Regular expression for ISO 8601 datetimes
datetimere = re.compile(r'(\d\d\d\d).(\d\d).(\d\d) (\d\d).(\d\d).(\d\d)\s*')

class File:
   def __init__(self, path, is_1d=False):
      self.f = io.open(path, 'r')
      self.is_1d = is_1d

   def read(self):
      for iline, line in enumerate(self.f):
         if line.startswith('#'):
            continue
         datematch = datetimere.match(line)
         if datematch is None:
            raise Exception('Line %i does not start with time (yyyy-mm-dd hh:mm:ss). Line contents: %s' % (iline + 1, line))
         refvals = map(int, datematch.group(1, 2, 3, 4, 5, 6)) # Convert matched strings into integers
         curtime = datetime.datetime(*refvals)
         data = line[datematch.end():].rstrip('\n').split()
         if self.is_1d:
            # Depth-explicit variable (on each line: time, depth, value)
            if len(data) != 2:
               raise Exception('Line %i does not contain two values (depth, observation) after the date + time, but %i values.' % (iline + 1, len(data)))
            z = float(data[0])
            if not numpy.isfinite(z):
               raise Exception('Depth on line %i is not a valid number: %s.' % (iline + 1, data[0]))
            yield curtime, z, float(data[2])
         else:
            # Depth-independent variable (on each line: time, value)
            #if len(data) != 1:
            #   raise Exception('Line %i does not contain one value (observation) after the date + time, but %i values.' % (iline + 1, len(data)))
            yield curtime, float(data[1])

class ObservationHandler:
   def __init__(self, state_layout_path='da_variables.dat'):
      # Communicator for observation handler only - verify there is only 1
      MPI_COMM_obs = MPI.COMM_WORLD.Split(color=color_obs)
      assert MPI_COMM_obs.Get_size() == 1, 'There can only be one instance of the observation handler active.'

      # Inter-model and inter-filter communicators - we do not participate in those, but still have to call split
      MPI_COMM_model = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)
      MPI_COMM_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)

      # Communicator between observation handler (1) and models (multiple)
      # The observation handler will get a new rank of 0 as we provide the lowest key
      self.MPI_COMM_obs_model = MPI.COMM_WORLD.Split(color=color_obs + color_model, key=-1)
      assert self.MPI_COMM_obs_model.Get_rank() == 0
      size_obs_model = self.MPI_COMM_obs_model.Get_size()
      self.nmodel = size_obs_model - 1
      print(' obs() connected with {} models'.format(self.nmodel))

      # Communicator between observation handler (1) and filter (1)
      # The observation handler will get a new rank of 0 as we provide the lowest key
      self.MPI_COMM_obs_filter = MPI.COMM_WORLD.Split(color=color_obs + color_filter, key=-1)
      assert self.MPI_COMM_obs_filter.Get_rank() == 0
      self.have_filter = self.MPI_COMM_obs_filter.Get_size() > 1

      # Communicator between filter and models - we do not participate in this one, but still have to call split
      MPI_COMM_model_filter = MPI.COMM_WORLD.Split(color=MPI.UNDEFINED)

      # Communicators are created

      MPI.COMM_WORLD.Barrier()

      # Everything is initialized

      self.MPI_COMM_obs_model.Barrier()

      print('Parsing memory map %s' % state_layout_path)
      self.memory_map = {}
      with io.open(state_layout_path, 'r') as f:
         for l in f:
            name, long_name, units, category, dims, start, length = l.split('\t')
            istart, istop = int(start), int(start) + int(length) - 1
            print('  %s @ %i:%i - %s (%s)' % (name, istart, istop, long_name, units))
            self.memory_map[name] = (istart, istop)

      self.datasets = []

   def add_observations(self, variable, path):
      idx = None
      if '[' in variable and variable[-1] == ']':
         variable, idx = variable[:-1].split('[')
         idx = int(idx)
      assert variable in self.memory_map
      istart, istop = self.memory_map[variable]
      if idx is None:
         assert istart == istop
         idx = istart
      elif idx < 0:
         idx = istop + idx + 1
      else:
         idx = istart + idx
      self.datasets.append((variable, idx, File(path)))

   def observations(self):
      assert self.datasets
      variable, idx, obsfile = self.datasets[0]
      for time, value in obsfile.read():
         yield time, numpy.array([idx], dtype='i4'), numpy.array([value], dtype=float)

   def start(self):
      reqs = []
      for obs_time, iobs, obs in self.observations():
         assert iobs.size == obs.size
         assert iobs.size != 0

         # Make sure all previous requests have been received
         MPI.Request.Waitall(reqs)

         reqs = []
         if self.nmodel:
            # Send new time to ensemble members
            strtime = obs_time.strftime('%Y-%m-%d %H:%M:%S')
            print(' obs(-> model) {}'.format(strtime))
            for dest in range(1, self.nmodel + 1):
               reqs.append(self.MPI_COMM_obs_model.Issend([strtime.encode('ascii'), MPI.CHARACTER], dest=dest, tag=tag_timestr))

         if self.have_filter:
            # Send new observations to filter
            nobs = iobs.size
            print(' obs(-> filter) {}'.format(nobs))
            reqs.append(self.MPI_COMM_obs_filter.Issend(numpy.array(nobs, dtype='i4'), dest=1, tag=tag_nobs))
            if nobs > 0:
               reqs.append(self.MPI_COMM_obs_filter.Issend(iobs, dest=1, tag=tag_iobs))
               reqs.append(self.MPI_COMM_obs_filter.Issend(obs, dest=1, tag=tag_obs))

      if self.have_filter:
         nobs = -1
         self.MPI_COMM_obs_filter.Send(numpy.array(nobs, dtype='i4'), dest=1, tag=tag_nobs)
      
      if self.nmodel:
         for dest in range(1, self.nmodel + 1):
            self.MPI_COMM_obs_model.Send([b'0000-00-00 00:00:00', MPI.CHARACTER], dest=dest, tag=tag_timestr)

if __name__ == '__main__':
   handler = ObservationHandler()
   parser = argparse.ArgumentParser()
   parser.add_argument('-o', '--obs', nargs=2, action='append', default=[])
   args = parser.parse_args()
   assert args.obs, 'At least one dataset with observations must be provided with -o/--obs'
   for variable, path in args.obs:
      handler.add_observations(variable, path)
   handler.start()