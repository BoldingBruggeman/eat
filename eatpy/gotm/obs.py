#!/usr/bin/env python

import datetime
import argparse
import re
import logging
from typing import Optional

import numpy
from mpi4py import MPI

from .. import shared

# Regular expression for ISO 8601 datetimes
datetimere = re.compile(r'(\d\d\d\d).(\d\d).(\d\d) (\d\d).(\d\d).(\d\d)\s*')

class File:
   def __init__(self, path: str, is_1d: bool=False):
      self.f = open(path, 'r')
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
            yield curtime, z, float(data[1])
         else:
            # Depth-independent variable (on each line: time, value)
            #if len(data) != 1:
            #   raise Exception('Line %i does not contain one value (observation) after the date + time, but %i values.' % (iline + 1, len(data)))
            yield curtime, float(data[0])

class ObservationHandler:
   def __init__(self, state_layout_path: str='da_variables.dat', start: Optional[datetime.datetime]=None, stop: Optional[datetime.datetime]=None, logger: Optional[logging.Logger]=None):
      self.logger = logger or logging.getLogger('obs')

      self.start_time = start
      self.stop_time = stop

      comm, self.MPI_COMM_obs_model, self.MPI_COMM_obs_filter, _ = shared.setup_mpi(shared.COLOR_OBS)
      assert comm.size == 1, 'There can only be one instance of the observation handler active.'
      assert self.MPI_COMM_obs_model.rank == 0
      assert self.MPI_COMM_obs_filter.rank == 0
      self.have_filter = self.MPI_COMM_obs_filter.size > 1

      size_obs_model = self.MPI_COMM_obs_model.size
      self.nmodel = size_obs_model - comm.size
      self.logger.info(' obs() connected with {} models'.format(self.nmodel))

      # Wait for model to generate the file thta describes the memory layout
      self.MPI_COMM_obs_model.Barrier()

      self.logger.info('Parsing memory map %s' % state_layout_path)
      self.memory_map = {}
      for name, metadata in shared.parse_memory_map(state_layout_path):
         self.memory_map[name] = (metadata['start'] + 1, metadata['start'] + metadata['length'])

      self.datasets = []

   def add_observations(self, variable: str, path: str):
      idx = None
      if '[' in variable and variable[-1] == ']':
         variable, idx = variable[:-1].split('[')
         idx = int(idx)
      assert variable in self.memory_map
      istart, istop = self.memory_map[variable]
      if idx is None:
         assert istart == istop, 'Variable %s can only be assimilated at a single depth. Specify a depth index like this: %s[<INDEX>], using for instance <INDEX>=0 for bottom or <INDEX>=-1 for surface.' % (variable, variable)
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
         if (self.start_time is None or time >= self.start_time) and (self.stop_time is None or time <= self.stop_time):
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
            self.logger.info('(-> model) {}'.format(strtime))
            for dest in range(1, self.nmodel + 1):
               reqs.append(self.MPI_COMM_obs_model.Issend([strtime.encode('ascii'), MPI.CHARACTER], dest=dest, tag=shared.TAG_TIMESTR))

         if self.have_filter:
            # Send new observations to filter
            reqs.append(self.MPI_COMM_obs_filter.Issend([strtime.encode('ascii'), MPI.CHARACTER], dest=1, tag=shared.TAG_TIMESTR))
            nobs = iobs.size
            self.logger.info('(-> filter) {}'.format(nobs))
            reqs.append(self.MPI_COMM_obs_filter.Issend(numpy.array(nobs, dtype='i4'), dest=1, tag=shared.TAG_NOBS))
            if nobs > 0:
               reqs.append(self.MPI_COMM_obs_filter.Issend(iobs, dest=1, tag=shared.TAG_IOBS))
               reqs.append(self.MPI_COMM_obs_filter.Issend(obs, dest=1, tag=shared.TAG_OBS))

      if self.have_filter:
         self.MPI_COMM_obs_filter.Send([b'0000-00-00 00:00:00', MPI.CHARACTER], dest=1, tag=shared.TAG_TIMESTR)
         self.MPI_COMM_obs_filter.Send(numpy.array(-1, dtype='i4'), dest=1, tag=shared.TAG_NOBS)
      
      if self.nmodel:
         for dest in range(1, self.nmodel + 1):
            self.MPI_COMM_obs_model.Send([b'0000-00-00 00:00:00', MPI.CHARACTER], dest=dest, tag=shared.TAG_TIMESTR)

def main():
   logging.basicConfig(level=logging.INFO)

   parser = argparse.ArgumentParser()
   parser.add_argument('-o', '--obs', nargs=2, action='append', default=[])
   parser.add_argument('--start')
   parser.add_argument('--stop')
   args = parser.parse_args()
   assert args.obs, 'At least one dataset with observations must be provided with -o/--obs'
   if args.start is not None:
      args.start = datetime.datetime.strptime(args.start, '%Y-%m-%d %H:%M:%S')
   if args.stop is not None:
      args.stop = datetime.datetime.strptime(args.stop, '%Y-%m-%d %H:%M:%S')
   handler = ObservationHandler(start=args.start, stop=args.stop)
   for variable, path in args.obs:
      handler.add_observations(variable, path)
   handler.start()

if __name__ == '__main__':
   main()