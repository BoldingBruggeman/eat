#!/usr/bin/env python

import datetime
import argparse
import re
import logging
from typing import Optional, List, Tuple

import numpy

# Regular expression for ISO 8601 datetimes
datetimere = re.compile(r'(\d\d\d\d).(\d\d).(\d\d) (\d\d).(\d\d).(\d\d)\s*')

class File:
   def __init__(self, path: str, is_1d: bool=False, sd: Optional[float]=None, offset: int=0):
      self.f = open(path, 'r')
      self.is_1d = is_1d
      self.iline = 0
      self.next = ''
      self.sd = sd
      self.offset = offset
      self.read()

   def read(self):
      if self.next is None:
         # EOF reached previously
         return

      line = None
      while line is None:
         line = self.f.readline()
         if line == '':
            # EOF
            self.next = None
            return
         self.iline += 1
         line = line.split('#', 1)[0].strip()
         if line == '':
            # skip commented lines
            line = None

      datematch = datetimere.match(line)
      if datematch is None:
         raise Exception('Line %i does not start with time (yyyy-mm-dd hh:mm:ss). Line contents: %s' % (self.iline, line))
      refvals = map(int, datematch.group(1, 2, 3, 4, 5, 6)) # Convert matched strings into integers
      curtime = datetime.datetime(*refvals)
      data = line[datematch.end():].rstrip('\n').split()
      if self.is_1d:
         # Depth-explicit variable (on each line: time, depth, value)
         if len(data) not in (2, 3):
            raise Exception('Line %i must contain two values (depth, observation) or three values (depth, observation, sd) after the date + time, but it contains %i values.' % (self.iline, len(data)))
         if len(data) == 2 and self.sd is None:
            raise Exception('Line %i must contain three values (depth, observation, sd) after the date + time, since the standard deviation of observations has not been prescribed separately. However, the line contains %i values.' % (self.iline, len(data)))
         z = float(data[0])
         if not numpy.isfinite(z):
            raise Exception('Depth on line %i is not a valid number: %s.' % (self.iline, data[0]))
         value = float(data[1])
         sd = self.sd if self.sd is not None else float(data[2])
         self.next = curtime, z, value, sd
      else:
         # Depth-independent variable (on each line: time, value)
         if len(data) not in (1, 2):
            raise Exception('Line %i must contain one value (observation) or two values (observation, sd) after the date + time, but it contains %i values.' % (self.iline, len(data)))
         value = float(data[0])
         sd = self.sd if self.sd is not None else float(data[1])
         self.next = curtime, None, value, sd

class ObservationHandler:
   def __init__(self, start: Optional[datetime.datetime]=None, stop: Optional[datetime.datetime]=None, verbose: bool=False, logger: Optional[logging.Logger]=None):
      self.logger = logger or logging.getLogger('obs')
      if verbose:
         self.logger.setLevel(logging.DEBUG)

      self.start_time = start
      self.stop_time = stop
      self.time = None

      self.datasets: List[Tuple[str, File]] = []

   def add_observations(self, variable: str, path: str, is_1d: bool, sd: Optional[float]=None):
      offset = 0
      if '[' in variable and variable[-1] == ']':
         variable, depth_index = variable[:-1].split('[')
         offset = int(depth_index)
      self.datasets.append((variable, File(path, is_1d=is_1d, sd=sd, offset=offset)))

   def observations(self):
      assert self.datasets
      while True:
         # Find time of next observation
         self.time = None
         for _, obsfile in self.datasets:
            if obsfile.next is not None:
               current_time = obsfile.next[0]
               if self.time is None or current_time < self.time:
                  self.time = current_time

         if self.time is None:
            # No more observations
            break

         if (self.start_time is None or self.time >= self.start_time) and (self.stop_time is None or self.time <= self.stop_time):
            self.logger.debug('next observation time: %s' % self.time.strftime('%Y-%m-%d %H:%M:%S'))
            yield self.time
         else:
            self.logger.debug('skipping observations at time %s' % self.time.strftime('%Y-%m-%d %H:%M:%S'))
            for _, obsfile in self.datasets:
               while obsfile.next is not None and obsfile.next[0] == self.time:
                  obsfile.read()

   def collect_observations(self):
      self.depth_map = []
      self.values = []
      self.sds = []
      self.depth_indices = []
      istart = 0
      for variable, obsfile in self.datasets:
         zs = []
         while obsfile.next is not None and obsfile.next[0] == self.time:
            _, z, value, sd = obsfile.next
            self.logger.debug('- %s%s = %s (sd = %s)' % (variable, '' if z is None else (' @ %.2f m' % z), value, sd))
            zs.append(z if obsfile.is_1d else obsfile.offset)
            self.values.append(value)
            self.sds.append(sd)
            obsfile.read()
         self.depth_map.append((variable, istart, numpy.array(zs), obsfile.is_1d))
         istart = len(self.values)
      self.values = numpy.array(self.values, dtype=float)
      self.sds = numpy.array(self.sds, dtype=float)
      self.depth_indices = numpy.full(self.values.shape, -10000, dtype='i4')

   def get_observations(self, model_variables, model_states, variables):
      z_model = model_states[0, model_variables['z']['start']:model_variables['z']['start'] + model_variables['z']['length']]
      for model_variable, istart, zs, is_1d in self.depth_map:
         assert model_variable in variables, '%s not present in model state (after processing by plugins, if any)' % model_variable
         if is_1d:
            # zs contains depth coordinate (depth above mean sea level, typically negative)
            zs = numpy.abs(z_model[numpy.newaxis, :] - zs[:, numpy.newaxis]).argmin(axis=1)
         else:
            # zs contains depth index, which may be negative to indicate distance from surface (e.g., z=-1 for surface)
            zs = zs % variables[model_variable]['length']
         istop = istart + len(zs)
         self.depth_indices[istart:istop] = variables[model_variable]['start'] + zs
         self.logger.debug('observed %s: %s (sd %s) @ %s' % (model_variable, self.values[istart:istop], self.sds[istart:istop], self.depth_indices[istart:istop])) 
      return self.depth_indices, self.values, self.sds

   def add_arguments(self, parser: argparse.ArgumentParser):
      parser.add_argument('-o', '--obs', nargs=3, action='append', default=[])
      parser.add_argument('--start')
      parser.add_argument('--stop')

   def process_arguments(self, args):
      assert args.obs, 'At least one dataset with observations must be provided with -o/--obs'
      if args.start is not None:
         self.start_time = datetime.datetime.strptime(args.start, '%Y-%m-%d %H:%M:%S')
      if args.stop is not None:
         self.stop_time = datetime.datetime.strptime(args.stop, '%Y-%m-%d %H:%M:%S')
      for type, variable, path in args.obs:
         assert type in ('0', '1'), 'First argument after --obs must be the observation type: 0 for scalar, 1 for depth-explicit'
         type = int(type)
         self.add_observations(variable, path, is_1d=type==1)
