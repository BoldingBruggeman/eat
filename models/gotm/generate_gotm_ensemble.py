#!/usr/bin/env python

import os.path
import fnmatch
import shutil
import argparse

import numpy.random
import netCDF4

def perturb_restart(path, n, variable='*', postfix='_%04i', sigma=0.2):
    name, ext = os.path.splitext(path)
    outpaths = [name + postfix % (i + 1) + ext for i in range(n)]
    for outpath in outpaths:
        shutil.copyfile(path, outpath)
        with netCDF4.Dataset(outpath, 'r+') as nc:
            for name in fnmatch.filter(nc.variables.keys(), variable):
                ncvar = nc.variables[name]
                values = ncvar[...]
                scale_factor = numpy.random.lognormal(sigma=sigma, size=values.shape)
                ncvar[...] = scale_factor * values

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('restart', help='path to restart file')
    parser.add_argument('N', type=int, help='ensemble size')
    args = parser.parse_args()

    perturb_restart(args.restart, args.N)