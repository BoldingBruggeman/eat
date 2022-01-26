#!/usr/bin/env python

import os.path
import fnmatch
import shutil
import argparse
from typing import Iterable

import numpy.random
import netCDF4

def perturb_restart(path, n: int, variable: str='*', postfix: str='_%04i', sigma: float=0.1, exclude: Iterable[str]=()):
    name, ext = os.path.splitext(path)
    outpaths = [name + postfix % (i + 1) + ext for i in range(n)]
    for outpath in outpaths:
        shutil.copyfile(path, outpath)
        with netCDF4.Dataset(outpath, 'r+') as nc:
            for name in fnmatch.filter(nc.variables.keys(), variable):
                ncvar = nc.variables[name]
                if name in exclude or (name in nc.dimensions and ncvar.ndim == 1):
                    continue
                values = ncvar[...]
                scale_factor = numpy.random.lognormal(sigma=sigma, size=values.shape)
                ncvar[...] = scale_factor * values

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('restart', help='path to restart file')
    parser.add_argument('N', type=int, help='ensemble size')
    parser.add_argument('--sigma', type=float, help='standard deviation of ln scale factor for log-normally distributed perturbations', default=0.1)
    parser.add_argument('-e', '--exclude', action='append', help='variable to exclude from perturbation', default=[])
    args = parser.parse_args()

    perturb_restart(args.restart, args.N, sigma=args.sigma, exclude=args.exclude)

if __name__ == '__main__':
    main()