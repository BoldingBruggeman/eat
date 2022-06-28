#!/usr/bin/env python

import os.path
import fnmatch
import shutil
import argparse
from typing import Iterable, Tuple

import numpy.random
import netCDF4
import yaml


# Hack into yaml parser to preserve order of yaml nodes,
# represent NULL by emoty string, skip interpretation of on/off as Boolean
import collections
def dict_representer(dumper, data):
    return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items())
def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))
def none_representer(self, _):
    return self.represent_scalar('tag:yaml.org,2002:null', '')
yaml_loader = yaml.SafeLoader
yaml_dumper = yaml.SafeDumper
del yaml_loader.yaml_implicit_resolvers['o']
del yaml_loader.yaml_implicit_resolvers['O']
yaml.add_representer(collections.OrderedDict, dict_representer, Dumper=yaml_dumper)
yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor, Loader=yaml_loader)
yaml.add_representer(type(None), none_representer, Dumper=yaml_dumper)


def perturb_restart(
    path,
    n: int,
    variable: str = "*",
    postfix: str = "_%04i",
    sigma: float = 0.1,
    exclude: Iterable[str] = (),
):
    name, ext = os.path.splitext(path)
    outpaths = [name + postfix % (i + 1) + ext for i in range(n)]
    for outpath in outpaths:
        print("Writing %s..." % outpath)
        shutil.copyfile(path, outpath)
        with netCDF4.Dataset(outpath, "r+") as nc:
            for name in fnmatch.filter(nc.variables.keys(), variable):
                ncvar = nc.variables[name]
                if name in exclude or (name in nc.dimensions and ncvar.ndim == 1):
                    continue
                values = ncvar[...]
                scale_factor = numpy.random.lognormal(sigma=sigma, size=values.shape)
                ncvar[...] = scale_factor * values


def perturb_yaml(
    path, n: int, params: Iterable[Tuple[str, float]], postfix: str = "_%04i"
):
    with open(path) as f:
        info = yaml.load(f, yaml_loader)
    par2value = {}
    for par, _ in params:
        root = info
        for comp in par.split("/"):
            assert comp in root, "Parameter %s not found in %s" % (par, path)
            root = root[comp]
        par2value[par] = float(root)

    name, ext = os.path.splitext(path)
    outpaths = [name + postfix % (i + 1) + ext for i in range(n)]

    for outpath in outpaths:
        print("Writing %s..." % outpath)
        for par, sigma in params:
            value = numpy.random.lognormal(sigma=float(sigma)) * par2value[par]
            root = info
            comps = par.split("/")
            for comp in comps[:-1]:
                root = root[comp]
            print("  %s: %s" % (par, value))
            root[comps[-1]] = value
        with open(outpath, "w") as f:
            yaml.dump(info, f, yaml_dumper)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd")

    parser_restart = subparsers.add_parser("restart")
    parser_restart.add_argument("file", help="path to restart file")
    parser_restart.add_argument("N", type=int, help="ensemble size")
    parser_restart.add_argument(
        "--sigma",
        type=float,
        help=(
            "standard deviation of ln scale factor for log-normally distributed"
            " perturbations"
        ),
        default=0.1,
    )
    parser_restart.add_argument(
        "-e",
        "--exclude",
        action="append",
        help="variable to exclude from perturbation",
        default=[],
    )

    parser_yaml = subparsers.add_parser("yaml")
    parser_yaml.add_argument("file", help="path to yaml file")
    parser_yaml.add_argument("N", type=int, help="ensemble size")
    parser_yaml.add_argument(
        "-p",
        "--param",
        nargs=2,
        action="append",
        help=(
            "parameter to perturb: path to parameter in yaml file, standard deviation"
            " of ln scale factor for log-normally distributed perturbations"
        ),
        default=[],
    )

    args = parser.parse_args()
    if args.cmd == "restart":
        perturb_restart(args.file, args.N, sigma=args.sigma, exclude=args.exclude)
    else:
        perturb_yaml(args.file, args.N, args.param)


if __name__ == "__main__":
    main()
