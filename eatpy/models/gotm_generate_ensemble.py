#!/usr/bin/env python

import os.path
import argparse

import numpy.random

from .gotm import YAMLEnsemble, RestartEnsemble


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="cmd", required=True)

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
    parser_yaml.add_argument(
        "-f",
        "--file",
        action="append",
        help=("parameter containing a file path to perturb"),
        dest="files",
        default=[],
    )

    args = parser.parse_args()
    if args.cmd == "restart":
        with RestartEnsemble(args.file, args.N) as f:
            for name, ncvar in f.template.items():
                if name in args.exclude or name in ncvar.dimensions:
                    continue
                shape = (f.n,) + ncvar.shape
                scale_factor = numpy.random.lognormal(sigma=args.sigma, size=shape)
                f[name] = scale_factor * f[name]
    else:
        with YAMLEnsemble(args.file, args.N) as f:
            for par, sigma in args.param:
                scale_factor = numpy.random.lognormal(sigma=float(sigma), size=f.n)
                f[par] = scale_factor * f.template[par]
            for par in args.files:
                filename, ext = os.path.splitext(f.template[par])
                f[par] = [filename + f.postfix % (i + 1) + ext for i in range(f.n)]


if __name__ == "__main__":
    main()
