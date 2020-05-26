# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
fiberassign.scripts.merge
============================

High-level functions for merging target catalogs with output files.

"""
from __future__ import absolute_import, division, print_function

import os
import argparse

import numpy as np

from ..utils import Logger

from ..assign import (join_results, result_tiles)

def parse_join(optlist=None):
    """Parse joining options.

    This parses either sys.argv or a list of strings passed in.  If passing
    an option list, you can create that more easily using the
    :func:`option_list` function.

    Args:
        optlist (list, optional): Optional list of arguments to parse instead
            of using sys.argv.

    Returns:
        (namespace):  an ArgumentParser namespace.

    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--targets", type=str, required=True, nargs="+",
                        help="Input file with targets of any type.  This "
                        "argument can be specified multiple times (for "
                        "example if standards / skies / science targets are "
                        "in different files).")

    parser.add_argument("--sky", type=str, required=False, nargs="+",
                        help="Input file with sky or 'bad sky' targets.  "
                        "This option exists in order to treat main-survey"
                        " sky target files as valid for other survey types."
                        "  If you are running a main survey assignment, you"
                        " can just pass the sky file to the --targets list.")

    parser.add_argument("--dir", type=str, required=True, default=None,
                        help="Directory containing fiberassign results.")

    parser.add_argument("--prefix", type=str, required=False,
                        default="fba-",
                        help="Prefix of each file (before the <tile>.fits).")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Output directory for the joine files.  Default"
                        " is the directory containing the fiberassign output.")

    parser.add_argument("--out_prefix", type=str, required=False,
                        default="summary_fba_",
                        help="Prefix of the joined file.")




    args = None
    if optlist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(optlist)

    if args.sky is None:
        args.sky = list()

    return args

def run_join_init(args, comm=None):
    """Initialize joining inputs.

    This uses the previously parsed options to load the input files needed.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        (tuple):  The (Tiles, columns) needed to run the merging.

    """
    log = Logger.get()
    # Check directory
    if (comm is None) or (comm.rank == 0):
        if not os.path.isdir(args.dir):
            log.error("Results directory {} does not exist".format(args.dir))
            if comm is not None:
                comm.Abort()
    return 

def run_join(args):
    """Run output joining.

    This uses the previously parsed options to read input data and perform
    merging of the input catalogs.  This runs on one node.

    Args:
        args (namespace): The parsed arguments.

    Returns:
        None

    """
    tiles, columns = run_merge_init(args)
    join_results(args.targets, args.sky, result_dir=args.dir,
                  result_prefix=args.prefix, 
                  out_dir=args.out, out_prefix=args.out_prefix
    return
