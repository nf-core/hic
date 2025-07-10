#!/usr/bin/env python

## nf-core-hic
## Copyright (c) 2020 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

"""
Script to merge any files with the same template
"""

import argparse
import sys
import glob
import os
from collections import OrderedDict


def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


if __name__ == "__main__":
    ## Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="List of input file(s)", type=str, nargs="+")
    parser.add_argument("-v", "--verbose", help="verbose mode", action="store_true")
    args = parser.parse_args()

    infiles = args.files
    li = len(infiles)

    if li > 0:
        if args.verbose:
            print("## merge_statfiles.py")
            print("## Merging " + str(li) + " files")

        ## Reading first file to get the template
        template = OrderedDict()
        if args.verbose:
            print("## Use " + infiles[0] + " as template")
        with open(infiles[0]) as f:
            for line in f:
                if not line.startswith("#"):
                    lsp = line.strip().split("\t")
                    data = map(num, lsp[1 : len(lsp)])
                    template[str(lsp[0])] = list(data)

        if len(template) == 0:
            print("Cannot find template files !")
            sys.exit(1)

        ## Int are counts / Float are percentage
        for fidx in list(range(1, li)):
            with open(infiles[fidx]) as f:
                for line in f:
                    if not line.startswith("#"):
                        lsp = line.strip().split("\t")
                        if lsp[0] in template:
                            for i in list(range(1, len(lsp))):
                                if isinstance(num(lsp[i]), int):
                                    template[lsp[0]][i - 1] += num(lsp[i])
                                else:
                                    template[lsp[0]][i - 1] = round((template[lsp[0]][i - 1] + num(lsp[i])) / 2, 3)
                        else:
                            sys.stderr.write(
                                "Warning : '" + lsp[0] + "' not found in template [" + infiles[fidx] + "]\n"
                            )

        ## Print template
        for x in template:
            sys.stdout.write(x)
            for y in template[x]:
                sys.stdout.write("\t" + str(y))
            sys.stdout.write("\n")

    else:
        print("No files to merge - stop")
        sys.exit(1)
