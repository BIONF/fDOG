# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2026 Vinh Tran
#
#  This file is part of fDOG tool https://github.com/BIONF/fDOG
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: tran@bio.uni-frankfurt.de
#
#######################################################################


import argparse
import os

from importlib.metadata import version
from fdog.libs.isoform_filter import filter_isoforms

def main(argv=None):
    fdog_version = version("fdog")

    parser = argparse.ArgumentParser(
        description=(
            "You are running fDOG version "
            + str(fdog_version)
            + ". Filter output files to retain only "
            "the best-scoring isoforms selected from a pp file."
        ),
        epilog=(
            "For more information on certain options, "
            "please refer to the wiki pages on github: "
            "https://github.com/BIONF/fDOG/wiki"
        )
    )

    parser.add_argument(
        "-p", "--pp",
        required=True,
        help="PhyloProfile file"
    )

    parser.add_argument(
        "-s", "--fasta",
        required=False, default=None,
        help="FASTA file containing all protein sequences"
    )

    parser.add_argument(
        "-f", "--forward",
        required=False, default=None,
        help="Forward domain annotation file"
    )

    parser.add_argument(
        "-r", "--reverse",
        required=False, default=None,
        help="Reverse domain annotation file"
    )

    parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Output directory. "
            "If not specified, input files are overwritten."
        )
    )

    args = parser.parse_args(argv)

    if args.output:
        os.makedirs(
            args.output,
            exist_ok=True
        )

    kept_seq = filter_isoforms(
        args.pp,
        args.fasta,
        args.forward,
        args.reverse,
        args.output
    )
    print(f"Done. {len(kept_seq)} orthologs retained.")

if __name__ == "__main__":
    main()