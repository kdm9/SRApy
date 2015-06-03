#!/usr/bin/env python
from __future__ import print_function
from Bio import Entrez
from docopt import docopt
from os import path
import sys

import srapy
from srapy import (
    get_sample_runs,
    download_run,
)


CLI_USAGE = """
USAGE:
    get-project-sras.py [-e EMAIL -d OUTDIR -F FMT] -p PROJECT_ID

OPTIONS:
    -e EMAIL        Your email, to provide to Bio.Entrez
                    [default: ] # defaults to empty string
    -d OUTDIR       Output directory, must exist. [default: .]
    -F FMT          Filename format. Fields 'name', 'id', and 'acc' are
                    recognised. Use python string formatting syntax.
                    [default: {acc}~{name}.sra]
    -p PROJECT_ID   BioProject ID
"""


def main(argv=sys.argv[1:]):
    opts = docopt(CLI_USAGE, argv=argv)
    proj_id = int(opts['-p'])
    Entrez.email = opts['-e']
    outdir = opts['-d']
    fname_fmt = opts['-F']

    if not path.isdir(outdir):
        print("ERROR: output directory '{}' doesn't exitst".format(outdir),
              file=sys.stderr)
        exit(1)

    runs = get_sample_runs(proj_id)
    # add a line between getting run list & downloading
    print(file=sys.stderr)
    for idx, sra_id in enumerate(runs):
        print("Downloading run", idx + 1, "of", len(runs), file=sys.stderr)
        download_run(sra_id, outdir=outdir, fmt=fname_fmt)
        print(file=sys.stderr)  # Extra newline


if __name__ == "__main__":
    print("SRApy version", srapy.__version__, file=sys.stderr, end='\n\n')
    main()
