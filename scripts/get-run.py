#!/usr/bin/env python
from __future__ import print_function
from Bio import Entrez
from docopt import docopt
from os import path
import sys

import srapy
from srapy import (
    accession_to_id,
    download_run,
)


CLI_USAGE = """
USAGE:
    get-run.py [-e EMAIL -d OUTDIR -F FMT -a] (-i SRA_ID | -f FILE)

OPTIONS:
    -e EMAIL        Your email, to provide to Bio.Entrez
                    [default: ] # defaults to empty string
    -d OUTDIR       Output directory, must exist. [default: .]
    -a              The IDs are accessions. Unless '-a' is specified, the
                    identifier is assumed to be an SRA ID if it is numeric,
                    otherwise is interpreted as an accession. This option
                    forces the ID to be interpreted as an accession.
                    [default: False]
    -F FMT          Filename format. Fields 'name', 'id', and 'acc' are
                    recognised. Use python string formatting syntax.
                    [default: {acc}~{name}.sra]
    -i SRA_ID       A single identifier to download. See above for
                    interpretation
    -f FILE         New-line delimited list of identifiers
"""


def main(argv=sys.argv[1:]):
    opts = docopt(CLI_USAGE, argv=argv)
    Entrez.email = opts['-e']
    outdir = opts['-d']
    fname_fmt = opts['-F']

    if not path.isdir(outdir):
        print("ERROR: output directory '{}' doesn't exitst".format(outdir),
              file=sys.stderr)
        exit(1)

    ids = []

    if opts['-i']:
        ids.append(opts['-i'])
    else:
        with open(opts['-f']) as fh:
            for line in fh:
                ids.append(line.strip())

    uids = []
    bad_ids = []
    print("Getting accession UIDs ...", end='\r', file=sys.stderr)
    for idx, id_or_acc in enumerate(ids):
        try:
            uids.append(accession_to_id(id_or_acc, force=opts['-a']))
            print("Getting accession UIDs ...", idx + 1, "of", len(ids),
                  end='\r', file=sys.stderr)
        except ValueError:
            print("\nWARNING: Skipping non-existant accession", id_or_acc,
                  end='\n', file=sys.stderr)
            bad_ids.append(id_or_acc)
    print("\33[2KGetting accession UIDs ... Done", end='\n\n', file=sys.stderr)

    for idx, sra_id in enumerate(uids):
        print("Downloading run", idx + 1, "of", len(uids), file=sys.stderr)
        download_run(sra_id, outdir=outdir, fmt=fname_fmt)
        print(file=sys.stderr)  # Extra newline between records

    print("Finished!", file=sys.stderr)
    if len(bad_ids) > 0:
        print("The following accessions or IDs were skipped:", file=sys.stderr)
        for id in bad_ids:
            print("\t", id, file=sys.stderr)


if __name__ == "__main__":
    print("SRApy version", srapy.__version__, file=sys.stderr, end='\n\n')
    main()
