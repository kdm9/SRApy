from __future__ import print_function
from docopt import docopt
from srapy import (
    get_sample_runs,
    download_run,
)

CLI_USAGE = """
USAGE:
    get-project.py [-e EMAIL] -p PROJECT_ID

OPTIONS:
    -e EMAIL        Your email, to provide to Bio.Entrez
                    [default: '']
    -p PROJECT_ID   BioProject ID
"""


def main(proj_id, email):
    Entrez.email = email

    for sra_id in get_sample_runs(proj_id):
        download_run(sra_id)


if __name__ == "__main__":
    opts = docopt(CLI_USAGE)
    main(int(opts['-p']), opts['-e'])
