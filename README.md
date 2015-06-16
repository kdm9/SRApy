SRApy
=====

A python library and scripts to make working with NCBI's SRA less arcane.


Installation
============

    pip install srapy

is all anyone on a modern system with pip should need. Dependencies are:

- `BioPython`
- `lxml`
- `docopt`
- `progressbar2`

`lxml` depends on the libxml and libxslt C libraries, so you may need to
install them with `sudo aptitude install libxml2-dev libxslt1-dev` or similar.

License
=======

SRApy is licensed under the GNU General Public License, version 3 or any later
version. See `./LICENSE` or https://www.gnu.org/licenses/gpl-3.0.en.html


Tools
=====

`get-project-sras.py`
----------------------

Script to download ALL sra run files for a given BioProject.

### Usage:

```
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
```

### Example:

Search for a project on the BioProject search engine:
http://www.ncbi.nlm.nih.gov/bioproject
, e.g, the [1001 genomes project](http://1001genomes.org)
([search link](http://www.ncbi.nlm.nih.gov/bioproject/?term=1001+genomes)).

Copy the ID of the BioProject of interest, after manual curation, e.g. `30811`
for Joe Ecker's contribution to the 1001 genomes project.

To download all SRA files, run:

    get-project-sras.py -d /path/to/sras -e you@example.com -p 30811

This will fetch the project, and search SRA for the metadata and run
accessions. It will download the SRA files, naming them by their SRA accession
and the submitter's sample label.


`get-run.py`
------------

Given an ID or accession, or list thereof, download the run file.

### Usage:

```
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
```

### Example:

    echo ERR605369 >sra_accessions.txt
    echo ERR612613 >>sra_accessions.txt
    echo ERR612614 >>sra_accessions.txt

    # Download a single accession
    get-run.py -e you@example.com -i ERR605369

    # Download a single accession (ERR612613), by run ID
    get-run.py -e you@example.com -i 1018875

    # Download all accessions
    get-run.py -e you@example.com -f sra_accessions.txt
