SRApy
=====

A python library and scripts to make working with NCBI's SRA less arcane.


Tools
=====


`get-project-sras.py`
----------------------

Script to download ALL sra run files for a given BioProject.

### Usage:

```
USAGE:
    get-project-sras.py [-e EMAIL -d OUTDIR] -p PROJECT_ID

OPTIONS:
    -e EMAIL        Your email, to provide to Bio.Entrez
                    [default: '']
    -d OUTDIR       Output directory, must exist. [default: .]
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

    get-projects-sras.py -d /path/to/sras -e you@example.com -p 30811

This will fetch the project, and search SRA for the metadata and run
accessions. It will download the SRA files, naming them by their SRA accession
and the submitter's sample label.
