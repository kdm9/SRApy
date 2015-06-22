from __future__ import print_function
from Bio import Entrez
from lxml import etree
from os import path
from progressbar import (
    Bar,
    Counter,
    ETA,
    FileTransferSpeed,
    Percentage,
    ProgressBar,
)
from sys import stderr
from six.moves import urllib


# SRA XML xpath
XPATHS = {
    'name':
        '//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/SAMPLE/@alias',
    'accession':
        '//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN/@accession'
}


def human_readable_size(size, suffix='B'):
    '''Returns a human-readable representation of the size of something'''
    for unit in ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z']:
        if abs(size) < 1000.0:
            return "{sz:3.3f}{unt}{suf}".format(sz=size, unt=unit, suf=suffix)
        size /= 1000.0
    return "{sz:3.3f}{unt}{suf}".format(sz=size, unt='Y', suf=suffix)


def urlretrieve(url, filename, stream=stderr, retries=3):
    '''Downloads ``url`` to path ``filename``, silently if 'stream' is None'''

    class HumanReadableCounter(Counter):
        '''Counter() which prints a human readable size'''
        def update(self, pbar):
            return human_readable_size(pbar.currval)

    if stream is not None:
        print('Downloading ', filename, file=stream, end=', ')

    # Establish connection to server
    exc = None
    urlhandle = None
    while retries > 0 and urlhandle is None:
        try:
            urlhandle = urllib.request.urlopen(url, timeout=10)
        except urllib.error.URLError as exc:
            retries -= 1
            urlhandle = None
            exc = exc
    if urlhandle is None:
        if exc is not None:
            raise exc
        else:
            raise RuntimeError("Download of", url, "failed")

    # Get metadata, check if we need to download
    meta = urlhandle.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    downloaded_size = 0
    # Check file size
    if path.exists(filename):
        downloaded_size = path.getsize(filename)
    if downloaded_size == file_size:
        if stream is not None:
            print('already downloaded!', file=stream)
        return

    # Set up the progress bar etc.
    if stream is not None:
        # On previous line
        print(human_readable_size(file_size), file=stream)
        # Widgets for the ProgressBar below
        widgets = [
            HumanReadableCounter(), ' ',
            FileTransferSpeed(), ' ',
            Bar(left='[', right=']'), ' ',
            Percentage(), ' ',
            ETA()
        ]
        pbar = ProgressBar(widgets=widgets, maxval=file_size)
        pbar.start()

    # Download the file
    file_size_dl = 0
    block_sz = 262144  # 256K
    with open(filename, 'wb') as fh:
        while True:
            buffer = urlhandle.read(block_sz)
            if not buffer:
                break
            fh.write(buffer)
            file_size_dl += len(buffer)
            if stream is not None:
                pbar.update(file_size_dl)
    if stream is not None:
        pbar.finish()


def accession_to_id(accession, force=False):
    '''Converts ``accession`` to being an id'''
    try:
        # Force the Entrez lookup by raising ValueError
        if force:
            raise ValueError
        # Try converting directly to an int, falling back to getting ID from
        # Entrez
        return int(accession)
    except ValueError:
        # Lookup numeric ID from accession
        term = '{}[Accession]'.format(accession)
        ids = esearch_ids(db='sra', term=term)
        if len(ids) == 1:
            return ids[0]
        else:
            raise ValueError("Couldn't get SRA UID of accession", accession)


def download_run(sra_id, outdir='.', stream=stderr, fmt='{acc}~{name}.sra'):
    '''Downloads run with SRA run uid ``sra_id`` to ``outdir``'''
    if stream is not None:
        print('Retrieving run info for run', sra_id, end='... ', file=stream)
    run_info = parse_run(sra_id)
    if stream is not None:
        print('done.', file=stream)
    url_template = ('ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/'
                    'reads/ByRun/sra/{leading3}/{leading6}/{all}/{all}.sra')
    run = run_info['accession']
    run_url = url_template.format(leading3=run[:3], leading6=run[:6], all=run)
    outfile = fmt.format(acc=run, name=run_info['name'], id=sra_id)
    outpath = path.join(outdir, outfile)
    urlretrieve(run_url, outpath)


def get_sample_runs(proj_id, stream=stderr):
    '''Gets a list of SRA run IDs from a given BioProject ID'''
    if stream is not None:
        print('Finding SRA runs for project', proj_id, end='... ', file=stream)
    term = '{}[BioProject]'.format(proj_id)
    sra_id_list = esearch_ids(db='sra', term=term)
    if stream is not None:
        print('found', len(sra_id_list), file=stream)
    return sra_id_list


def parse_run(sra_id):
    """
    Given an SRA id, gets the human-readable name and SRA Accession as a dict:

    {'id': sra_id, 'name': "sample_name", 'accession': 'SRR12345678'}

    sra_id should be a valid sra ID, as an int
    """
    sra_id = str(sra_id)
    handle = Entrez.efetch(db='sra', id=sra_id)
    xml_tree = etree.parse(handle)
    name = xml_tree.xpath(XPATHS['name'])[0]
    accession = xml_tree.xpath(XPATHS['accession'])[0]
    return {'id': sra_id, 'name': name, 'accession': accession}


def esearch_ids(**kwargs):
    '''
    Gets a list of IDs from an esearch query, handling mulitple requests.
    kwargs are passed directly to Bio.Entrez.esearch(), so consult their
    documentation, but 'db' and 'term' are required.
    '''
    handle = Entrez.esearch(**kwargs)
    rec = Entrez.read(handle)
    count = int(rec['Count'])
    id_list = rec['IdList']
    done = len(id_list)  # Num done so far
    while done < count:
        handle = Entrez.esearch(retstart=done, **kwargs)
        rec = Entrez.read(handle)
        id_list.extend(rec['IdList'])
        done += len(rec['IdList'])
    # Convert str-encoded ints to actual ints.
    id_list = map(int, id_list)
    return id_list


# versioneer
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
