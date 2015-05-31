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
import urllib2


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


def urlretrieve(url, filename, silent=False):
    '''Downloads ``url`` to path ``filename``, silently if ``silent`` is True'''

    class HumanReadableCounter(Counter):
        '''Counter() which prints a human readable size'''
        def update(self, pbar):
            return human_readable_size(pbar.currval)

    if not silent:
        print('Downloading ', filename, file=stderr, end=', ')
    urlhandle = urllib2.urlopen(url)
    meta = urlhandle.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    downloaded_size = 0
    if path.exists(filename):
        downloaded_size = path.getsize(filename)
    if downloaded_size == file_size:
        if not silent:
            print('already downloaded!', file=stderr)
        return
    # Check file size
    if not silent:
        print(human_readable_size(file_size), file=stderr)
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
    file_size_dl = 0
    block_sz = 262144  # 256K
    with open(filename, 'wb') as fh:
        while True:
            buffer = urlhandle.read(block_sz)
            if not buffer:
                break
            fh.write(buffer)
            file_size_dl += len(buffer)
            if not silent:
                pbar.update(file_size_dl)
            pct_complete = file_size_dl * 100. / file_size
    if not silent:
        pbar.finish()


def download_run(sra_id, outdir='.', silent=False):
    '''Downloads run with SRA run uid ``sra_id`` to ``outdir``'''
    if not silent:
        print('Retrieving run info for run', sra_id, end='... ', file=stderr)
    run_info = parse_run(sra_id)
    if not silent:
        print('done.', file=stderr)
    url_template = ('ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/'
                    'reads/ByRun/sra/{leading3}/{leading6}/{all}/{all}.sra')
    run = run_info['accession']
    run_url = url_template.format(leading3=run[:3], leading6=run[:6], all=run)
    outfile = '{}~{}.sra'.format(run, run_info['name'])
    outpath = path.join(outdir, outfile)
    urlretrieve(run_url, outpath)


def get_sample_runs(proj_id, silent=False):
    '''Gets a list of SRA run IDs from a given BioProject ID'''
    if not silent:
        print('Finding SRA runs for project', proj_id, end='... ', file=stderr)
    term = '{}[BioProject]'.format(proj_id)
    sra_id_list = esearch_ids(db='sra', term=term)
    if not silent:
        print('found', len(sra_id_list), file=stderr)
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
