from os import path
from shutil import rmtree
from tempfile import mkdtemp

import srapy
from nose.tools import (
    assert_raises,
)

OUTDIR = None


def setup():
    '''Setup function for this module'''
    srapy.Entrez.email = 'kevin.murray@anu.edu.au'

    global OUTDIR
    OUTDIR = mkdtemp('srapy-test')


def teardown():
    '''Tear-down function for this module'''
    global OUTDIR
    rmtree(OUTDIR)
    OUTDIR = None


def test_accession_to_id():
    '''Test srapy.acccession_to_id'''

    # accession
    sra_id = srapy.accession_to_id('SRR1999343')
    assert sra_id == 1466194, sra_id

    # str(id)
    sra_id = srapy.accession_to_id('1466194')
    assert sra_id == 1466194, sra_id

    # int(id)
    sra_id = srapy.accession_to_id(1466194)
    assert sra_id == 1466194, sra_id

    # bad accession
    with assert_raises(ValueError):
        sra_id = srapy.accession_to_id('notacc')


def test_esearch_ids():
    '''Test srapy.esearch_ids'''

    # Search with fewer IDs than 'retmax'
    id_list = sorted(srapy.esearch_ids(db='bioproject', term='1001 genomes'))
    assert id_list == [30811, 197428, 273563], id_list

    # Search with more IDs than 'retmax'
    sra_id_list = [
        1508680, 1508681, 1508682, 1508683, 1508685, 1508686, 1508687, 1508688,
        1508689, 1508690, 1508691, 1508692, 1508693, 1508694, 1508695, 1508696,
    ]
    id_list = sorted(srapy.esearch_ids(db='sra', term='PRJDB3021'))
    assert id_list == sra_id_list, id_list

    # Search with no results
    id_list = srapy.esearch_ids(db='sra', term='abasdkjfslakdjflskdjflasjdf')
    assert id_list == [], id_list


def test_urlretrieve():
    '''Test srapy.urlretrieve'''

    # Test downloading a file silently
    outfile = path.join(OUTDIR, 'test_1mb.dat')
    srapy.urlretrieve('http://ftp.iinet.net.au/1mb.dat', outfile, stream=None)
    assert path.exists(outfile), 'Download failed'
