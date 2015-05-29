from Bio import Entrez
from lxml import etree
from docopt import docopt
import urllib2

CLI_USAGE = """
USAGE:
    get-project.py [-e EMAIL] -p PROJECT_ID

OPTIONS:
    -e EMAIL        Your email, to provide to Bio.Entrez
                    [default: '']
    -p PROJECT_ID   BioProject ID
"""

# SRA XML xpath
XPATHS = {
    'name':
        '//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/SAMPLE/@alias',
    'accession':
        '//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN/@accession'
}


def esearch_ids(**kwargs):
    handle = Entrez.esearch(**kwargs)
    rec = Entrez.read(handle)
    count = int(rec['Count'])
    id_list = rec['IdList']
    done = int(rec['RetMax'])  # Num done so far
    while done < count:
        handle = Entrez.esearch(retstart=done, **kwargs)
        rec = Entrez.read(handle)
        id_list.extend(rec['IdList'])
        done += int(rec['RetMax'])
    id_list = map(int, id_list)
    return id_list


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


def urlretrieve(url, filename):
    urlhandle = urllib2.urlopen(url)
    meta = urlhandle.info()
    file_size_dl = 0
    block_sz = 1048576  # 1MB
    file_size = int(meta.getheaders("Content-Length")[0])
    print 'Downloading {}, {}B'.format(filename, file_size)
    with open(filename, 'wb') as fh:
        while True:
            buffer = urlhandle.read(block_sz)
            if not buffer:
                break
            file_size_dl += len(buffer)
            fh.write(buffer)
            pct_complete = file_size_dl * 100. / file_size
            status = '{:3.2f}% Downloaded\r'.format(pct_complete)
            print status,
    print "Done"


def get_run(sra_id):
    run_info = parse_run(sra_id)
    url_template = ('ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/'
                    'reads/ByRun/sra/{leading3}/{leading6}/{all}/{all}.sra')
    run_url = url_template.format(leading3=run[:3], leading6=run[:6], all=run)
    outfile = '{}~{}.sra'.format(run_info['accession'], run_info['name'])
    urlretrieve(run_url, outfile)


def get_sample_runs(proj_id):
    term = '{}[BioProject]'.format(proj_id)
    sra_id_list = esearch_ids(db='sra', term=term)
    return sra_id_list


def main(proj_id, email):
    Entrez.email = email

    for sra_id in get_sample_runs(proj_id):
        get_run(sra_id)


if __name__ == "__main__":
    opts = docopt(CLI_USAGE)
    main(int(opts['-p']), opts['-e'])
