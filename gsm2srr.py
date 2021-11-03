from bs4 import BeautifulSoup
import requests
import time


def get_uid(gsm):
    r = requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
        params={
            'db': 'sra',
            'retmode': 'json',
            'term': gsm
        })
    return r.json()['esearchresult']['idlist']


def get_record(*uid):
    r = requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
        params={
            'db': 'sra',
            'id': ','.join(uid)
        })
    return r.text


def get_meta(record):
    soup = BeautifulSoup(record)
    for run in soup.findall('run'):
        srr = run['accession']
        srx = run.find('experiment_ref')['accession']
        gsm = run.find('experiment_ref')['refname']
        yield gsm, srx, srr


def main():
    gsm = 'GSM5214839'
    uid = get_uid(gsm)
    record = get_record(*uid)
    with open('test.xml', 'w') as f:
        f.write(record)
    time.sleep(1)


if __name__ == '__main__':
    main()
