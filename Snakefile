# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: Snakefile
# time: 2021/06/11
from bs4 import BeautifulSoup
from pathlib import Path
from io import StringIO
import pandas as pd
import requests
import shutil
import json
import time
import os
import re

CWD = os.getcwd()


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


rule fetch_one_gsm:
    output:
        CWD + '/gsm2srr/{gsm}.xml'
    resources:
        newwork=1
    run:
        uid = get_uid(wildcards.gsm)
        record = get_record(*uid)
        with open(output[0], 'w') as f:
            f.write(record)
        time.sleep(1)


def get_dump(record):
    with open(record, 'r') as f:
        soup = BeautifulSoup(f, 'lxml')
    try:
        gse = soup.find('study_ref')['refname']
    except TypeError:
        return
    for run in soup.find_all('run'):
        srr = run['accession']
        srx = run.find('experiment_ref')['accession']
        gsm = run.find('experiment_ref')['refname']
        srafile = run.find('srafile', semantic_name="run")
        md5 = srafile['md5']
        url = srafile.find('alternatives', org='AWS')['url']
        yield gse, gsm, srr, url, md5


def gsm_file(wildcards):
    with open(f'{CWD}/gsm.list') as f:
        return [
            f'{CWD}/gsm2srr/{gsm}.xml' for gsm in f.read().strip().splitlines()
        ]


checkpoint gsm2srr:
    output:
        CWD + '/gsm2srr.txt'
    input:
        gsm_file
    run:
        with open(output[0], 'w') as f:
            f.write('gsm\tsrx\tsrr\n')
            for file in input:
                for record in get_meta(file):
                    f.write('{}\t{}\t{}\n'.format(*record))


rule fetch_one_srr:
    output:
        CWD + '/filereport/{srr}.json'
    resources:
        network=1
    shell:
        'curl \'https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession={wildcards.srr}&format=json&fields=study_alias,sample_alias,run_accession,library_strategy,library_layout,fastq_ftp,fastq_aspera,fastq_bytes,fastq_md5,sra_ftp\''
        ' -o {output}'


def srr_file(wildcards):
    with checkpoints.gsm2srr.get().output[0].open() as f:
        gsm2srr = pd.read_table(f)
        return [
            f'{CWD}/filereport/{srr}.json' for srr in gsm2srr['srr']
        ]


def srr_records(files):
    for file in files:
        with open(file) as f:
            data = json.load(f)
        for record in data:
            yield record


checkpoint filereport:
    output:
        CWD + '/filereport.txt',
        CWD + '/download.list'
    input:
        srrlist=rules.gsm2srr.output,
        srrfiles=srr_file
    run:
        gsm2srr = pd.read_table('gsm2srr.txt')

        records = list(srr_records(input.srrfiles))
        records = pd.DataFrame.from_records(list(records))
        records = gsm2srr.merge(records, how='left', left_on='srr', right_on='run_accession')
        records.to_csv(output[0], sep='\t', index=False)

        files = records['fastq_ftp'].str.split(';').apply(pd.Series).stack().rename('ftp').reset_index()
        files['ascp'] = records['fastq_aspera'].str.split(';').apply(pd.Series).stack().rename('ascp').values
        files['md5'] = records['fastq_md5'].str.split(';').apply(pd.Series).stack().rename('md5').values
        files = files.merge(records[['gsm', 'srx', 'srr']], left_on='level_0', right_index=True)
        files['basename'] = files['ftp'].apply(os.path.basename)
        files['local'] = f'{CWD}/fastq/' + files['gse'] + '/' + files['gsm'] + '/' + files['srr'] + '/' + files['basename']
        files[['local', 'ftp', 'ascp', 'md5']].to_csv(output[1], sep='\t', index=False)


def fastq_ascp(wildcards):
    files = pd.read_table(f'{CWD}/download.list', index_col=0)
    fastq = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}_{r}.fastq.gz'.format(CWD=CWD, **wildcards)
    return files.loc[fastq, 'ascp']


def fastq_ftp(wildcards):
    files = pd.read_table(f'{CWD}/download.list', index_col=0)
    fastq = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}_{r}.fastq.gz'.format(CWD=CWD, **wildcards)
    return files.loc[fastq, 'ftp']


def fastq_md5(wildcards):
    files = pd.read_table(f'{CWD}/download.list', index_col=0)
    fastq = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}_{r}.fastq.gz'.format(CWD=CWD, **wildcards)
    return files.loc[fastq, 'md5']


# rule download_one_fastq:
#     output:
#         CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}_{r}.fastq.gz'
#     input:
#         ancient(rules.filereport.output)
#     resources:
#         network=1
#     run:
#         md5 = fastq_md5(wildcards)
#         downloader = config.get('downloader', 'wget')
#         if downloader == 'ascp':
#             url = fastq_ascp(wildcards)
#             ascp = shutil.which('ascp')
#             if ascp is None:
#                 raise ValueError('`ascp` is not installed, install it by `conda install -c hcc aspera-cli` or using ftp by add arg `--config downloader=ftp`')
#             key = Path(ascp).parents[1] / 'etc' / 'asperaweb_id_dsa.openssh'
#             download = f'ascp -l 300m -P 33001 -i {key} era-fasp@{url} {output[0]}.tmp'
#         elif downloader == 'wget':
#             url = fastq_ftp(wildcards)
#             download = f'wget -c -w 30 -T 60 -O {output[0]}.tmp ftp://{url}'
#         else:
#             raise ValueError(f'unknow downloader {downloader}')
#         download += (
#             f' && echo "{md5}  {output[0]}.tmp" > {output[0]}.md5'
#             f' && md5sum -c {output[0]}.md5'
#             f' && mv {output[0]}.tmp {output[0]}'
#             f' && rm {output[0]}.md5'
#             f' || sleep 60'
#         )
#         shell(download)


def fastq_file(wildcards):
    with checkpoints.filereport.get().output[1].open() as f:
        return [line.split('\t')[0] for line in f.read().strip().splitlines()[1:]]


rule download_all_fastq:
    input:
        rules.filereport.output,
        fastq_file


def get_dump(record):
    with open(record, 'r') as f:
        soup = BeautifulSoup(f, 'lxml')
    try:
        gse = soup.find('study_ref')['refname']
    except TypeError:
        return
    for run in soup.find_all('run'):
        srr = run['accession']
        srx = run.find('experiment_ref')['accession']
        gsm = run.find('experiment_ref')['refname']
        srafile = run.find('srafile', semantic_name="run")
        md5 = srafile['md5']
        url = srafile.find('alternatives', org='AWS')['url']
        yield gse, gsm, srr, url, md5


checkpoint dump_list:
    output:
        CWD + '/dump.list'
    input:
        gsm_file
    run:
        buffer = StringIO()
        buffer.write('gse\tgsm\tsrr\turl\tmd5\n')
        for file in input:
            for record in get_dump(file):
                buffer.write('{}\t{}\t{}\t{}\t{}\n'.format(*record))
        buffer.seek(0)
        files = pd.read_table(buffer)
        files['local'] = f'{CWD}/fastq/' + files['gse'] + '/' + files['gsm'] + '/' + files['srr'] +  '/' + files['srr'] + '.sra'
        files = files.set_index('local')
        files.to_csv(output[0], sep='\t')


def sra_url(wildcards):
    files = pd.read_table(f'{CWD}/dump.list', index_col=0)
    sra = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}.sra'.format(CWD=CWD, **wildcards)
    return files.loc[sra, 'url']


def sra_md5(wildcards):
    files = pd.read_table(f'{CWD}/dump.list', index_col=0)
    sra = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}.sra'.format(CWD=CWD, **wildcards)
    return files.loc[sra, 'md5']
        


rule download_one_sra:
    output:
        temp(CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}.sra')
    input:
        ancient(rules.dump_list.output)
    resources:
        network=1
    run:
        url = sra_url(wildcards)
        md5 = sra_md5(wildcards)
        download = (
            f'wget -c {url} -O {output[0]}.tmp'
            f' && echo "{md5}  {output[0]}.tmp" > {output[0]}.md5'
            f' && md5sum -c {output[0]}.md5'
            f' && mv {output[0]}.tmp {output[0]}'
            f' && rm {output[0]}.md5'
            f' || sleep 60'
        )
        shell(download)


rule dump_one_fastq:
    output:
        CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}_1.fastq.gz',
        CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}_2.fastq.gz'
    input:
        CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}.sra'
    priority:
        10
    params:
        outdir=CWD + '/fastq/{gse}/{gsm}/{srr}'
    shell:
        'fastq-dump --split-files --gzip -O {params.outdir} {input}'


def sra_files(wildcards):
    with checkpoints.dump_list.get().output[0].open() as f:
        return [line.split('\t')[0].strip('.sra') + '_1.fastq.gz' for line in f.read().strip().splitlines()[1:]]

rule dump_all_fastq:
    input:
        rules.dump_list.output,
        sra_files
