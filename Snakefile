# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: Snakefile
# time: 2021/06/11
from bs4 import BeautifulSoup
from pathlib import Path
import pandas as pd
import shutil
import json
import time
import os
import re

CWD = os.getcwd()


rule fetch_one_gsm:
    output:
        CWD + '/sra/{gsm}.xml'
    resources:
        newwork=1
    shell:
        'esearch -db sra -query {wildcards.gsm} | efetch -format native > {output} && sleep 10 || sleep 30'


def get_meta(record):
    if os.path.getsize(record) == 0:
        raise OSError(f'{record} is a empty file')
    with open(record, 'r') as f:
        soup = BeautifulSoup(f, 'lxml')
    try:
        gse = soup.find('study')['alias']
        paired = soup.find('library_layout').paired
    except TypeError:
        return
    for run in soup.find_all('run'):
        gsm = run.find('experiment_ref')['refname']
        srr = run['accession']
        srafile = run.find('srafile', semantic_name="run")
        url = sorted(
            (alt['org'], alt['url'])
            for alt in srafile.find_all('alternatives', access_type='anonymous')
        )[0][1]
        md5 = srafile['md5']
        sra = f'{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}.sra'
        if paired:
            fastq=f'{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}_1.fastq.gz,{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}_2.fastq.gz'
        else:
            fastq=f'{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}.fastq.gz'
        yield sra, gse, gsm, srr, url, md5, fastq


def gsm_files(wildcards):
    with open(f'{CWD}/gsm.list') as f:
        return [
            f'{CWD}/sra/{gsm}.xml' for gsm in f.read().strip().splitlines()
        ]


checkpoint metadata:
    output:
        CWD + '/metadata.list'
    input:
        gsm_files
    run:
        with open(output[0], 'w') as f:
            f.write('sra\tgse\tgsm\tsrr\turl\tmd5\tfastq\n')
            for file in input:
                for record in get_meta(file):
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(*record))


def sra_url(wildcards):
    files = pd.read_table(f'{CWD}/metadata.list', index_col=0)
    sra = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}.sra'.format(CWD=CWD, **wildcards)
    return files.loc[sra, 'url']


def sra_md5(wildcards):
    files = pd.read_table(f'{CWD}/metadata.list', index_col=0)
    sra = '{CWD}/fastq/{gse}/{gsm}/{srr}/{srr}.sra'.format(CWD=CWD, **wildcards)
    return files.loc[sra, 'md5']


rule download_one_sra:
    output:
        temp(CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}.sra')
    input:
        ancient(rules.metadata.output)
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


rule dump_paired_fastq:
    output:
        temp(CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}_1.fastq'),
        temp(CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}_2.fastq')
    input:
        CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}.sra'
    params:
        workdir=CWD + '/fastq/{gse}/{gsm}',
        outdir=CWD + '/fastq/{gse}/{gsm}/{srr}'
    shell:
        'cd {params.workdir}'
        ' && fasterq-dump -O {params.outdir} -t {params.outdir} -e 8 -S {wildcards.srr}'


rule dump_singal_fastq:
    output:
        temp(CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}.fastq')
    input:
        CWD + '/fastq/{gse}/{gsm}/{srr}/{srr}.sra'
    params:
        workdir=CWD + '/fastq/{gse}/{gsm}',
        outdir=CWD + '/fastq/{gse}/{gsm}/{srr}'
    shell:
        'cd {params.workdir}'
        ' && fasterq-dump -O {params.outdir} -t {params.outdir} -e 8 -S {wildcards.srr}'


rule comprese_fastq:
    output:
        CWD + '/fastq/{gse}/{gsm}/{srr}/{fastq}.fastq.gz'
    input:
        CWD + '/fastq/{gse}/{gsm}/{srr}/{fastq}.fastq'
    shell:
        'pigz -p 4 -c {input} > {output}.tmp && mv {output}.tmp {output}'


def fastq_files(wildcards):
    with checkpoints.metadata.get().output[0].open() as f:
        return sum([line.split('\t')[-1].split(',') for line in f.read().strip().splitlines()[1:]], [])


rule dump_all_fastq:
    input:
        rules.metadata.output,
        fastq_files
