"""Query the cancergenomeinterpreter (CGI) via it's Web API"""


import glob
import logging
import os
import os.path
import sys
import time
from datetime import date

import click
import pandas as pd

import json
from zipfile import ZipFile

import httplib2 as http
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def hg_assembly(genome):
    """
    Use correct assembly name

    :param genome
    """

    if genome == "GRCh37":
        genome = "hg19"
    if genome == "GRCh38":
        genome = "hg38"
    return genome

def prepare_cgi_query_file(input, output):
    """
    Prepare a minimal file for upload
    VCF/TSV like structure with: CHROM POS ID REF ALT ...
    Attach sample name to get an output with a sample name
    """

    cgi_file = pd.read_csv(input, sep="\t")
    cgi_file = cgi_file.iloc[:,:5].drop_duplicates()
    cgi_file['SAMPLE'] = [output] * len(cgi_file.iloc[:,0])
    cgi_file.to_csv(output + '.cgi_input.tsv', sep="\t", index=False)
    cgi_path = output + '.cgi_input.tsv'
    return cgi_path



def submit_query_cgi(input, genome, cancer, headers, logger):
    """
    Function that submits the query to the REST API of CGI

    :param input: Input cgi file
    :param genome: CGI takes hg19 or hg38 str
    :param email:  To query cgi a user account is needed str
    :param cancer_type: str
    :param token: user token for CGI str
    :param logger: prints info to console
    """

    logger.info('Querying REST API')

    genome = hg_assembly(genome)

    payload = {'cancer_type': cancer.name, 'title': 'CGI_query', 'reference': genome}

    try:
        #TODO: ADD CNAS + TRANSLOCATIONS
        # submit query
        r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                    headers=headers,
                    files={ 'mutations': open(input, 'rb') },
                    data=payload)
        r.raise_for_status()

        job_id = r.json()
        url = 'https://www.cancergenomeinterpreter.org/api/v1/' + job_id
        return url

    except requests.exceptions.HTTPError as err:
        raise SystemExit(err)


def status_done(url, headers):
    """
    Check query status

    :param url:
    :param headers:
    :raises Exception
    :return
    """

    counter = 0
    payload={'action':'logs'}

    retry_strategy = Retry(
        total=3,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)

    r = requests.get(url, headers=headers, params=payload, timeout=5)

    log = r.json()
    while 'Analysis done' not in ''.join(log['logs']):
        if  counter == 5:
            print("Query took too looong :-(")
            break
        time.sleep(60)
        try:
            r = requests.get(url, headers=headers, params=payload, timeout=5)
            r.raise_for_status()
            counter +=1
            log = r.json()
        except requests.exceptions.HTTPError as err:
            raise SystemExit(err)

    return True


def download_cgi(url, headers, output):
    """
    Download query results from cgi
    :param url: str with job_id
    :param headers: str email + token
    :param output: str
    :raises Exception
    """

    # download results cgi
    try:
        payload={'action':'download'}
        r = requests.get(url, headers=headers, params=payload)
        with open(output + '.cgi_results.zip', 'wb') as fd:
            fd.write(r._content)
    except Exception:
        print("Ooops, sth went wrong with the download. Sorry for the inconvenience.")


def add_cgi_metadata(url, output):
    """
    Attach metadata to cgi query
    :param url:
    :param output:
    """

    ZipFile(output + '.cgi_results.zip').extractall(output + '.cgi_results')

    # create additional file with metadata
    with open(output + '.cgi_results' + '/metadata.txt', 'w') as f:
        f.write('CGI query date: ' +  str(date.today()))
        f.write('\nAPI version: ' + url[:-20])
        f.close()


def query_cgi(input, genome, cancer, headers, logger, output):
    """
    Actual query to cgi

    :param cgi_file: file
    :param genome: str
    :param cancer: str
    :param headers: str
    :param logger:
    :param output: str

    :return
    :raises Exception
    """
    cgi_file = prepare_cgi_query_file(input, output)
    url = submit_query_cgi(cgi_file, genome, cancer, headers, logger)
    done = status_done(url, headers)
    logger.info("CGI Query finished")
    if done:
        logger.info("Downloading CGI results")
        download_cgi(url, headers, output)
        add_cgi_metadata(url, output)
