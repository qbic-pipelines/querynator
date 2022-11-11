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
from biothings_client import get_client
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

def prepare_cgi_query_file(infile):
    """
    Prepare a minimal file for upload
    Either hgvs or CHROM POS ID REF ALT
    Attach sample name to get an output with a sample name
    """

    file_name = os.path.basename(infile).split(".")[0]
    cgi_file = pd.read_csv(infile, sep="\t")
    cgi_file = cgi_file.iloc[:,:5].drop_duplicates()
    cgi_file['SAMPLE'] = [file_name] * len(cgi_file.iloc[:,0])
    cgi_file.head()
    cgi_file.to_csv(file_name + '.cgi_input.tsv', sep="\t", index=False)
    cgi_path = file_name + '.cgi_input.tsv'
    return cgi_path



def submit_query_cgi(infile, genome, cancer, headers, logger):
    """
    Function that submits the query to the REST API of CGI

    :param infile
    :param genome: CGI takes hg19 or hg38 str
    :param email:  To query cgi a user account is needed str
    :param cancer_type: str
    :param token: user token for CGI str
    :param logger: prints info to console
    """

    logger.info('Querying CGI directly')

    genome = hg_assembly(genome)

    #'susanne.jodoin@qbic.uni-tuebingen.de 4b1f20ee84976f38f1b6', 'LVB' (liver)

    payload = {'cancer_type': cancer, 'title': 'CGI_query', 'reference': genome}

    try:
        # submit query
        r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                    headers=headers,
                    files={ 'mutations': open(infile, 'rb') },
                    data=payload)

        job_id = r.json()
        url = 'https://www.cancergenomeinterpreter.org/api/v1/' + job_id
        return url

    except Exception:
        print("An exception occurred")


def status_done(url, headers):
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
            # print('status: ', log['logs'])
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


def query_cgi(infile, genome, cancer, headers, logger, output):
    """
    Actual query to cgi

    :param infile: file
    :param genome: str
    :param cancer: str
    :param headers: str
    :param logger:
    :param output: str

    :return
    :raises Exception
    """
    prepare_cgi_query_file(infile)
    url = submit_query_cgi(infile, genome, cancer, headers, logger)
    done = status_done(url, headers)
    logger.info("CGI Query finished")
    if done:
        logger.info("Downloading CGI results")
        download_cgi(url, headers, output)
        add_cgi_metadata(url, output)

def add_cgi_metadata(url, output):

    ZipFile(output + '.cgi_results.zip').extractall('cgi_results')

    # create additional file with metadata
    with open('cgi_results' + '/metadata.txt', 'w') as f:
        f.write('CGI query date: ' +  str(date.today()))
        f.write('\nAPI version: ' + url[:-20])
        f.close()
