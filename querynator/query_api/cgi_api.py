"""Query the cancergenomeinterpreter (CGI) via it's Web API"""


import json
import logging
import os
import os.path
import sys
import time
from datetime import date
from zipfile import BadZipfile, ZipFile

import click
import httplib2 as http
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def hg_assembly(genome):
    """
    Use correct assembly name

    :param genome: Genome build version, defaults to hg38
    :type genome: str
    :return: genome
    :rtype: str

    """

    if genome == "GRCh37":
        genome = "hg19"
    if genome == "GRCh38":
        genome = "hg38"
    return genome


def submit_query_cgi(mutations, cnas, translocations, genome, cancer, headers, logger):
    """
    Function that submits the query to the REST API of CGI

    :param mutations: Variant file (vcf,tsv,gtf,hgvs)
    :type mutations: str
    :param cnas: File with copy number alterations
    :type cnas: str
    :param translocations: File with translocations
    :type translocations: str
    :param genome: CGI takes hg19 or hg38
    :type genome: str
    :param email:  To query cgi a user account is needed
    :type email: str
    :param cancer: Cancer type from cancertypes.js
    :type cancer: str
    :param token: user token for CGI
    :type token: str
    :param logger: prints info to console
    :return: API url with job_id
    :rtype: str

    """

    logger.info("Querying REST API")

    genome = hg_assembly(genome)

    payload = {"cancer_type": cancer.name, "title": "CGI_query", "reference": genome}

    files = {"mutations": mutations, "cnas": cnas, "translocations": translocations}

    input_files = {k: open(v, "rb") for k, v in files.items() if v is not None}

    try:
        # submit query
        r = requests.post(
            "https://www.cancergenomeinterpreter.org/api/v1",
            headers=headers,
            files=input_files,
            data=payload,
        )
        r.raise_for_status()

        job_id = r.json()
        url = "https://www.cancergenomeinterpreter.org/api/v1/" + job_id
        return url

    except requests.exceptions.HTTPError as err:
        raise SystemExit(err)


def status_done(url, headers, logger):
    """
    Check query status

    :param url: API url with job_id
    :type url: str
    :param headers: Valid headers for API query
    :type headers: dict
    :raises: HTTPError
    :return: True if query performed successfully
    :rtype: bool

    """

    counter = 0
    payload = {"action": "logs"}

    retry_strategy = Retry(
        total=3, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)

    r = requests.get(url, headers=headers, params=payload, timeout=5)

    log = r.json()
    while "Analysis done" not in "".join(log["logs"]):
        if log["status"] == "Error":
            logger.exception("An Error has occurred with your request. Please check your input format")
            raise SystemExit()
        if counter == 5:
            print("Query took too looong :-(")
            break
        time.sleep(60)
        try:
            r = requests.get(url, headers=headers, params=payload, timeout=5)
            r.raise_for_status()
            counter += 1
            log = r.json()
        except requests.exceptions.RequestException as err:
            logger.exception("An Error has occurred with your request. Please check your input format")
            raise SystemExit(err)
        except requests.exceptions.HTTPError as err:
            raise SystemExit(err)
        except requests.exceptions.ConnectionError as err:
            logger.exception("Please check your internet connection.")
            raise SystemExit(err)
        except Exception as err:
            logger.exception("An unexpected error has occured " + type(err))
            raise SystemExit(err)

    return True


def download_cgi(url, headers, output, logger):
    """
    Download query results from cgi

    :param url: API url with job_id
    :type url: str
    :param headers: Valid headers for API query
    :type headers: dict
    :param output: sample name
    :type output: str
    :raises: Exception

    """

    # download results cgi
    try:
        payload = {"action": "download"}
        r = requests.get(url, headers=headers, params=payload)
        r.raise_for_status()
        with open(output + ".cgi_results.zip", "wb") as fd:
            fd.write(r._content)
    except requests.exceptions.HTTPError as err:
        raise SystemExit(err)
    except Exception:
        logger.exception("Ooops, sth went wrong with the download. Sorry for the inconvenience.")


def add_cgi_metadata(url, output):
    """
    Attach metadata to cgi query

    :param url: API url with job_id
    :type url: str
    :param output: sample name
    :type output: str
    :return: None
    :raises: BadZipfile

    """

    try:
        ZipFile(output + ".cgi_results.zip").extractall(output + ".cgi_results")
        # create additional file with metadata
        with open(output + ".cgi_results" + "/metadata.txt", "w") as f:
            f.write("CGI query date: " + str(date.today()))
            f.write("\nAPI version: " + url[:-20])
            f.close()
    except BadZipfile:
        logger.exception("Oops, sth went wrong with the zip archive. Please check your input format.")


def query_cgi(mutations, cnas, translocations, genome, cancer, headers, logger, output):
    """
    Actual query to cgi

    :param mutations: Variant file (vcf,tsv,gtf,hgvs)
    :type mutations: str
    :param cnas: File with copy number alterations
    :type cnas: str
    :param translocations: File with translocations
    :type translocations: str
    :param genome: Genome build version
    :type genome: str
    :param email:  To query cgi a user account is needed
    :type email: str
    :param cancer: Cancer type from cancertypes.js
    :type cancer: str
    :param logger: prints info to console
    :param output: sample name
    :type output: str

    """

    url = submit_query_cgi(mutations, cnas, translocations, genome, cancer, headers, logger)
    done = status_done(url, headers, logger)
    logger.info("CGI Query finished")
    if done:
        logger.info("Downloading CGI results")
        download_cgi(url, headers, output, logger)
        add_cgi_metadata(url, output)
