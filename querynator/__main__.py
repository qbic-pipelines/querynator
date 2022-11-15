"""Command line interface querynator."""
import click
import sys
import logging
import requests
from enum import Enum
import json
import os

import querynator
import querynator.query_api.cgi.cgi_api as c

# Create logger
logger = logging.getLogger('Querynator')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


class EnumType(click.Choice):
    def __init__(self, enum, case_sensitive=False):
        self.__enum = enum
        super().__init__(choices=[item.value for item in enum], case_sensitive=case_sensitive)

    def convert(self, value, param, ctx):
        converted_str = super().convert(value, param, ctx)
        return self.__enum(converted_str)


def make_enum(values):
    _k = _v = None
    class CancerType(Enum):
        nonlocal _k, _v
        for _k, _v in values.items():
            locals()[_k] = _v
    return CancerType

def Cancer():
    """
    Function to create instance of Cancer Enum
    source of file: https://www.cancergenomeinterpreter.org/js/cancertypes.js
    """
    with open('./querynator/query_api/cgi/cancertypes.js') as dataFile:
        data = dataFile.read()
        obj = data[data.find(' {') : data.rfind('};')+1]
        jsonObj = json.loads(obj)
        Cancer_enum = make_enum(jsonObj)

    return Cancer_enum


def run_querynator():
    print("\n                                           __ ")
    print("  ____ ___  _____  _______  ______  ____ _/ /_____  _____")
    print(" / __ `/ / / / _ \/ ___/ / / / __ \/ __ `/ __/ __ \/ ___/")
    print("/ /_/ / /_/ /  __/ /  / /_/ / / / / /_/ / /_/ /_/ / /")
    print("\__, /\__,_/\___/_/   \__, /_/ /_/\__,_/\__/\____/_/")
    print("  /_/                /____/\n\n")

    # Launch the click cli
    querynator_cli()


@click.group()
@click.version_option(version=querynator.__version__)
def querynator_cli():
    """
    querynator is a command-line tool to query cancer variant databases.
    It can query the web APIs or local instances of the databases
    """

    logger.info('Start')


# querynator cgi_api
@querynator_cli.command()
@click.option('-i', '--input', help='Please provide the path to the PASS filtered vcf/tsv:\n',
                required=True, type=click.Path(exists=True))
@click.option('-o', '--output', required=True, type=click.STRING,
                help='Output name for output files - i.e. sample name. Extension filled automatically')
@click.option('-c', '--cancer', help='Please enter the cancer type to be searched',
                type=EnumType(Cancer()), show_default=False)
@click.option('-g', '--genome', type=click.Choice(['hg19', 'GRCh37', 'hg38', 'GRCh38'], case_sensitive=True),
                help='Please enter the genome version', required=True, default = 'hg38')
@click.option('-t', '--token', help='Please provide your token for CGI database', required=True,
                type=click.STRING, default=None)
@click.option('-e', '--email', help='Please provide your user email address for CGI', required=False, type=click.STRING, default=None)
def query_api_cgi(input, cancer, genome, token, email, output):
    """
    Command to query the cancergenomeinterpreter API
    """

    try:
        logger.info("Query the cancergenomeinterpreter (CGI)")
        headers = {'Authorization': email + ' ' + token }
        c.query_cgi(input, genome, cancer, headers, logger, output)

    except FileNotFoundError:
        print('Cannot find file on disk. Please try another path.')


if __name__ == "__main__":
    run_querynator()
