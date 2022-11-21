"""Command line interface querynator."""
import json
import logging
import os
import sys
from enum import Enum

import click
import requests

import querynator
from querynator.query_api import query_cgi

# Create logger
logger = logging.getLogger("Querynator")
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


class EnumType(click.Choice):
    """
    This is a class for a click.Choice of type EnumType

    """

    def __init__(self, enum, case_sensitive=False):
        self.__enum = enum
        super().__init__(choices=[item.value for item in enum], case_sensitive=case_sensitive)

    def convert(self, value, param, ctx):
        converted_str = super().convert(value, param, ctx)
        return self.__enum(converted_str)


def make_enum(values):
    """
    Function to create an EnumType from a dict

    :param values: json/dict like object with {key: value} pairs
    :type values: dict
    :return: enumeration
    :rtype: Enum

    """

    _k = _v = None

    class CancerType(Enum):
        nonlocal _k, _v
        for _k, _v in values.items():
            locals()[_k] = _v

    return CancerType


def Cancer():
    """
    Function to create instance of click.Choice EnumType with cancer types\n
    source: https://www.cancergenomeinterpreter.org/js/cancertypes.js

    :return: Enumeration of cancer types
    :rtype: click.Choice EnumType

    """

    directory_path = os.path.dirname(os.path.abspath(__file__))
    new_path = os.path.join(directory_path, "query_api/cancertypes.js")
    with open(new_path) as dataFile:
        data = dataFile.read()
        obj = data[data.find(" {") : data.rfind("};") + 1]
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

    logger.info("Start")


# querynator cgi_api
@querynator_cli.command()
@click.option(
    "-m",
    "--mutations",
    help="Please provide the path to the PASS filtered variant file (vcf,tsv,gtf,hgvs):\nSee more info here: https://www.cancergenomeinterpreter.org/faq#q22",
    type=click.Path(exists=True),
)
@click.option(
    "-a",
    "--cnas",
    help="Please provide the path to the copy number alterations file:\n",
    type=click.Path(exists=True),
)
@click.option(
    "-l",
    "--translocations",
    help="Please provide the path to the file containing translocations:\n",
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.STRING,
    help="Output name for output files - i.e. sample name. Extension filled automatically",
)
@click.option(
    "-c", "--cancer", help="Please enter the cancer type to be searched", type=EnumType(Cancer()), show_default=False
)
@click.option(
    "-g",
    "--genome",
    type=click.Choice(["hg19", "GRCh37", "hg38", "GRCh38"], case_sensitive=True),
    help="Please enter the genome version",
    required=True,
    default="hg38",
)
@click.option(
    "-t", "--token", help="Please provide your token for CGI database", required=True, type=click.STRING, default=None
)
@click.option(
    "-e",
    "--email",
    help="Please provide your user email address for CGI",
    required=False,
    type=click.STRING,
    default=None,
)
def query_api_cgi(mutations, cnas, translocations, cancer, genome, token, email, output):
    """
    Command to query the cancergenomeinterpreter API

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

    if mutations is None and cnas is None and translocations is None:
        raise click.UsageError(
            "No input file provided. Please provide at least one of [mutations/cnas/translocations] as input."
        )

    try:
        logger.info("Query the cancergenomeinterpreter (CGI)")
        headers = {"Authorization": email + " " + token}
        query_cgi(mutations, cnas, translocations, genome, cancer, headers, logger, output)

    except FileNotFoundError:
        print("Cannot find file on disk. Please try another path.")


if __name__ == "__main__":
    run_querynator()
