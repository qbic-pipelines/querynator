"""Command line interface querynator."""
import click
import sys
import logging

import querynator
import querynator.cgi_api

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

    print("blub")


# querynator cgi_api
@querynator_cli.command()
@click.option('-i', '--input', help='Enter path to PASS filtered vcf/tsv:\n',
                required=True, type=click.Path(exists=True))
@click.option('-o', '--output', required=True, type=click.STRING,
                help='Output name for output files - i.e. sample name. Extension filled automatically')
@click.option('-c', '--cancer', help='Enter cancer type to be searched',
                type=click.Choice(['liver', 'LVB'], case_sensitive=False), required=False, default=None)
@click.option('-g', '--genome', type=click.Choice(['hg19', 'GRCh37', 'hg38', 'GRCh38'], case_sensitive=True),
                help='Please enter the genome version', required=True)
@click.option('-t', '--token', help='Enter token for CGI database', required=False,
                type=click.STRING, default=None)
def query_api_cgi(input, cancer, genome, token, output):
    """
    Command to query the cancergenomeinterpreter API
    """

    try:
        print("bla")
    except FileNotFoundError:
        print('Cannot find file on disk. Please try another path.')


if __name__ == "__main__":
    run_querynator()
