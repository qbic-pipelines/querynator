"""Command line interface querynator."""
import json
import logging
import os
import random
import shutil
from enum import Enum

import click
import vcf
from vcf.parser import _Info as VcfInfo
from vcf.parser import field_counts as vcf_field_counts

import querynator
from querynator.helper_functions import gunzip_compressed_files, gzipped
from querynator.query_api import query_cgi, query_civic, vcf_file
from querynator.report_scripts import (
    add_tiers_and_scores_to_df,
    combine_cgi,
    combine_cgi_civic,
    combine_civic,
    create_report_htmls,
)

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
logger.propagate = False


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


def filter_vcf_by_vep(vcf_path, logger):
    """
    Function to filter given vcf to remove synonymous and low impact variants based on VEP annotation

    :param vcf_path: Variant Call Format (VCF) file (Version 4.2)
    :type vcf_path: str
    :return: list of lists of pyVCF3 records (input file, removed, filtered)
    :rtype: list

    """

    if not vcf_file(vcf_path):
        logger.error("Can only filter variants in vcf files.")
        exit(1)

    if gzipped(vcf_path):
        vcf_path = gunzip_compressed_files(vcf_path, logger)

    # read vcf file in pyVCF
    in_vcf = vcf.Reader(open(vcf_path))

    # creates dictionary with VEP info names as keys and index in list as columns
    # Name must be VEPs default "CSQ"
    if "CSQ" in in_vcf.infos:
        logger.info("Filtering vcf file")

        vep_dict = {name: pos for pos, name in enumerate(in_vcf.infos["CSQ"].desc.split(":")[1].strip().split("|"))}
        """
        Exemplary for nf-core/sarek (https://nf-co.re/sarek) output
        {'Allele': 0,
        'Consequence': 1,
        'IMPACT': 2,
        'SYMBOL': 3,
        'Gene': 4,
        'Feature_type': 5,
        'Feature': 6,
        'BIOTYPE': 7,
        'EXON': 8,
        'INTRON': 9,
        'HGVSc': 10,
        'HGVSp': 11,
        'cDNA_position': 12,
        'CDS_position': 13,
        'Protein_position': 14,
        'Amino_acids': 15,
        'Codons': 16,
        'Existing_variation': 17,
        'DISTANCE': 18,
        'STRAND': 19,
        'FLAGS': 20,
        'VARIANT_CLASS': 21,
        'SYMBOL_SOURCE': 22,
        'HGNC_ID': 23,
        'CANONICAL': 24,
        'MANE_SELECT': 25,
        'MANE_PLUS_CLINICAL': 26,
        'TSL': 27,
        'APPRIS': 28,
        'CCDS': 29,
        'ENSP': 30,
        'SWISSPROT': 31,
        'TREMBL': 32,
        'UNIPARC': 33,
        'UNIPROT_ISOFORM': 34,
        'GENE_PHENO': 35,
        'SIFT': 36,
        'PolyPhen': 37,
        'DOMAINS': 38,
        'miRNA': 39,
        'AF': 40,
        'AFR_AF': 41,
        'AMR_AF': 42,
        'EAS_AF': 43,
        'EUR_AF': 44,
        'SAS_AF': 45,
        'AA_AF': 46,
        'EA_AF': 47,
        'gnomAD_AF': 48,
        'gnomAD_AFR_AF': 49,
        'gnomAD_AMR_AF': 50,
        'gnomAD_ASJ_AF': 51,
        'gnomAD_EAS_AF': 52,
        'gnomAD_FIN_AF': 53,
        'gnomAD_NFE_AF': 54,
        'gnomAD_OTH_AF': 55,
        'gnomAD_SAS_AF': 56,
        'MAX_AF': 57,
        'MAX_AF_POPS': 58,
        'FREQS': 59,
        'CLIN_SIG': 60,
        'SOMATIC': 61,
        'PHENO': 62,
        'PUBMED': 63,
        'MOTIF_NAME': 64,
        'MOTIF_POS': 65,
        'HIGH_INF_POS': 66,
        'MOTIF_SCORE_CHANGE': 67,
        'TRANSCRIPTION_FACTORS': 68}
        """

        to_remove = []
        to_keep = []
        for record in in_vcf:
            if len(record.INFO["CSQ"]) > 1:
                # filter out low impact & synonymous variants from variants with multiple vep annotations
                if all(
                    info_list[vep_dict["IMPACT"]] == "LOW"
                    and info_list[vep_dict["Consequence"]] == "synonymous_variant"
                    for info_list in [i.split("|") for i in record.INFO["CSQ"]]
                ):
                    to_remove.append(record)
                else:
                    to_keep.append(record)
            else:
                # filter out low impact & synonymous variants from variants with single vep annotations
                if (
                    record.INFO["CSQ"][0].split("|")[vep_dict["IMPACT"]] == "LOW"
                    and record.INFO["CSQ"][0].split("|")[vep_dict["Consequence"]] == "synonymous_variant"
                ):
                    to_remove.append(record)
                else:
                    to_keep.append(record)

        return [in_vcf, to_keep, to_remove]

    else:
        logger.error("vcf file does not include required VEP INFO fields (key must be default 'CSQ')")
        exit(1)


def write_vcf(vcf_template, vcf_record_list, out_name):
    """
    Function to write a vcf file from list of pyvcf3 records to result directory

    :param vcf_header: pysam header object from input vcf
    :type vcf_header: pysam header object
    :param vcf_record_list: list of pysam records
    :type vcf_record_list: list
    :param out_name: name for the created vcf file
    :type out_name: str
    :return: None
    :rtype: None
    """

    # add querynator_id info to header
    vcf_template.infos["QID"] = VcfInfo("QID", ".", "String", "Querynator ID", ".", ".", ".")
    writer = vcf.Writer(open(f"{out_name}", "w"), vcf_template, lineterminator="\n")

    for record in vcf_record_list:
        # add querynator_id to record
        record.add_info("QID", random.randint(1000000, 9999999))
        writer.write_record(record)


def get_unique_querynator_dir(querynator_output):
    """
    add index if "querynator_results" already exists in user given out dir

    :param querynator_output: path to store querynator results
    :type querynator_output: str
    :return: unique result directory
    :rtype: str
    """
    filename, extension = os.path.splitext(querynator_output)

    count = 1
    while os.path.exists(querynator_output):
        querynator_output = filename + f"_{count}"
        count += 1

    return querynator_output


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
    It can query web APIs or local instances of the databases
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
    "--outdir",
    required=True,
    type=click.STRING,
    help="i.e. sample name. Directory in which results will be stored.",
)
@click.option(
    "-c",
    "--cancer",
    help="Please enter the cancer type to be searched. You must use quotation marks.",
    type=EnumType(Cancer()),
    show_default=False,
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
@click.option(
    "-f",
    "--filter_vep",
    help="If set, filters out synoymous and low impact variants based on VEP annotation",
    is_flag=True,
    show_default=True,
    default=False,
)
def query_api_cgi(mutations, cnas, translocations, cancer, genome, token, email, outdir, filter_vep):
    if mutations is None and cnas is None and translocations is None:
        raise click.UsageError(
            "No input file provided. Please provide at least one of [mutations/cnas/translocations] as input."
        )

    try:
        result_dir = get_unique_querynator_dir(f"{outdir}")
        dirname, basename = os.path.split(result_dir)
        original_input = {"mutations": mutations, "translocations": translocations, "cnas": cnas}
        # filter vcf file if required
        if mutations is not None and filter_vep:
            in_vcf_header, candidate_variants, removed_variants = filter_vcf_by_vep(mutations, logger)

            # create result directories
            os.makedirs(f"{result_dir}/vcf_files")
            write_vcf(in_vcf_header, removed_variants, f"{result_dir}/vcf_files/{basename}.removed_variants.vcf")
            write_vcf(in_vcf_header, candidate_variants, f"{result_dir}/vcf_files/{basename}.filtered_variants.vcf")

            # create and set new input file for cgi query
            mutations = f"{result_dir}/vcf_files/{basename}.filtered_variants.vcf"

        logger.info("Query the cancergenomeinterpreter (CGI)")
        headers = {"Authorization": email + " " + token}
        # run analysis
        query_cgi(
            mutations, cnas, translocations, genome, cancer, headers, logger, result_dir, original_input, filter_vep
        )

        # move downloaded results to result dir
        if dirname != "":
            if not filter_vep:
                os.makedirs(result_dir)
            shutil.move(f"{dirname}/{basename}.cgi_results", result_dir)
            shutil.move(f"{dirname}/{basename}.cgi_results.zip", result_dir)
        else:
            if not filter_vep:
                os.makedirs(result_dir)
            shutil.move(f"{basename}.cgi_results", result_dir)
            shutil.move(f"{basename}.cgi_results.zip", result_dir)

    except FileNotFoundError:
        print("Cannot find file on disk. Please try another path.")


# querynator civic_api
@querynator_cli.command()
@click.option(
    "-v",
    "--vcf",
    help="Please provide the path to a Variant Call Format (VCF) file (Version 4.2)",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.STRING,
    help="i.e. sample name. Directory in which results will be stored.",
)
@click.option(
    "-g",
    "--genome",
    type=click.Choice(["GRCh37", "GRCh38", "NCBI36"], case_sensitive=True),
    help="Please enter the reference genome version",
    required=True,
    default="GRCh38",
)
@click.option(
    "-f",
    "--filter_vep",
    help="if set, filters out synoymous and low impact variants based on VEP annotation",
    is_flag=True,
    show_default=True,
    default=False,
)
def query_api_civic(vcf, outdir, genome, filter_vep):
    try:
        result_dir = get_unique_querynator_dir(f"{outdir}")
        dirname, basename = os.path.split(result_dir)
        if filter_vep:
            in_vcf_header, candidate_variants, removed_variants = filter_vcf_by_vep(vcf, logger)
            # create result directories
            os.makedirs(f"{result_dir}/vcf_files")
            write_vcf(in_vcf_header, removed_variants, f"{result_dir}/vcf_files/{basename}.removed_variants.vcf")
            write_vcf(in_vcf_header, candidate_variants, f"{result_dir}/vcf_files/{basename}.filtered_variants.vcf")

            logger.info("Query the Clinical Interpretations of Variants In Cancer (CIViC)")
            # run analysis
            query_civic(candidate_variants, result_dir, logger, vcf, genome, filter_vep)

        else:
            logger.info("Query the Clinical Interpretations of Variants In Cancer (CIViC)")
            query_civic(vcf, result_dir, logger, vcf, genome, filter_vep)

    except FileNotFoundError:
        print("The provided file cannot be found. Please try another path.")


# querynator create report
@querynator_cli.command()
@click.option(
    "-c",
    "--cgi_path",
    help="Path to a CGI result folder generated using the querynator",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-j",
    "--civic_path",
    help="Path to a CIViC result folder generated using the querynator",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.STRING,
    help="Name of new directory in which reports will be stored.",
)
def create_report(cgi_path, civic_path, outdir):
    # create outdir
    report_dir = get_unique_querynator_dir(outdir)
    dirname, basename = os.path.split(report_dir)
    os.makedirs(f"{report_dir}/combined_files")
    os.makedirs(f"{report_dir}/report/variant_reports")
    os.makedirs(f"{report_dir}/report/plots")

    # combine the results
    combine_civic(civic_path, report_dir, logger)
    combine_cgi(cgi_path, report_dir, logger)
    combine_cgi_civic(report_dir, logger)

    # add tiers & ranking-score to merged results
    add_tiers_and_scores_to_df(report_dir, logger)

    # create report
    create_report_htmls(report_dir, basename, civic_path, logger)


if __name__ == "__main__":
    run_querynator()
