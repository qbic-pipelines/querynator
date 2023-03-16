""" Query the Clinical Interpretations of Variants In Cancer (CIViC) API via its python tool CIViCPY"""

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)
import gzip
import io
import os
import shutil
from datetime import date

import civicpy
import numpy as np
import pandas as pd
import pysam
from civicpy import civic

# load the civic cache (necessary to run bulk analysis)
civic.load_cache()


def check_vcf_input(vcf_path, logger):
    """
    Checks whether input is vcf-file.

    :param vcf_path: Variant Call Format (VCF) file (Version 4.2)
    :type vcf_path: str
    :return: None
    """
    header = False
    needed_cols = ["chrom", "pos", "ref", "alt"]

    if not vcf_path.endswith(".vcf"):
        logger.error("Given File does not end with '.vcf'")
        exit(1)
    else:
        for line in open(vcf_path, "r"):
            if line.startswith("#") and not line.startswith("##"):
                header = True
                missing_cols = [i for i in needed_cols if i not in line.lower()]
                if len(missing_cols) != 0:
                    logger.error(f"vcf file is missing crucial columns: {', '.join(missing_cols).upper()}")
                    exit(1)
        if not header:
            logger.error(f"vcf file requires header column!")
            exit(1)


def get_coordinates_from_vcf(vcf_path, build, logger):
    """
    Read in vcf file using "pysam",
    creates CoordinateQuery objects for each variant .
    This function does find (ref-alt):
    SNPs (A-T)
    DelIns (AA-TT)
    Deletions (TTTCA -  AT)

    :param build: Genome build version, currently only GRCh37 allowed
    :type build: str
    :return: List of CoordinateQuery objects
    :rtype: list
    """

    coord_list = []
    for record in pysam.VariantFile(vcf_path):
        for alt_base in record.alts:
            # INSERTION
            if len(record.ref) < len(alt_base):
                coord_list.append(
                    civic.CoordinateQuery(
                        chr=str(record.chrom),
                        start=int(record.start) + 1,
                        stop=int(record.start) + 2,
                        alt=alt_base[1:],
                        ref="",
                        build=build,
                    )
                )
            # DELETION
            elif len(record.ref) > len(alt_base) and len(alt_base) == 1:
                coord_list.append(
                    civic.CoordinateQuery(
                        chr=str(record.chrom),
                        start=int(record.start) + 1,
                        stop=int(record.stop),
                        alt="",
                        ref=record.ref,
                        build=build,
                    )
                )
            # SNPs, DelIns
            else:
                coord_list.append(
                    civic.CoordinateQuery(
                        chr=str(record.chrom),
                        start=int(record.start) + 1,
                        stop=int(record.stop),
                        alt=alt_base,
                        ref=record.ref,
                        build=build,
                    )
                )

    return coord_list


def access_civic_by_coordinate(coord_list):
    """
    Query CIViC API for individual variants

    :param coord_list: List of CoordinateQuery objects
    :type coord_list: list
    :return: List of CIViC variant objects of successfully queried variants
    :rtype: list
    """
    coord_list.sort()

    # bulk search to quickly focus on variants found in the civic-db
    # time-intensive search must then only be done for variants that will be hits
    bulk = civic.bulk_search_variants_by_coordinates(coord_list, search_mode="exact")

    # actual search for each variant
    variant_list = []
    for coord_obj in bulk.keys():
        variant = civic.search_variants_by_coordinates(coord_obj, search_mode="exact")
        if len(coord_obj) > 0:
            variant_list.append([coord_obj, variant])

    return variant_list


def get_variant_information_from_variant(variant_obj):
    """
    Get all variant information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Dictionary with variant information for respective CIViC variant object
    :rtype: dict
    """
    return {
        "variant_name": variant_obj.name,
        "variant_aliases": variant_obj.aliases,
        "variant_type": [i.name for i in variant_obj.types],
        "variant_clinvar_entries": variant_obj.clinvar_entries,
        "variant_entrez_id": variant_obj.entrez_id,
        "variant_entrez_name": variant_obj.entrez_name,
        "variant_hgvs_expressions": variant_obj.hgvs_expressions,
        "variant_groups": [i.name for i in variant_obj.variant_groups],
    }


def get_molecular_profile_information_from_variant(variant_obj):
    """
    Get all molecular profile information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Dictionary with molecular profile information for respective CIViC variant object
    :rtype: dict
    """
    try:
        mol_profile = variant_obj.molecular_profiles[0]
        mol_profile_dict = {
            "mol_profile_name": mol_profile.name,
            "mol_profile_definition": mol_profile.description,
            "mol_profile_score": mol_profile.molecular_profile_score,
        }
    except IndexError:
        mol_profile_dict = {"mol_profile_name": None, "mol_profile_definition": None, "mol_profile_score": None}
    return mol_profile_dict


def get_gene_information_from_variant(variant_obj):
    """
    Get all gene information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Dictionary with gene information for respective CIViC variant object
    :rtype: dict
    """
    gene = variant_obj.gene
    return {
        "gene_name": gene.name,
        "gene_aliases": gene.aliases,
        "gene_description": gene.description,
        "gene_entrez_id": gene.entrez_id,
        "gene_source": [i.name for i in gene.sources],
    }


def get_assertion_information_from_variant(variant_obj):
    """
    Get all assertion information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Dictionary with assertion information for respective CIViC variant object
    :rtype: dict
    """
    try:
        assertion = variant_obj.molecular_profiles[0].assertions[0]
        assertion_dict = {
            "assertion_name": assertion.name,
            "assertion_acmg_codes": assertion.acmg_codes,
            "assertion_amp_level": assertion.amp_level,
            "assertion_direction": assertion.assertion_direction,
            "assertion_type": assertion.assertion_type,
            "assertion_description": assertion.description,
            "assertion_disease": assertion.disease,
            "assertion_phenotypes": [i.name for i in assertion.phenotypes],
            "assertion_significance": assertion.significance,
            "assertion_status": assertion.status,
            "assertion_summary": assertion.summary,
            "assertion_therapies": [i.name for i in assertion.therapies],
            "assertion_therapy_interaction_type": assertion.therapy_interaction_type,
            "assertion_variant_origin": assertion.variant_origin,
        }
    except IndexError:
        assertion_dict = {
            "assertion_name": np.nan,
            "assertion_acmg_codes": np.nan,
            "assertion_amp_level": np.nan,
            "assertion_direction": np.nan,
            "assertion_type": np.nan,
            "assertion_description": np.nan,
            "assertion_disease": np.nan,
            "assertion_phenotypes": np.nan,
            "assertion_significance": np.nan,
            "assertion_status": np.nan,
            "assertion_summary": np.nan,
            "assertion_therapies": np.nan,
            "assertion_therapy_interaction_type": np.nan,
            "assertion_variant_origin": np.nan,
        }
    return assertion_dict


def get_evidence_information_from_variant(variant_obj):
    """
    Get all evidence from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Dictionary with evidence information for respective CIViC variant object
    :rtype: dict
    """
    try:
        evidence = variant_obj.molecular_profiles[0].assertions[0].evidence[0]
        evidence_dict = {
            "evidence_name": evidence.name,
            "evidence_description": evidence.description,
            "evidence_disease": evidence.disease,
            "evidence_level": evidence.evidence_level,
            "evidence_support": evidence.evidence_direction,
            "evidence_type": evidence.evidence_type,
            "evidence_phenotypes": [i.name for i in evidence.phenotypes],
            "evidence_rating": evidence.rating,
            "evidence_significance": evidence.significance,
            "evidence_source": evidence.source,
            "evidence_status": evidence.status,
            "evidence_therapies": [i.name for i in evidence.therapies],
            "evidence_therapy_interaction_type": evidence.therapy_interaction_type,
        }
    except IndexError:
        evidence_dict = {
            "evidence_name": np.nan,
            "evidence_description": np.nan,
            "evidence_disease": np.nan,
            "evidence_level": np.nan,
            "evidence_support": np.nan,
            "evidence_type": np.nan,
            "evidence_phenotypes": np.nan,
            "evidence_rating": np.nan,
            "evidence_significance": np.nan,
            "evidence_source": np.nan,
            "evidence_status": np.nan,
            "evidence_therapies": np.nan,
            "evidence_therapy_interaction_type": np.nan,
        }
    return evidence_dict


def get_positional_information_from_coord_obj(coord_obj):
    """
    Get information about the position of the variant in the genome

    :param coord_obj: CoordinateQuery Object to respective variant object
    :type coord_obj: CIViC CoordinateQuery Object
    :return: Dictionary with positional information for respective CIViC variant object
    :rtype: dict
    """
    return {"chr": coord_obj[0], "start": coord_obj[1], "stop": coord_obj[2], "ref": coord_obj[4], "alt": coord_obj[3]}


def concat_dicts(coord_obj, variant_obj):
    """
    Create and combine different dictionaries created for single CIViC variant object

    :param coord_obj: CoordinateQuery Object to respective variant object
    :type coord_obj: CIViC CoordinateQuery Object
    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Dictionary with all information for respective CIViC variant object
    :rtype: dict
    """
    coordinates_info = get_positional_information_from_coord_obj(coord_obj)
    variant_info = get_variant_information_from_variant(variant_obj[0])
    gene_info = get_gene_information_from_variant(variant_obj[0])
    mol_profile_info = get_molecular_profile_information_from_variant(variant_obj[0])
    assertion_info = get_assertion_information_from_variant(variant_obj[0])
    evidence_info = get_evidence_information_from_variant(variant_obj[0])

    return coordinates_info | variant_info | gene_info | mol_profile_info | assertion_info | evidence_info


def create_civic_results(variant_list, out_path, logger):
    """
    Combine result dictionaries of all CIViC variant objects
    to a table and write it to user-specified file

    :param variant_list: List of CIViC variant objects of successfully queried variants
    :type variant_list: list
    :param out_path: Name for directory in which result-table will be stored
    :type out_path: str
    """
    civic_result_df = pd.DataFrame()
    for coord_obj, variant in variant_list:
        civic_result_df = civic_result_df.append(concat_dicts(coord_obj, variant), ignore_index=True)

    logger.info("CIViC Query finished")
    logger.info("Creating Results")
    try:
        civic_result_df.to_csv(f"{out_path}/civic_results.tsv", sep="\t", index=False)
    except OSError:
        os.mkdir(out_path)
        civic_result_df.to_csv(f"{out_path}/civic_results.tsv", sep="\t", index=False)


def sort_coord_list(coord_list):
    """
    Sort the input list to the bulk search

    :param coord_list: List of CoordinateQuery objects
    :type coord_list: list
    :return: sorted coord_list
    :rtype: list
    """
    return sorted(coord_list, key=lambda x: (int(x[0]) if x[0] != "X" else np.inf, x[1], x[2]))


def add_civic_metadata(out_path):
    """
    Attach metadata to civic query

    :param out_path: Name of directory in which results are stored
    :type out_path: str
    :return: None
    """

    with open(out_path + "/metadata.txt", "w") as f:
        f.write("CIViC query date: " + str(date.today()))
        f.write("\nAPI version: " + str(civicpy.version()))
        f.close()


def query_civic(vcf_path, out_path, logger):
    """
    Command to query the CIViC API

    :param vcf_path: Variant Call Format (VCF) file (Version 4.2)
    :type vcf_path: str
    :param out_path: Name for directory in which result-table will be stored
    :type out_path: str

    """

    coord_list = get_coordinates_from_vcf(vcf_path, "GRCh37", logger)

    # list needs to be sorted for bulk search
    sort_coord_list(coord_list)

    # create result table
    create_civic_results(access_civic_by_coordinate(coord_list), out_path, logger)
    add_civic_metadata(out_path)

    logger.info("CIViC Analysis done")
