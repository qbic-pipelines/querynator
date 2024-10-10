""" Query the Clinical Interpretations of Variants In Cancer (CIViC) API via its python tool CIViCPY """

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

import os
import random
from datetime import date
from os.path import abspath, dirname

import civicpy
import numpy as np
import pandas as pd
import vcf
from civicpy import civic

from querynator.helper_functions import (
    get_num_from_chr,
    gunzip_compressed_files,
    gzipped,
    ontology,
)


def check_vcf_input(vcf_path, logger):
    """
    Checks whether input is vcf-file with all necessary columns.

    :param vcf_path: Variant Call Format (VCF) file (Version >= 4.0)
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


def vcf_file(vcf_path):
    """
    Checks whether input is vcf-file.

    :param vcf_path: Variant Call Format (VCF) file (Version 4.2)
    :type vcf_path: str
    :return: None
    """
    if vcf_path.endswith(".vcf") or vcf_path.endswith(".vcf.gz"):
        return True
    else:
        return False


def get_coordinates_from_vcf(input, build, logger):
    """
    Read in vcf file using "pyVCF3",
    creates CoordinateQuery objects for each variant.
    This function does find (ref-alt):
    SNPs (A-T)
    DelIns (AA-TT)
    Deletions (TTTCA -  AT)

    :param input: list of pyVCF3 records or vcf file to query
    :type input: list or str
    :param build: reference genome
    :type build: str
    :return: CoordinateQuery objects
    :rtype: list
    """
    if type(input) == list:
        variant_file = input
    else:
        if vcf_file(input):
            if gzipped(input):
                variant_file = vcf.Reader(open(gunzip_compressed_files(input, logger)))
            else:
                variant_file = vcf.Reader(open(input))

    coord_dict = {}
    for record in variant_file:
        if "QID" in record.INFO.keys():
            querynator_id = record.INFO["QID"]
        else:
            # filter_vep not applied and no rerun with filtered vcf, generate random QID for following steps which will not be reported in results
            querynator_id = random.randint(1000000, 9999999)
        for alt_base in record.ALT:
            # INSERTION
            if len(record.REF) < len(alt_base):
                coord_dict.update(
                    {
                        civic.CoordinateQuery(
                            chr=get_num_from_chr(record.CHROM),
                            start=int(record.start) + 1,
                            stop=int(record.start) + 2,
                            alt=str(alt_base)[1:],
                            ref="",
                            build=build,
                        ): querynator_id
                    }
                )
            # DELETION
            elif len(record.REF) > len(alt_base) and len(alt_base) == 1:
                coord_dict.update(
                    {
                        civic.CoordinateQuery(
                            chr=get_num_from_chr(record.CHROM),
                            start=int(record.start) + 1,
                            stop=int(record.end),
                            alt="",
                            ref=record.REF,
                            build=build,
                        ): querynator_id
                    }
                )
            # SNPs, DelIns
            else:
                coord_dict.update(
                    {
                        civic.CoordinateQuery(
                            chr=get_num_from_chr(record.CHROM),
                            start=int(record.start) + 1,
                            stop=int(record.end),
                            alt=str(alt_base),
                            ref=record.REF,
                            build=build,
                        ): querynator_id
                    }
                )

    return coord_dict


def append_to_dict(dict1, dict2):
    """
    appends values of a dictionary to another dictionary with lists as values
    :param dict1: dictionary to append to
    :type dict1: dict
    :param dict2: dictionary with values to append
    :type: dict2: dict
    :return: appended dict
    :rtype: dict
    """
    for key, value in dict2.items():
        if key in dict1:
            if type(dict1[key]) == list:
                dict1[key].append(value)
            else:
                if dict1[key] == "":
                    dict1[key] = [dict2[key]]
                else:
                    dict1[key] = [dict1[key], dict2[key]]

    return dict1


def smoothen_dict(dict, s):
    """
    makes string out of lists

    :param dict: dict with lists as values
    :type dict: dict
    :param s: True if string and special string character needed
    :type s: bool
    :return: dict with strings as values
    :rtype: dict
    """
    for key in dict:
        if type(dict[key]) == list:
            filtered_list = [i if i is not None else "" for i in dict[key]]
            if s:
                if key in ["evidence_source", "evidence_description"]:
                    dict[key] = "|".join(str(i) for i in filtered_list)
                else:
                    dict[key] = ",".join(str(i) for i in filtered_list)
            else:
                dict[key] = ",".join(str(i) for i in filtered_list)
    return dict


def access_civic_by_coordinate(coord_dict, logger, build):
    """
    Query CIViC API for individual variants

    :param coord_list: List of CoordinateQuery objects
    :type coord_list: list
    :param build: reference genome
    :type build: str
    :return: CIViC variant objects of successfully queried variants
    :rtype: list
    """

    # bulk search to quickly focus on variants found in the civic-db
    # time-intensive search must then only be done for variants that will be hits
    # only possible for GRCh37 ref genome
    if build == "GRCh37":
        bulk = civic.bulk_search_variants_by_coordinates(list(coord_dict.keys()), search_mode="exact")
        # reconnect coordinates that passed bulk search and respective IDs
        input_dict = {key: coord_dict[key] for key in bulk}
    else:
        input_dict = coord_dict

    # actual search for each variant
    variant_list = []
    for coord_obj, querynator_id in input_dict.items():
        variant = civic.search_variants_by_coordinates(coord_obj, search_mode="exact")

        if variant != None and len(variant) > 0:
            for variant_obj in variant:
                variant_list.append([{coord_obj: querynator_id}, [variant_obj]])

    # break if no variants are found
    if variant_list == None:
        logger.error("No hits are found in CIViC for this vcf file")
        exit(1)

    return variant_list


def get_variant_information_from_variant(variant_obj):
    """
    Get all variant information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Variant information for respective CIViC variant object
    :rtype: dict
    """
    return {
        "variant_name": variant_obj.name,
        "variant_aliases": ", ".join(variant_obj.aliases),
        "variant_type": ", ".join([i.name for i in variant_obj.types]),
        "variant_clinvar_entries": ", ".join(variant_obj.clinvar_entries),
        "variant_entrez_id": variant_obj.entrez_id,
        "variant_entrez_name": variant_obj.entrez_name,
        "variant_hgvs_expressions": ", ".join(variant_obj.hgvs_expressions),
        "variant_groups": ", ".join([i.name for i in variant_obj.variant_groups]),
    }


def get_molecular_profile_information_from_variant(variant_obj):
    """
    Get all molecular profile information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Molecular profile information for respective CIViC variant object
    :rtype: dict
    """
    mol_profile_dict = {"mol_profile_name": "", "mol_profile_definition": "", "mol_profile_score": ""}

    for mol_profile in variant_obj.molecular_profiles:
        try:
            # mol_profile = variant_obj.molecular_profiles[0]
            new_dict = {
                "mol_profile_name": mol_profile.name,
                "mol_profile_definition": mol_profile.description,
                "mol_profile_score": mol_profile.molecular_profile_score,
            }
        except IndexError:
            new_dict = {"mol_profile_name": None, "mol_profile_definition": None, "mol_profile_score": None}
        mol_profile_dict = append_to_dict(mol_profile_dict, new_dict)

    return smoothen_dict(mol_profile_dict, False)


def get_gene_information_from_variant(variant_obj):
    """
    Get all gene information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Gene information for respective CIViC variant object
    :rtype: dict
    """
    gene = variant_obj.gene
    return {
        "gene_name": gene.name,
        "gene_aliases": ", ".join(gene.aliases),
        "gene_description": gene.description,
        "gene_entrez_id": gene.entrez_id,
        "gene_source": ", ".join([i.name for i in gene.sources]),
    }


def get_assertion_information_from_variant(variant_obj):
    """
    Get all assertion information from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :return: Assertion information for respective CIViC variant object
    :rtype: dict
    """
    assertion_dict = {
        "assertion_name": "",
        "assertion_acmg_codes": "",
        "assertion_acmg_codes_description": "",
        "assertion_amp_level": "",
        "assertion_direction": "",
        "assertion_type": "",
        "assertion_description": "",
        "assertion_disease_name": "",
        "assertion_disease_doid": "",
        "assertion_disease_url": "",
        "assertion_disease_aliases": "",
        "assertion_phenotypes": "",
        "assertion_significance": "",
        "assertion_status": "",
        "assertion_summary": "",
        "assertion_therapies_name": "",
        "assertion_therapies_ncit_id": "",
        "assertion_therapies_aliases": "",
        "assertion_therapies_interaction_type": "",
        "assertion_variant_origin": "",
    }
    for mol_prof in variant_obj.molecular_profiles:
        for assertion in mol_prof.assertions:
            try:
                new_dict = {
                    "assertion_name": assertion.name,
                    "assertion_acmg_codes": ", ".join([i.code for i in assertion.acmg_codes]),
                    "assertion_acmg_codes_description": ", ".join([i.description for i in assertion.acmg_codes]),
                    "assertion_amp_level": assertion.amp_level,
                    "assertion_direction": assertion.assertion_direction,
                    "assertion_type": assertion.assertion_type,
                    "assertion_description": assertion.description,
                    "assertion_disease_name": ", ".join([i.name for i in assertion.disease]),
                    "assertion_disease_doid": ", ".join([i.doid for i in assertion.disease]),
                    "assertion_disease_url": ", ".join([i.disease_url for i in assertion.disease]),
                    "assertion_disease_aliases": ", ".join([i.aliases for i in assertion.disease]),
                    "assertion_phenotypes": ", ".join([i.name for i in assertion.phenotypes]),
                    "assertion_significance": assertion.significance,
                    "assertion_status": assertion.status,
                    "assertion_summary": assertion.summary,
                    "assertion_therapies_name": ", ".join([i.name for i in assertion.therapies]),
                    "assertion_therapies_ncit_id": ", ".join([i.ncit_id for i in assertion.therapies]),
                    "assertion_therapies_aliases": ", ".join([", ".join(i.aliases) for i in assertion.therapies]),
                    "assertion_therapies_interaction_type": assertion.therapy_interaction_type,
                    "assertion_variant_origin": assertion.variant_origin,
                }
            except IndexError:
                new_dict = {
                    "assertion_name": np.nan,
                    "assertion_acmg_codes": np.nan,
                    "assertion_acmg_codes_description": np.nan,
                    "assertion_amp_level": np.nan,
                    "assertion_direction": np.nan,
                    "assertion_type": np.nan,
                    "assertion_description": np.nan,
                    "assertion_disease_name": np.nan,
                    "assertion_disease_doid": np.nan,
                    "assertion_disease_url": np.nan,
                    "assertion_disease_aliases": np.nan,
                    "assertion_phenotypes": np.nan,
                    "assertion_significance": np.nan,
                    "assertion_status": np.nan,
                    "assertion_summary": np.nan,
                    "assertion_therapies_name": np.nan,
                    "assertion_therapies_ncit_id": np.nan,
                    "assertion_therapies_aliases": np.nan,
                    "assertion_therapies_interaction_type": np.nan,
                    "assertion_variant_origin": np.nan,
                }
            assertion_dict = append_to_dict(assertion_dict, new_dict)
    return smoothen_dict(assertion_dict, False)


def get_evidence_information_from_variant(variant_obj, diseases, evidence_filters):
    """
    Get all evidence from a single CIViC variant object

    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :param diseases: the patients cancer type as Disease Ontology Name or id and a list of allowed diseases that will count as matches
    :type diseases: tuple (Ontology.Term, list)
    :param evidence_filters: evidence filters
    :type evidence_filters: dict
    :return: Evidence information for respective CIViC variant object
    :rtype: dict
    """

    evidence_dict = {
        "evidence_name": "",
        "evidence_description": "",
        "evidence_disease": "",
        "evidence_level": "",
        "evidence_support": "",
        "evidence_type": "",
        "evidence_phenotypes": "",
        "evidence_rating": "",
        "evidence_significance": "",
        "evidence_source": "",
        "evidence_status": "",
        "evidence_therapies": "",
        "evidence_therapy_interaction_type": "",
    }

    for mol_prof in variant_obj.molecular_profiles:
        for evidence in mol_prof.evidence:
            if filter_evidence(evidence, evidence_filters):
                try:
                    # evidence = variant_obj.molecular_profiles[0].evidence[0]
                    new_dict = {
                        "evidence_name": evidence.name,
                        "evidence_description": evidence.description,
                        "evidence_disease": evidence.disease.name if type(evidence.disease) != dict else np.nan,
                        "evidence_level": evidence.evidence_level,
                        "evidence_support": evidence.evidence_direction,
                        "evidence_type": evidence.evidence_type,
                        "evidence_phenotypes": "+".join([i.name for i in evidence.phenotypes]),
                        "evidence_rating": evidence.rating,
                        "evidence_significance": evidence.significance,
                        "evidence_source": evidence.source.name,
                        "evidence_status": evidence.status,
                        "evidence_therapies": "+".join([i.name for i in evidence.therapies]),
                        "evidence_therapy_interaction_type": evidence.therapy_interaction_type,
                    }
                except IndexError:
                    new_dict = {
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

                # check special rules for evidences and cancer types
                disease, allowed_diseases = diseases
                if disease and allowed_diseases and evidence.disease and evidence.disease.name not in allowed_diseases:
                    if evidence.evidence_level == "A":
                        # AMP/ASCO/CAT: level A evidence for other tumor is level C for this tumor
                        new_dict["evidence_level"] = "C"
                        new_dict["evidence_disease"] = disease.name
                    else:
                        # evidence item is irrelevant for this tumor
                        continue

                evidence_dict = append_to_dict(evidence_dict, new_dict)

    return smoothen_dict(evidence_dict, True)


def filter_evidence(evidence, filters) -> bool:
    """Check if the evidence item is relevant and passes the filters

    :param evidence: CIViC evidence object
    :type evidence: civicpy.evidence.Evidence
    :param filters: dict of property:[accepted-values] pairs to filter the evidence items
    :type filters: dict
    :return: True if the evidence item passes the filters, False otherwise
    :rtype: bool
    """
    for prop, values in filters.items():
        if prop == "level".casefold():
            if evidence.evidence_level.casefold() > min(values):
                # lower level has higher char value A < B < C
                return False
        elif prop == "type".casefold():
            if evidence.evidence_type.casefold() not in values:
                return False
        elif prop == "significance".casefold():
            if evidence.significance.casefold() not in values:
                return False
        elif prop == "rating".casefold():
            if evidence.rating.casefold() not in values:
                return False
        elif prop == "status".casefold():
            if evidence.status.casefold() not in values:
                return False
        elif prop == "direction".casefold():
            if evidence.evidence_direction.casefold() not in values:
                return False
    return True


def get_positional_information_from_coord_obj(coord_obj):
    """
    Get information about the position of the variant in the genome

    :param coord_obj: CoordinateQuery Object to respective variant object
    :type coord_obj: CIViC CoordinateQuery Object
    :return: Positional information for respective CIViC variant object
    :rtype: dict
    """
    return {"chr": coord_obj[0], "start": coord_obj[1], "stop": coord_obj[2], "ref": coord_obj[4], "alt": coord_obj[3]}


def get_querynator_id(querynator_id):
    """
    Get the querynator id in dict format

    :param querynator_id: Querynator id
    :type querynator_id: str
    :return: Querynator id for respective CIViC variant object
    :rtype: dict
    """
    return {"querynator_id": querynator_id}


def concat_dicts(coord_id_dict, variant_obj, diseases, filter_vep, evidence_filters):
    """
    Create and combine different dictionaries created for single CIViC variant object

    :param coord_obj: CoordinateQuery Object to respective variant object
    :type coord_obj: CIViC CoordinateQuery Object
    :param variant_obj: single CIViC variant object
    :type variant_ob: CIViC variant object
    :param diseases: the patients cancer type as Disease Ontology Name or id and a list of allowed diseases that will count as matches
    :type diseases: tuple (Ontology.Term, list)
    :param filter_vep: flag whether VEP based filtering should be performed
    :type filter_vep: bool
    :param evidence_filters: evidence filters
    :type evidence_filters: dict
    :return: All information for respective CIViC variant object
    :rtype: dict
    """
    coordinates_info = get_positional_information_from_coord_obj(list(coord_id_dict.keys())[0])
    variant_info = get_variant_information_from_variant(variant_obj[0])
    gene_info = get_gene_information_from_variant(variant_obj[0])
    mol_profile_info = get_molecular_profile_information_from_variant(variant_obj[0])
    assertion_info = get_assertion_information_from_variant(variant_obj[0])
    evidence_info = get_evidence_information_from_variant(variant_obj[0], diseases, evidence_filters)
    if filter_vep:
        querynator_id_info = get_querynator_id(coord_id_dict[list(coord_id_dict.keys())[0]])
        return {
            **coordinates_info,
            **querynator_id_info,
            **variant_info,
            **gene_info,
            **mol_profile_info,
            **assertion_info,
            **evidence_info,
        }
    else:
        return {**coordinates_info, **variant_info, **gene_info, **mol_profile_info, **assertion_info, **evidence_info}


def create_civic_results(variant_list, out_path, disease, logger, filter_vep, evidence_filters):
    """
    Combine result dictionaries of all CIViC variant objects
    to a table and write it to user-specified file

    :param variant_list: List of CIViC variant objects of successfully queried variants
    :type  variant_list: list
    :param out_path: Name for directory in which result-table will be stored
    :type  out_path: str
    :param disease: the patients cancer type as Disease Ontology Name
    :type  disease: str
    :param filter_vep: flag whether VEP based filtering should be performed
    :type  filter_vep: bool
    :param evidence_filters: evidence filters
    :type  evidence_filters: dict
    :return: None
    :rtype: None
    """
    civic_result_df = pd.DataFrame()

    if disease:
        proj_root = dirname(dirname(abspath(__file__)))
        doid_file = os.path.join(proj_root, "helper_functions/doid.obo")
        doid = ontology.DiseaseOntology(doid_file)
        diseases = doid.get(disease), doid.get_all_ancestors(disease)
        logger.info(f"Mapped specified disease {disease} to Disease Ontology (DO) {str(diseases[0])}")
    else:
        diseases = None, None

    for coord_id_dict, variant in variant_list:
        civic_result_df = civic_result_df.append(
            concat_dicts(coord_id_dict, variant, diseases, filter_vep, evidence_filters), ignore_index=True
        )

    logger.info("CIViC Query finished")
    logger.info("Creating Results")
    try:
        civic_result_df.to_csv(f"{out_path}/{os.path.basename(out_path)}.civic_results.tsv", sep="\t", index=False)
    except OSError:
        os.makedirs(out_path)
        civic_result_df.to_csv(f"{out_path}/{os.path.basename(out_path)}.civic_results.tsv", sep="\t", index=False)


def sort_coord_list(coord_dict):
    """
    Sort the input list to the bulk search

    :param coord_list: List of CoordinateQuery objects
    :type coord_list: list
    :return: sorted coordinates
    :rtype: list
    """
    return {
        key: value
        for key, value in sorted(
            coord_dict.items(),
            key=lambda x: (
                int(x[0][0]) if x[0][0] != "X" and x[0][0] != "Y" and x[0][0] != "M" else sort_rules(x[0][0]),
                x[0][1],
                x[0][2],
            ),
        )
    }


def sort_rules(s):
    """
    Set rules to correctly sort chromosomes X,Y,M

    :param s: "string" chromosome (X,Y,M)
    :type s: str
    :return: integer to sort by
    :rtype: int
    """
    if s == "X":
        return 100
    elif s == "Y":
        return 1000
    elif s == "M":
        return 10000


def add_civic_metadata(out_path, input_file, search_mode, genome, filter_vep):
    """
    Attach metadata to civic query

    :param out_path: Name of directory in which results are stored
    :type out_path: str
    :param input_file: path of original input file
    :type input_file: str
    :param search_mode: search mode used in CIViC Query
    :type search_mode: str
    :param filter_vep: flag whether VEP based filtering should be performed
    :type filter_vep: bool
    :return: None
    :rtype: None
    """

    with open(out_path + "/metadata.txt", "w") as f:
        f.write("CIViC query date: " + str(date.today()))
        f.write("\nCIViCpy version: " + str(civicpy.version()))
        f.write("\nSearch mode: " + str(search_mode))
        f.write("\nReference genome: " + str(genome))
        if filter_vep:
            f.write("\nFiltered out synonymous & low impact variants based on VEP annotation")
        f.write("\nInput File: " + str(input_file))
        f.close()


def query_civic(vcf, out_path, logger, input_file, genome, disease, filter_vep, evidence_filters):
    """
    Command to query the CIViC API

    :param vcf: Variant Call Format (VCF) file (Version 4.2) or list of pyVCF3 variant records
    :type vcf: str or list
    :param out_path: Name for directory in which result-table will be stored
    :type out_path: str
    :param input_file: path of original input file
    :type input_file: str
    :param disease: the patients cancer type as Disease Ontology Name
    :type disease: str
    :param filter_vep: flag whether VEP based filtering should be performed
    :type filter_vep: bool
    :param evidence_filters: evidence filters
    :type evidence_filters: dict
    :return: None
    :rtype: None
    """
    # necessary for bulk run
    logger.info("Updating CIViCpy Cache")
    civic.load_cache()

    logger.info("Querying")

    coord_dict = get_coordinates_from_vcf(vcf, genome, logger)

    # coordinates needs to be sorted for bulk search
    coord_dict = sort_coord_list(coord_dict)

    # create result table
    create_civic_results(
        access_civic_by_coordinate(coord_dict, logger, genome), out_path, disease, logger, filter_vep, evidence_filters
    )
    add_civic_metadata(out_path, input_file, "exact", genome, filter_vep)

    logger.info("CIViC Analysis done")
