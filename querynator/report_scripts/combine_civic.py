""" Combine the results of the CIViC query with the initial VEP annotation """

import os

import numpy as np
import pandas as pd
import vcf

from querynator.helper_functions import flatten, get_num_from_chr


def read_filtered_vcf(filtered_vcf):
    """
    Create a table containing the VEP annotation of each variant

    :param filtered_vcf: Path to the project's VEP filtered vcf
    :type filtered_vcf: str
    :return: vep table
    :rtype: pandas DataFrame
    """
    # read in filtered and normalized vcf file
    reader = vcf.Reader(open(filtered_vcf))

    vep_headers = reader.infos["CSQ"].desc.split(":")[1].strip().split("|")
    vep_headers.insert(0, "querynator_id")
    vep_headers.insert(1, "chr")
    vep_headers.insert(2, "pos")
    vep_headers.insert(3, "ref")
    vep_headers.insert(4, "alt")

    record_info = []
    for record in reader:
        vep_list = []
        # mulitple entries
        if len(record.INFO["CSQ"]) > 1:
            for vep_anno in record.INFO["CSQ"]:
                # first entry
                if len(vep_list) == 0:
                    vep_list = vep_anno.split("|")
                # all later entries
                else:
                    for i, anno in enumerate(vep_anno.split("|")):
                        if anno not in vep_list[i]:
                            vep_list[i] = ",".join([i for i in sorted(flatten([vep_list[i], anno]))])
            # remove "empty" strings (",")
            vep_list = [i if any(i.split(",")) else "" for i in vep_list]

        else:  # just one entry
            vep_list = record.INFO["CSQ"][0].split("|")

        # add coords & querynator ID
        vep_list.insert(0, int("".join(record.INFO["QID"])))
        vep_list.insert(1, get_num_from_chr(record.CHROM))
        vep_list.insert(2, record.POS)
        vep_list.insert(3, record.REF)
        vep_list.insert(4, "".join(str(i) for i in record.ALT))

        # add vep_list to final dataframe list
        record_info.append(vep_list)

    vep_df = pd.DataFrame(record_info, columns=vep_headers)

    # add VEP suffix to vep columns
    vep_df = vep_df.add_suffix("_VEP")

    vep_df.insert(0, "querynator_id", vep_df.pop("querynator_id_VEP"))

    return vep_df


def read_civic_results(civic_results):
    """
    Read in the project's CIViC annotation created by the querynator

    :param civic_results: Path to the project's CIViC annotation
    :type civic_results: str
    :return: DataFrame of CIViC resu
    :rtype: pandas DataFrame
    """
    civic_df = pd.read_csv(civic_results, sep="\t")

    # add civic suffix to civic columns
    civic_df = civic_df.add_suffix("_CIVIC")
    # add querynator_id as first column
    civic_df.insert(0, "querynator_id", civic_df.pop("querynator_id_CIVIC"))

    # add "duplicated" column --> true for duplicated rows, used later when creating report
    civic_df["duplicated_CIVIC"] = civic_df.duplicated(subset="querynator_id", keep=False)

    return civic_df


def merge_civic_vep(vep_df, civic_df):
    """
    merge vep and civic annotation for each variant based on the Querynator ID

    :param vep_df: DataFrame of variants and their VEP annotation
    :type vep_df: pandas DataFrame
    :param civic_df: DataFrame of variants and their CIViC annotation
    :type civic_df: pandas DataFrame
    :return: merged DataFrame of variants and their VEP & CIViC annotation
    :rtype: pandas DataFrame
    """
    # merge vep df into civic df (size of merge is equal to size of civic df)
    return civic_df.merge(vep_df, on="querynator_id", suffixes=("_vep", "_civic"), how="left")


def combine_civic(civic_path, outdir, logger):
    """
    Command to combine the civic results with the vcf's VEP annotation

    :param civic_path: Path to a CIViC result folder generated using the querynator
    :type civic_path: str
    :param outdir: Path to report directory
    :type outdir: str
    :return: None
    :rtype: None
    """
    logger.info("Combining CIViC & VEP")

    # get necessary files from result path
    dirname, basename = os.path.split(civic_path)
    filtered_vcf = f"{civic_path}/vcf_files/{basename}.filtered_variants.vcf"
    civic_results = f"{civic_path}/{basename}.civic_results.tsv"

    vep_df = read_filtered_vcf(filtered_vcf)
    civic_df = read_civic_results(civic_results)
    # combine results
    merged_df = merge_civic_vep(vep_df, civic_df)

    # add merged df to outdir
    merged_df.to_csv(f"{outdir}/combined_files/civic_vep.tsv", sep="\t", index=False)
