""" Combine CIViC-VEP with CGI-VEP """

import pandas as pd
import vcf
import os



def merge_civic_cgi(alterations_vep, civic_vep):
    """
    merge CIViC and CGI alterations annotations for each variant based on the similar variant VEP annotation

    :param alterations_vep: DataFrame of variants and their VEP and CGI alterations annotations
    :type alterations_vep: pandas DataFrame
    :param civic_vep: DataFrame of variants and their VEP and CIViC alterations annotations
    :type civic_vep: pandas DataFrame
    :return: merged DataFrame of variants and their VEP & CIViC & CGI alterations annotations
    :rtype: pandas DataFrame
    """
    # drop cols created for earlier steps
    alterations_vep = alterations_vep.drop(["chr_merge_VEP", "pos_merge_VEP", "ref_merge_VEP", "alt_merge_VEP"], axis=1)

    # CIViC (merge VEP in CIViC) might be all nan for cols in which alterations does have items (CGI in VEP). 
    # so, all _VEP cols must be of the type as in alterations_vep.tsv
    for col in alterations_vep.columns:
        if col.endswith("_VEP"):
            civic_vep[col] = civic_vep[col].astype(alterations_vep[col].dtypes.name)

    # merge
    vep_civic_cgi_merge = alterations_vep.merge(civic_vep, on=list(civic_vep.columns[civic_vep.columns.str.endswith("_VEP")]), suffixes=("_cgi", "_civic"), how='left')

    return vep_civic_cgi_merge






def combine_cgi_civic(outdir, logger):
    """
    Combine cgi-vep table with the civic-vep table

    :param outdir: Path to report directory
    :type outdir: str
    :return: None
    :rtype: None
    """
    logger.info("Combining CIViC-VEP & CGI-VEP")

    # read in the files
    alterations_vep = pd.read_csv(f"{outdir}/combined_files/alterations_vep.tsv", sep="\t")
    civic_vep = pd.read_csv(f"{outdir}/combined_files/civic_vep.tsv", sep="\t")

    vep_civic_cgi_merge = merge_civic_cgi(alterations_vep, civic_vep)

    # write combined results to result dir
    vep_civic_cgi_merge.to_csv(f"{outdir}/combined_files/civic_cgi_vep.tsv", sep="\t", index=False)


    logger.info("Files combined")
