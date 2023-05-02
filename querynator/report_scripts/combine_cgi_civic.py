""" Combine CIViC-VEP with CGI-VEP """

import pandas as pd


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
    # so, all _VEP cols on which the merging is performed must be of the same type
    for i in ["chr_VEP", "pos_VEP", "ref_VEP", "alt_VEP"]:
        civic_vep[i] = civic_vep[i].astype(str)
        alterations_vep[i] = alterations_vep[i].astype(str)

    # merge
    vep_civic_cgi_merge = alterations_vep.merge(
        civic_vep, on=["chr_VEP", "pos_VEP", "ref_VEP", "alt_VEP"], suffixes=("_cgi", "_civic"), how="left"
    )

    # remove all unnecessary VEP cols
    vep_civic_cgi_merge = vep_civic_cgi_merge[
        [col for col in vep_civic_cgi_merge.columns if not col.endswith("_VEP_civic")]
    ]

    # rename VEP columns
    vep_civic_cgi_merge.columns = [
        col.replace("_VEP_cgi", "_VEP") if col.endswith("_VEP_cgi") else col for col in vep_civic_cgi_merge.columns
    ]

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
