""" Combine the results of the CGI query with the initial VEP annotation """

import os

import numpy as np
import pandas as pd
import vcf

pd.options.mode.chained_assignment = None

from querynator.helper_functions import flatten, get_num_from_chr


def remove_prefix(s, prefix):
    """
    removes prefix of a string

    :param s: string which prefix should be removed
    :type s: str
    :param prefix: prefix to remove from string
    :type prefix: str
    :return: string with prefix removed
    :rtype: str
    """
    return s[len(prefix) :] if s.startswith(prefix) else s


def read_filtered_vcf(filtered_vcf):
    """
    Create a table containing the VEP annotation of each variant and positional information to connect to the alterations.tsv

    :param filtered_vcf: Path to the project's VEP filtered vcf
    :type filtered_vcf: str
    :return: vep table
    :rtype: pandas DataFrame
    """
    reader = vcf.Reader(open(filtered_vcf))
    # variant information from vcf
    vep_headers = reader.infos["CSQ"].desc.split(":")[1].strip().split("|")
    vep_headers.insert(0, "chr")
    vep_headers.insert(1, "pos")
    vep_headers.insert(2, "ref")
    vep_headers.insert(3, "alt")

    # adapted variant information to resemble cgi notation
    vep_headers.insert(4, "chr_merge")
    vep_headers.insert(5, "pos_merge")
    vep_headers.insert(6, "ref_merge")
    vep_headers.insert(7, "alt_merge")
    record_info = []
    count = 0
    for record in reader:
        vep_list = []
        counter = 0
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

        # add variant information
        vep_list.insert(0, get_num_from_chr(record.CHROM))
        vep_list.insert(1, record.POS)
        vep_list.insert(2, record.REF)
        vep_list.insert(3, "".join(str(i) for i in record.ALT))

        # add variant information for cgi merge
        chr_merge = ""
        pos_merge = ""
        ref_merge = ""
        alt_merge = ""

        # if DEL or INS, variant positional information is adapted to fit to CGI notation
        # DEL
        # if len(record.REF) > 1 & len(record.ALT[0]) == 1:
        # pos = pos+1, ref=ACGT -> CGT alt=A -> -
        if len(record.REF) > len(record.ALT[0]):
            chr_merge = get_num_from_chr(record.CHROM)
            pos_merge = record.POS + 1
            ref_merge = remove_prefix(str(record.REF), str(record.ALT[0]))
            alt_merge = "-"

        # INS
        # elif len(record.REF) == 1 & len(record.ALT[0]) > 1
        # pos = pos, ref=A -> -, alt=ACGT -> CGT
        elif len(record.REF) < len(record.ALT[0]):
            chr_merge = get_num_from_chr(record.CHROM)
            pos_merge = record.POS
            ref_merge = "-"
            alt_merge = remove_prefix(str(record.ALT[0]), str(record.REF))

        # larger "symmetric" InDels (AA -> GG) & SNPs
        else:
            chr_merge = get_num_from_chr(record.CHROM)
            pos_merge = record.POS
            ref_merge = record.REF
            alt_merge = str(record.ALT[0])

        vep_list.insert(4, chr_merge)
        vep_list.insert(5, pos_merge)
        vep_list.insert(6, ref_merge)
        vep_list.insert(7, alt_merge)

        # add vep_list to final dataframe list
        record_info.append(vep_list)

    vep_df = pd.DataFrame(record_info, columns=vep_headers)

    # add VEP suffix to vep columns
    vep_df = vep_df.add_suffix("_VEP")

    return vep_df


def extract_coords(row):
    """
    extracts coordinates from the hgvs notation provided in "alterations.tsv"

    :param row: row of a pandas DataFrame
    :type row: pandas Series
    :return: extracted coordinates
    :rtype: pandas Series
    """
    # get chromosome and strip "chr"
    chrom = row["Mutation"].split(":")[0][3:]
    init_pos = row["Mutation"].split(":")[1].split(" ")[0]
    ref, alt = row["Mutation"].split(":")[1].split(" ")[1].split(">")
    pos = ""

    # if sth other than SNPs, hgvs is the e.g. following: chr1:1234-1234 TTTCCA>-
    # DEL & INS & larger symmetric InDels
    if "-" in init_pos:
        pos = int(init_pos.split("-")[0])
    # SNPs (chr1:1234 T>A)
    else:
        pos = int(init_pos)

    return pd.Series([chrom, pos, ref, alt])


def read_modify_alterations(alterations_path):
    """
    reads in and adds positional information to alterations file

    :param alterations_path: Path to alterations file
    :type alterations_path: str
    :return: None
    :rtype: None
    """
    alterations_df = pd.read_csv(alterations_path, sep="\t")

    # extract positional information & add to df
    alterations_df[["chr", "pos", "ref", "alt"]] = alterations_df.apply(lambda x: extract_coords(x), axis=1)

    # rearrange cols
    alterations_df.insert(0, "chr", alterations_df.pop("chr"))
    alterations_df.insert(1, "pos", alterations_df.pop("pos"))
    alterations_df.insert(2, "ref", alterations_df.pop("ref"))
    alterations_df.insert(3, "alt", alterations_df.pop("alt"))

    # add suffix
    alterations_df = alterations_df.add_suffix("_CGI")

    return alterations_df


def merge_alterations_vep(vep_df, alterations_df):
    """
    merge vep and CGI alterations annotations for each variant based on positional information (chr, pos, ref, alt)

    :param vep_df: DataFrame of variants and their VEP annotation
    :type vep_df: pandas DataFrame
    :param alterations_df: DataFrame of variants and their CGI alterations annotations
    :type alterations_df: pandas DataFrame
    :return: merged DataFrame of variants and their VEP & CGI alterations annotations
    :rtype: pandas DataFrame
    """
    alterations_vep = vep_df.merge(
        alterations_df,
        left_on=["chr_merge_VEP", "pos_merge_VEP", "ref_merge_VEP", "alt_merge_VEP"],
        right_on=["chr_CGI", "pos_CGI", "ref_CGI", "alt_CGI"],
        suffixes=("_alterations", "_input"),
        how="left",
    )

    # drop duplicates, duplicates were seen to be exact duplicates here
    alterations_vep = alterations_vep.drop_duplicates(
        subset=["chr_merge_VEP", "pos_merge_VEP", "ref_merge_VEP", "alt_merge_VEP"]
    )

    return alterations_vep


def get_all_alterations(row):
    """
    extract only the alteration strings from the "Alterations" col in biomarkers.tsv

    :param row: row of a pandas DataFrame
    :type row: pandas Series
    :return: link of biomarker to all related alterations
    :rtype: list
    """
    alteration_links = [j.split(")")[0] for j in [i.split("(")[1] for i in row["Alterations"].split(", ")]]
    return alteration_links


def link_biomarkers(biomarkers_df):
    """
    add alteration-link column to "biomarkers.tsv"

    :param biomarkers_df: pd DataFrame of the projects "biomarkers.tsv"
    :type biomarkers_df: pandas DataFrame
    :return: DataFrame of biomarkers with additional alteration-link col
    :rtype: pandas DataFrame
    """
    biomarkers_df["alterations_link"] = biomarkers_df.apply(lambda x: get_all_alterations(x), axis=1)

    return biomarkers_df


def get_highest_evidence(row, biomarkers_linked):
    """
    get highest associated CGI evidence of the current alteration (A-D) from the biomarkers datafrane

    :param row: row of a pandas DataFrame
    :type row: pandas Series
    :param biomarkers_linked: pd DataFrame of the projects "biomarkers.tsv"
    :type biomarkers_linked: pandas DataFrame
    :return: highest associated evidence
    :rtype: str
    """
    if not pd.isnull(row["Protein Change_CGI"]):
        for evidence in ["A", "B", "C", "D"]:
            if not biomarkers_linked.loc[
                (biomarkers_linked["alterations_link"].str.contains(row["Protein Change_CGI"]))
                & (biomarkers_linked["Evidence"] == evidence)
            ].empty:
                return evidence
    else:
        return row["Protein Change_CGI"]  # return nan


def combine_cgi(cgi_path, outdir, logger):
    """
    Command to combine the cgi results with the vcf's VEP annotation

    :param cgi_path: Path to a CGI result folder generated using the querynator
    :type cgi_path: str
    :param outdir: Path to report directory
    :type outdir: str
    :return: None
    :rtype: None
    """
    logger.info("Combining CGI & VEP")

    # get necessary files from result path
    dirname, basename = os.path.split(cgi_path)
    filtered_vcf = f"{cgi_path}/vcf_files/{basename}.filtered_variants.vcf"
    alterations_path = f"{cgi_path}/{basename}.cgi_results/alterations.tsv"
    biomarkers_path = f"{cgi_path}/{basename}.cgi_results/biomarkers.tsv"

    # read biomarkers.tsv to pd df
    biomarkers_df = pd.read_csv(biomarkers_path, sep="\t")

    # combine cgi & vep
    vep_df = read_filtered_vcf(filtered_vcf)
    alterations_df = read_modify_alterations(alterations_path)
    merged_df = merge_alterations_vep(vep_df, alterations_df)

    # link alterations in biomarkers
    biomarkers_df = link_biomarkers(biomarkers_df)
    biomarkers_df.to_csv(f"{outdir}/combined_files/biomarkers_linked.tsv", sep="\t", index=False)

    # add CGI evidence col to merged_df

    # adapt biomarkers to only consider "complete" matches between alteration & biomarker
    biomarkers_df = biomarkers_df[biomarkers_df.BioM == "complete"]
    # biomarkers_linked["alterations_link"] = biomarkers_linked["alterations_link"].astype(str)
    biomarkers_df["alterations_link"] = biomarkers_df["alterations_link"].apply(str)
    # add CGI evidence col
    merged_df["evidence_CGI"] = merged_df.apply(lambda x: get_highest_evidence(x, biomarkers_df), axis=1)
    # write merged to report dir
    merged_df.to_csv(f"{outdir}/combined_files/alterations_vep.tsv", sep="\t", index=False)
