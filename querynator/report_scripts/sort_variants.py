""" Sort the variants into tiers (1-4) and provide a score for each variant to sort them within the tiers.
    The procedure is based on the Standards and Guideline provided by the AMP (Association for Molecular Pathology)
    and described by Li etal. 
    (Li MM etal. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. J Mol Diagn. 2017 Jan;19(1):4-23. doi: 10.1016/j.jmoldx.2016.10.002. PMID: 27993330; PMCID: PMC5707196.)      
"""

import numpy as np
import pandas as pd

# ============================================================================ #
#                      ASSIGN VARIANTS TO TIERS
# ============================================================================ #


def get_allele_freq_tiering(row):
    """
    checks if AF is < 0.01 to assign to tier 3 (true) or 4 (false)
    :param row: row of a pandas DataFrame
    :type row: pandas Series
    :return: True if AF < 0.01, False if AF > 0.01
    :rtype: bool
    """
    af, gnomad = row["AF_VEP"], row["gnomAD_AF_VEP"]
    bool_nan, nan_val, not_nan_val = check_nan_in_pair([af, gnomad])

    if bool_nan and nan_val == "both":
        return True
    elif bool_nan and type(nan_val) == int:
        if get_largest_af([af, gnomad][not_nan_val]) < 0.01:
            return True
        else:
            return False
    elif bool_nan == False:
        if get_largest_af(af) < 0.01 and get_largest_af(gnomad) < 0.01:
            return True
        else:
            return False


def subset_variants_into_tiers(row):
    """
    decides tier (1-4) for specific variant

    :param row: row of a pandas DataFrame
    :type row: pandas Series
    :return: variant's tier
    :rtype: str
    """
    # check whether civic already provides AMP-based Tiers
    if not pd.isnull(row["assertion_amp_level_CIVIC"]):
        if "TIER_I_" in row["assertion_amp_level_CIVIC"]:
            return "tier_1"
        elif "TIER_II_" in row["assertion_amp_level_CIVIC"]:
            return "tier_2"
        elif "TIER_III_" in row["assertion_amp_level_CIVIC"]:
            return "tier_3"
        elif "TIER_IV_" in row["assertion_amp_level_CIVIC"]:
            return "tier_4"
    else:
        if row["evidence_CGI"] in ["A", "B"] or row["evidence_level_CIVIC"] in ["A", "B"]:
            return "tier_1"
        elif row["evidence_CGI"] in ["C", "D"] or row["evidence_level_CIVIC"] in ["C", "D"]:
            return "tier_2"
        else:  # no evidence given (tier 3 or 4)
            # Tier 3 if oncogenic & MAF fits
            if row["evidence_level_CIVIC"] == "E":
                return "tier_3"
            elif (
                row["Oncogenic Summary_CGI"] not in ["non-oncogenic", "non-protein affecting"]
                and pd.isnull(row["Oncogenic Summary_CGI"]) == False
            ):
                if get_allele_freq_tiering(row):
                    return "tier_3"
                else:
                    return "tier_4"
            else:
                return "tier_4"


# ============================================================================ #
#                       ADD SCORE TO RANK WITHIN TIERS
# ============================================================================ #


def check_nan_in_pair(pair):
    """
    check if one of the values is nan

    :param pair: list of values that are either string or nan
    :return: True, index of the nan value and index of the non-nan value
    :rtype: list
    """
    if pd.isnull(pair[0]) and pd.isnull(pair[1]):
        return [True, "both", "both"]
    elif pd.isnull(pair[0]):
        return [True, 0, 1]
    elif pd.isnull(pair[1]):
        return [True, 1, 0]
    else:
        return [False, "neither", "neither"]


def get_largest_af(af_string):
    """
    extracts the largest allele frequency from a string of allele frequencies.
    If single str is given, gives it out as int

    :param af_string: comma-separated string of allele frequencies (af)
    :type af_string: str
    :return: largest af in string
    :rtype: int
    """
    if type(af_string) in [int, float]:
        return af_string
    else:
        freq_ints = [int(i) for i in af_string.split(",")]
        return max(freq_ints)


def extract_num(s):
    """
    extracts the number from a string with pattern (e.g.) "benign(0.001)"

    :param s: string to extract number from
    :type s: str
    :return: extracted number
    :rtype: float
    """
    return float(s.split("(")[1].split(")")[0]) if s != "" else np.nan


def get_largest_path_score(ps_string, ps_score):
    """
    extracts the largest pathogenicity score from a string of pathogenicity scores.
    If single str is given, gives it out as int

    :param af_string: comma-separated string of pathogenicity scores
    :type af_string: str
    :param ps_score: respective pathogenicity score (SIFT or PolyPhen2)
    :type ps_score: str
    :return: largest af in string
    :rtype: int
    """
    if ps_score == "sift":
        return np.nanmax([extract_num(i) for i in ps_string.split(",")])
    if ps_score == "polyphen":
        return np.nanmin([extract_num(i) for i in ps_string.split(",")])


def generate_allele_freq_score(af, gnomad):
    """
    generates the variantMTB score for the allele frequency of a specific variant

    :param af: a variant's associated allele frequencies
    :type af: str
    :param gnomad: a variant's associated gnomAD frequencies
    :type gnomad: str
    :return: allele frequency score
    :rtype: int
    """
    # check for nans and get which one is nan
    bool_nan, nan_val, not_nan_val = check_nan_in_pair([af, gnomad])

    if bool_nan and nan_val == "both":
        return 2
    # one is nan
    elif bool_nan and type(nan_val) == int:
        # very low af as proposed by VIC
        if get_largest_af([af, gnomad][not_nan_val]) < 0.005:
            return 2
        elif get_largest_af([af, gnomad][not_nan_val]) < 0.01:
            return 1
        elif get_largest_af([af, gnomad][not_nan_val]) >= 0.01:
            return -10
    # both are given
    elif bool_nan == False:
        if get_largest_af(af) < 0.005 and get_largest_af(gnomad) < 0.005:
            return 2
        elif get_largest_af(af) < 0.01 and get_largest_af(gnomad) < 0.01:
            return 1
        elif get_largest_af(af) >= 0.01 or get_largest_af(gnomad) >= 0.01:
            return -10


def get_consequence_score(consequence_str):
    """
    Returns the highest score of the associated CIViC variant types for a variant

    :param consequence_str: a variant's consequence
    :type consequence_str: str
    :return: consequence score
    :rtype: int
    """
    # consequence scoring lists
    high_impact = [
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
    ]
    moderate_impact = ["inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant"]
    low_impact = [
        "splice_region_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
    ]
    modifier_impact = [
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        "intergenic_variant",
    ]

    score_dict_civic = {
        "Frameshift Truncation": 2,
        "Gene Variant": 0,
        "Transcript Variant": 0,
        "Loss Of Function Variant": 2,
        "Gain Of Function Variant": 2,
        "Exon Variant": 0,
        "Transcript Fusion": 2,
        "Gene Fusion": 2,
        "Transcript Translocation": 2,
        "Feature Translocation": 2,
        "Transcript Fusion": 2,
        "Wild Type": -2,
        "Loss Of Heterozygosity": 2,
        "Copy Number Change": 2,
        "Exon Loss Variant": 2,
    }

    if consequence_str in high_impact:
        return 2
    elif consequence_str in moderate_impact:
        return 1
    elif consequence_str in low_impact:
        return -2
    elif consequence_str in modifier_impact:
        return 0
    else:
        if consequence_str in score_dict_civic:
            return score_dict_civic[consequence_str]
        # unknown consequence provided (usually CIViC variant type that is not given here: https://civic.readthedocs.io/en/latest/model/variants/types.html)
        else:
            print(consequence_str)
            print("whoops, I dont know this consequence")
            return 0


def get_civic_consequence_score(civic_consequence):
    """
    Translates the CIViC consequence (variant type) (https://civic.readthedocs.io/en/latest/model/variants/types.html) nomenclature
    in the CGI/VEP nomenclature when possible and scores the variant

    :param civic_consequence: a variant's consequence as given by CIViC
    :type civic_consequence: str
    :return: CIViC consequence score
    :rtype: int
    """
    translation_dict = {
        "Missense Variant": "missense_variant",
        "Stop Gained": "stop_gained",
        "Protein Altering Variant": "protein_altering_variant",
        "Inframe Deletion": "inframe_deletion",
        "Inframe Insertion": "inframe_insertion",
        "Splice Acceptor Variant": "splice_acceptor_variant",
        "Splice Donor Variant": "splice_donor_variant",
        "Transcript Amplification": "transcript_amplification",
        "Transcript Ablation": "transcript_ablation",
        "5 Prime UTR Variant": "5_prime_UTR_variant",
        "3_prime_UTR_variant": "5_prime_UTR_variant",
        "Synonymous Variant": "synonymous_variant",
        "Frameshift Truncation": "Frameshift Truncation",
        "Gene Variant": "Gene Variant",
        "Transcript Variant": "Transcript Variant",
        "Loss Of Function Variant": "Loss Of Function Variant",
        "Gain Of Function Variant": "Gain Of Function Variant",
        "Exon Variant": "Exon Variant",
        "Transcript Fusion": "Transcript Fusion",
        "Gene Fusion": "Gene Fusion",
        "Transcript Translocation": "Transcript Translocation",
        "Feature Translocation": "Feature Translocation",
        "Transcript Fusion": "Transcript Fusion",
        "Wild Type": "Wild Type",
        "Loss Of Heterozygosity": "Loss Of Heterozygosity",
        "Copy Number Change": "Copy Number Change",
        "Exon Loss Variant": "Exon Loss Variant",
    }

    if not pd.isnull(civic_consequence):
        return max(
            [
                get_consequence_score(translation_dict[i.strip()]) if i in translation_dict else 0
                for i in civic_consequence.split(",")
            ]
        )
    else:
        return 0  # unknown consequence provided (usually CIViC variant type that is not given here: https://civic.readthedocs.io/en/latest/model/variants/types.html)


def get_cgi_consequence_score(cgi_consequence):
    """
    Scores the CGI consequence

    :return: CGI consequence score
    :rtype: int
    :param civic_consequence: a variant's consequence as given by CIViC
    :type civic_consequence: str
    """
    if not pd.isnull(cgi_consequence):
        return get_consequence_score(cgi_consequence)
    else:
        return 0


def generate_pathogenicity_score_score(sift, polyphen):
    """
    generates the variantMTB score for the pathogenicity scores (SIFT & PolyPhen2) of a specific variant

    :param af: a variant's associated allele frequencies
    :type af: str
    :param gnomad: a variant's associated gnomAD frequencies
    :type gnomad: str
    :return: allele frequency score
    :rtype: int
    """
    # check for nans and get which one is nan
    bool_nan, nan_val, not_nan_val = check_nan_in_pair([sift, polyphen])

    if bool_nan and nan_val == "both":
        return 0

    # one is nan
    elif bool_nan and type(nan_val) == int:
        # SIFT is nan
        if not_nan_val == 0:
            if get_largest_path_score([sift, polyphen][not_nan_val], "polyphen") <= 0.05:
                return 1
            else:
                return 0
        # PolyPhen is nan
        elif not_nan_val == 1:
            if get_largest_path_score([sift, polyphen][not_nan_val], "sift") >= 0.85:
                return 1
            else:
                return 0
    # both are given
    elif bool_nan == False:
        if get_largest_path_score(sift, "sift") <= 0.05 and get_largest_path_score(polyphen, "polyphen") >= 0.85:
            return 1
        else:
            return 0


def generate_consequence_score(cgi_consequence, civic_consequence):
    """
    generates the variantMTB score for the consequence of a specific variant

    :param cgi_consequence: a variant's associated consequence provided by CGI
    :type cgi_consequence: str
    :param civic_consequence: a variant's consequence as given by CIViC
    :type civic_consequence: str
    :return: consequence score
    :rtype: int
    """
    return max([get_cgi_consequence_score(cgi_consequence), get_civic_consequence_score(civic_consequence)])


def generate_evidence_score(evidence_col):
    """
    generates the variantMTB score for the evidence level of one of the Knowledgebases of a specific variant

    :param evidence_col: evidence value of one of the Knowledgebases
    :return: evidence score
    :rtype: int
    """
    if evidence_col == "A":
        return 5
    elif evidence_col == "B":
        return 3
    elif evidence_col in ["C", "D"]:
        return 2
    elif evidence_col == "E":
        return 1
    else:
        return 0


def scoring_variants(row):
    """
    adds score to each variant to rank them within their tiers.

    :param row: row of a pandas DataFrame
    :type row: pandas Series
    :return: variant's score
    :rtype: int
    """
    score = 0

    # present in both KBs
    if (
        row["Oncogenic Summary_CGI"] not in ["non-oncogenic", "non-protein affecting"]
        and not pd.isnull(row["Oncogenic Summary_CGI"])
        and not pd.isnull(row["chr_CIVIC"])
    ):
        score += 3

    # multiple external  oncogenic annotations
    if not pd.isnull(row["External oncogenic annotation_CGI"]):
        score += len(row["External oncogenic annotation_CGI"].split(","))

    # drugs associated
    if (
        not pd.isnull(row["evidence_therapies_CIVIC"])
        or not pd.isnull(row["assertion_therapies_name_CIVIC"])
        or not pd.isnull(row["evidence_CGI"])
    ):
        score += 2

    # Evidence associated CIViC
    score += generate_evidence_score(row["evidence_level_CIVIC"])

    # Evidence associated CGI
    score += generate_evidence_score(row["evidence_CGI"])

    # consequence
    score += generate_consequence_score(row["Consequence_CGI"], row["variant_type_CIVIC"])

    # Allele Freq
    score += generate_allele_freq_score(row["AF_VEP"], row["gnomAD_AF_VEP"])

    # SIFT / PolyPhen2
    score += generate_pathogenicity_score_score(row["SIFT_VEP"], row["PolyPhen_VEP"])

    return score


# ============================================================================ #
#                       REPORT TIERS & SCORES
# ============================================================================ #


def add_tiers_and_scores_to_df(outdir, logger):
    """
    Assigning variants from combined CGI-CIViC-VEP table to tiers
    and give them a score to rank them in these tiers

    :param outdir: Path to report directory
    :type outdir: str
    :return: None
    :rtype: None
    """
    logger.info("Assigning variants to tiers")
    vep_civic_cgi_merge = pd.read_csv(f"{outdir}/combined_files/civic_cgi_vep.tsv", sep="\t")

    # add tiers
    vep_civic_cgi_merge["report_tier"] = vep_civic_cgi_merge.apply(lambda x: subset_variants_into_tiers(x), axis=1)

    # add ranking-score
    vep_civic_cgi_merge["ranking_score"] = vep_civic_cgi_merge.apply(lambda x: scoring_variants(x), axis=1)

    # write tiers & scores to result dir
    vep_civic_cgi_merge.to_csv(f"{outdir}/combined_files/civic_cgi_vep.tsv", sep="\t", index=False)

    logger.info("Tiers assigned")
