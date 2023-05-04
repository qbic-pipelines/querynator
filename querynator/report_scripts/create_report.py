""" Create one report of the querynator results and individual reports for each variant  """

import base64
import os
from io import BytesIO

import pandas as pd
from matplotlib import pyplot as plt
from pretty_html_table import build_table
from upsetplot import from_indicators, plot

from querynator.helper_functions import flatten

# ============================================================================ #
#                               Overall Report
# ============================================================================ #


def get_KB_count(df):
    """
    Get the count of variants for each Knowledgebase to add to the pieplot

    :param df: Variant dataframe
    :type df: pandas DataFrame
    :return: Count of variants for each Knowledgebase
    :rtype: pandas DataFrame
    """
    cgi_count = (df["sources"] == "cgi").sum() + (df["sources"] == "cgi,civic").sum()
    civic_count = sum(df["sources"] == "civic") + sum(df["sources"] == "cgi,civic")
    both_count = (df["sources"] == "cgi,civic").sum()
    none_count = (df["sources"] == "").sum()

    return [cgi_count, civic_count, both_count, none_count]


def get_sources(row):
    """
    Get the Knowledgebases containing information about the variant

    :param row: Row of the variant dataframe
    :type row: pandas DataFrame row
    """
    if not pd.isnull(row["Oncogenic Summary_CGI"]) and not pd.isnull(row["chr_CIVIC"]):
        return "cgi,civic"
    elif not pd.isnull(row["Oncogenic Summary_CGI"]) and pd.isnull(row["chr_CIVIC"]):
        return "cgi"
    elif pd.isnull(row["Oncogenic Summary_CGI"]) and not pd.isnull(row["chr_CIVIC"]):
        return "civic"
    else:
        return ""


def assign_comb_evidence_labels(row):
    """
    Assign the evidence labels for each Knowledgebase

    :param row: Row of the variant dataframe
    :type row: pandas DataFrame row
    :return: Evidence labels for each Knowledgebase
    :rtype: str
    """
    cgi_str = "-"
    civic_str = "-"
    if not pd.isnull(row["evidence_CGI"]):
        cgi_str = row["evidence_CGI"]

    if not pd.isnull(row["evidence_level_CIVIC"]):
        civic_str = row["evidence_level_CIVIC"]

    return f"{cgi_str}(cgi), {civic_str}(civic)"


def get_disease_names_CIViC(row):
    """
    Get CIViC's disease names for the variant

    :param row: Row of the variant dataframe
    :type row: pandas DataFrame row
    :return: Disease names for the variant
    :rtype: str
    """

    disease_rows = [
        row["evidence_phenotypes_CIVIC"],
        row["assertion_phenotypes_CIVIC"],
        row["evidence_disease_CIVIC"],
        row["assertion_disease_name_CIVIC"],
    ]
    disease_list = [i for i in disease_rows if not pd.isnull(i)]

    if disease_list != []:
        return ", ".join([i for i in disease_list])
    else:
        return ""


def get_therapy_names(row, civic_only):
    """
    Get the therapy names for the variant

    :param row: Row of the variant dataframe
    :type row: pandas DataFrame row
    """
    therapy_str = ""
    therapy_rows = [row["assertion_therapies_name_CIVIC"], row["evidence_therapies_CIVIC"]]
    therapy_list = [i.split(",") for i in therapy_rows if not pd.isnull(i)]
    if therapy_list != []:
        therapy_list = [i.strip() for i in flatten(therapy_list)]
        therapy_str = ", ".join([i for i in set(therapy_list)])
    if not pd.isnull(row["evidence_CGI"]) and not civic_only:
        if therapy_list != []:
            therapy_str = therapy_str + ", CGI (see Details)"
        else:
            therapy_str = "CGI (see Details)"
    else:
        return ""

    return therapy_str


def add_variant_name_report(df):
    """
    Adds a column with a name of the variant for the report to the df

    :param df: result df
    :type df: pandas DataFrame
    :return: List of variant names
    :rtype: list
    """
    result_list = []
    for index, row in df.iterrows():
        name_str = "chr{}-{}-{}-{}".format(row["chr_VEP"], row["pos_VEP"], row["ref_VEP"], row["alt_VEP"])
        count = 1
        while name_str in result_list:
            name_str = name_str + f"_{count}"
            count += 1
        result_list.append(name_str)

    return result_list


def create_link_col(row, report_path):
    """
    Creates a link to the individual report to display in the overall report

    :param row: row of the df
    :type row: pandas Series
    :param report_path: path to the directory in which the individual reports will be saved in
    :type report_path: str
    :return: link to the individual report
    :rtype: str"""
    return f'<a href="file://{report_path}/{row["report_name"]}.html">Variant: {row["report_name"]}</a>'


def save_plot(input, title, out_path):
    """
    Creates and saves a upsetplot figure as png

    :param input: input dataframe
    :type input: pandas.DataFrame
    :param title: title of the plot
    :type title: str
    :param out_path: output path
    :type out_path: str
    :return: matplotlib figure
    :rtype: matplotlib figure
    """

    fig = plt.figure(figsize=(6, 4))
    plot(input, subset_size="count", show_counts=True, facecolor="darkblue", element_size=None, fig=fig)
    plt.suptitle(title)
    plt.savefig(out_path)

    return fig


def create_upsetplots(df, out_path):
    """
    Create of (1) the number of variants per Knowledgebase and (2) the number of variants per tier

    :param df: Variant dataframe
    :type df: pandas.DataFrame
    :param out_path: Path to output directory
    :type out_path: str
    :return: List of Upsetplot figures
    :rtype: list
    """

    df_kb = pd.DataFrame()
    df_tiers = pd.DataFrame()

    # fill KB df & transform to upsetplot input format
    df_kb["cgi"] = df["Oncogenic Summary_CGI"].apply(lambda x: True if pd.isnull(x) == False else False)
    df_kb["civic"] = df["chr_CIVIC"].apply(lambda x: True if pd.isnull(x) == False else False)
    df_kb["none"] = ~(df_kb["cgi"] | df_kb["civic"])

    kb_input = from_indicators(df_kb)

    # fill tier df  & transform to upsetplot input format
    for i in ["tier_1", "tier_2", "tier_3", "tier_4"]:
        df_tiers[i] = df["report_tier"].apply(lambda x: True if x == i else False)

    tier_input = from_indicators(df_tiers)

    # create upsetplot figures
    fig_kb = save_plot(kb_input, "Number of variants per Knowledgebase", os.path.join(out_path, "kb_upsetplot.png"))
    fig_tiers = save_plot(tier_input, "Number of variants per Tier", os.path.join(out_path, "tier_upsetplot.png"))

    return [fig_kb, fig_tiers]


def encode_upsetplot(fig):
    """
    Encodes an upsetplot figure as base64

    :param fig: Upsetplot figure
    :type fig: matplotlib figure
    :return: Encoded upsetplot figure as string to add to report
    :rtype: str
    """
    tmpfile = BytesIO()
    fig.savefig(tmpfile, format="png")
    encoded = base64.b64encode(tmpfile.getvalue()).decode("utf-8")

    return "<img src='data:image/png;base64,{}'>".format(encoded)


def create_tier_table(df, tier, report_path):
    """
    Creates a table for the report for the given tier

    :param df: df containing all variants
    :type df: pandas DataFrame
    :param tier: tier to create the table for
    :type tier: int
    :param report_path: path to the directory in which the individual reports will be saved in
    :type report_path: str
    :return: df containing the variants of the given tier
    :rtype: pandas DataFrame
    """
    # extract tier from full df
    tier_x = df[df["report_tier"] == tier]

    # sort based on ranking_score
    tier_x = tier_x.sort_values("ranking_score", ascending=False)

    # assign "link col" to test
    tier_x["report_link"] = tier_x.apply(lambda x: create_link_col(x, report_path), axis=1)

    # create new df structure for report
    # source | HGVS | Oncogenic Summary | Disease Name | Therapy Name | Evidence_Level | link
    tier_report_subset = tier_x[
        [
            "sources",
            "Mutation_CGI",
            "Oncogenic Summary_CGI",
            "disease_names_report",
            "therapy_names_report",
            "evidence_levels_comb",
            "report_link",
        ]
    ]
    tier_report_subset = tier_report_subset.rename(
        columns={
            "sources": "Sources",
            "Mutation_CGI": "HGVS",
            "Oncogenic Summary_CGI": "Oncogenicity",
            "disease_names_report": "Disease Name",
            "therapy_names_report": "Therapies",
            "evidence_levels_comb": "Evidence Level",
            "report_link": "Details",
        }
    )
    html_table = build_table(
        tier_report_subset,
        color="blue_light",
        escape=False,
        width_dict=["100px", "300px", "200px", "200px", "200px", "100px", "auto"],
        text_align="center",
    )

    return html_table


def write_overall_report(template_html, report_html, fig_kb, fig_tiers, tier_table_list):
    """
    Create the overall report for the querynator results

    :param template_html: Path to the template html file
    :type template_html: str
    :param report_html: Path to the report html file
    :type report_html: str
    :param fig_kb: Path to upset plot of the kb distribution
    :type fig_kb: str
    :param fig_tiers: Path to upset plot of the tier distribution
    :type fig_tiers: str
    :param tier_table_list: List of pretty-html-tables of the tier tables
    :type tier_table_list: list
    :return: None
    :rtype: None
    """
    with open(template_html, "r") as f:
        with open(report_html, "w") as r:
            for line in f:
                if "TEMPLATE_1" in line:
                    line = line.replace("TEMPLATE_1", fig_kb)
                    line = line.replace("TEMPLATE_1", fig_tiers)
                    r.write(line)
                elif "TEMPLATE_2" in line:
                    line = line.replace("TEMPLATE_2", fig_tiers)
                    r.write(line)
                elif "TEMPLATE_TIER1" in line:
                    r.write(tier_table_list[0])
                elif "TEMPLATE_TIER2" in line:
                    r.write(tier_table_list[1])
                elif "TEMPLATE_TIER3" in line:
                    r.write(tier_table_list[2])
                else:
                    r.write(line)


# ============================================================================ #
#                               Individual Report
# ============================================================================ #


def check_if_nan(value):
    """
    Checks if a value is NaN and returns an empty string if it is.
    :param value: The value to check.
    :type value: str
    :return: The value if it is not NaN, otherwise an empty string.
    :rtype: str
    """
    if not pd.isnull(value):
        return value
    else:
        return "None"


def create_html_link(s):
    """
    Creates a clickable link from a string.
    Used to convert file path into clickable form.

    :param s: The string to create a link from.
    :type s: str
    :return: The string as a clickable link.
    :rtype: str
    """
    splitted = s.split("https")
    if "oncokb" in s and "CIVIC" in s:
        return s
    else:
        for i in splitted:
            if i.startswith("://"):
                path = f"https{i}"
                return f'<a href="{path}">{s}</a>'
        return s


def get_therapy_information_CGI(row, biomarkers_df, response, width_dict):
    """
    Gets all associated disease names of a specific variant.
    :param row: The row of the dataframe.
    :type row: pandas.Series
    :param biomarkers_df: The biomarkers dataframe from CGI.
    :type biomarkers_df: pandas.DataFrame
    :param response: The response to a specific drug.
    :type response: str
    :param width_dict: A list containing the width of each column.
    :type width_dict: list
    :return: A pandas DataFrame containing all Therapy & Drug related information provided by CGI for a specific Protein Change.
    :rtype: pandas.DataFrame
    """
    if not pd.isnull(row["Protein Change_CGI"]):
        alterations_df = biomarkers_df[biomarkers_df["alterations_link"].str.contains(row["Protein Change_CGI"])][
            ["Drugs", "Response", "Evidence", "Diseases", "Tumor type", "Source"]
        ]
        alterations_df.columns = alterations_df.columns = [
            "associated Diseases" if x == "Diseases" else x for x in alterations_df.columns
        ]
        alterations_df["Source"] = alterations_df["Source"].apply(lambda x: create_html_link(x))

        if response == "Responsive":
            return build_table(
                alterations_df[alterations_df["Response"] == "Responsive"].sort_values("Evidence", ascending=True),
                color="blue_light",
                escape=False,
                width_dict=width_dict,
            )

        else:
            return build_table(
                alterations_df[alterations_df["Response"] != "Responsive"].sort_values("Evidence", ascending=True),
                color="blue_light",
                escape=False,
                width_dict=width_dict,
            )

    else:
        return ""


def create_evidence_table(row):
    """
    Creates a table containing CIViC evidence information for a specific variant.
    :param row: The row of the dataframe.
    :type row: pandas.Series
    """
    evidence_subset = (
        row.loc[
            [
                "evidence_name_CIVIC",
                "evidence_type_CIVIC",
                "evidence_level_CIVIC",
                "evidence_rating_CIVIC",
                "evidence_significance_CIVIC",
                "evidence_support_CIVIC",
            ]
        ]
        .to_frame()
        .T
    )
    evidence_subset.columns = ["Name", "Level", "Rating", "Type", "Significance", "Direction"]
    if len(evidence_subset) > 0:
        return build_table(evidence_subset, color="blue_light", escape=False)
    else:
        return ""


def create_therapy_table(row, response, width_dict, biomarkers_df):
    """
    Creates a HTML table containing all Therapy & Drug related information provided by CGI for a specific Protein Change.
    :param row: The row of the dataframe.
    :type row: pandas.Series
    :param response: The response to a specific drug.
    :type response: str
    :param width_dict: A list containing the width of each column.
    :type width_dict: list
    :param biomarkers_df: The biomarkers dataframe from CGI.
    :type biomarkers_df: pandas.DataFrame
    :return: A HTML table containing all Therapy & Drug related information provided by CGI for a specific Protein Change.
    :rtype: HTML table
    """
    return build_table(
        get_therapy_information_CGI(row, biomarkers_df, response),
        color="blue_light",
        escape=False,
        width_dict=width_dict,
    )


def get_reference_build(metadata_path):
    """
    Gets the reference build from the metadata file.
    :param metadata_path: The path to the metadata file.
    :type metadata_path: str
    :return: The reference build.
    :rtype: str
    """
    with open(metadata_path, "r") as f:
        for line in f:
            if line.startswith("Reference genome"):
                return line.split(":")[1].strip()
        return ""


def retrieve_info_from_row(row, biomarkers_df, metadata_path):
    """
    This function retrieves the information from a row of the merged dataframe
    and returns a dictionary with the information for the report of a specific variant.

    :param row: row of the dataframe
    :type row: pandas.core.series.Series
    :param biomarkers_df: dataframe with all biomarkers linked to a specific variant
    :type biomarkers_df: pandas.core.frame.DataFrame
    :param metadata_path: The path to the metadata file.
    :type metadata_path: str
    :return: dictionary with the information for the report of a specific variant
    :rtype: dict
    """
    info_dict = {}
    info_dict["VARIANT_NAME"] = check_if_nan(row["report_name"])
    info_dict["TIER_ASSIGNED"] = check_if_nan(row["report_tier"])
    info_dict["CHROMOSOME"] = check_if_nan(row["chr_CGI"])
    info_dict["POSITION"] = check_if_nan(row["pos_CGI"])
    info_dict["REF"] = check_if_nan(row["ref_CGI"])
    info_dict["ALT"] = check_if_nan(row["alt_CGI"])
    info_dict["STRAND"] = check_if_nan(row["Strand_CGI"])
    info_dict["HGVS_NOTATION"] = check_if_nan(row["Mutation_CGI"])
    info_dict["VARIANT_CIVIC"] = check_if_nan(row["variant_name_CIVIC"])
    info_dict["ENTREZ_NAME"] = check_if_nan(row["variant_entrez_name_CIVIC"])
    info_dict["ENTREZ_ID"] = (
        int(check_if_nan(row["variant_entrez_id_CIVIC"]))
        if check_if_nan(row["variant_entrez_id_CIVIC"]) != "None"
        else check_if_nan(row["variant_entrez_id_CIVIC"])
    )
    info_dict["VARIANT_TYPE"] = check_if_nan(row["Type_CGI"])
    info_dict["CIVIC_GENE"] = check_if_nan(row["gene_name_CIVIC"])
    info_dict["CGI_GENE"] = check_if_nan(row["Gene_CGI"])
    info_dict["GENE_DESCRIPTION"] = check_if_nan(row["gene_description_CIVIC"])
    info_dict["BUILD"] = get_reference_build(metadata_path) if get_reference_build(metadata_path) != "" else "None"
    info_dict["SIFT_SCORE"] = check_if_nan(row["SIFT_VEP"])
    info_dict["PolyPhen2_SCORE"] = check_if_nan(row["PolyPhen_VEP"])
    info_dict["ALLELE_FREQ"] = check_if_nan(row["AF_VEP"])
    info_dict["GNOMAD_FREQ"] = check_if_nan(row["gnomAD_AF_VEP"])
    info_dict["CONSEQUENCE_CGI"] = check_if_nan(row["Consequence_CGI"])
    info_dict["CONSEQUENCE_CIVIC"] = check_if_nan(row["variant_type_CIVIC"])
    info_dict["ONCOGENICITY_CGI"] = check_if_nan(row["Oncogenic Summary_CGI"])
    info_dict["ONCOGENIC_PREDICTION"] = check_if_nan(row["Oncogenic Prediction_CGI"])
    info_dict["CLINVAR_ENTRIES"] = check_if_nan(row["variant_clinvar_entries_CIVIC"])
    info_dict["PROT_CHANGE"] = check_if_nan(row["Protein Change_CGI"])
    info_dict["CIVIC_CGI_PRESENT"] = check_if_nan(row["sources"])
    info_dict["EXT_ANNOS"] = check_if_nan(row["External oncogenic annotation_CGI"])
    info_dict["LINKED_DISEASES"] = get_disease_names_CIViC(row) if get_disease_names_CIViC(row) != "" else "None"
    info_dict["LINKED_THERAPIES_CIVIC"] = (
        get_therapy_names(row, civic_only=True) if get_therapy_names(row, civic_only=True) != "" else "None"
    )
    info_dict["EVIDENCE_LIST"] = create_evidence_table(row) if create_evidence_table(row) != "" else "None"
    info_dict["CIVIC_SUMMARY"] = check_if_nan(row["assertion_summary_CIVIC"])
    info_dict["ASSERTION_DESCRIPTION"] = check_if_nan(row["assertion_description_CIVIC"])
    info_dict["EVIDENCE_DESCRIPTION"] = check_if_nan(row["evidence_description_CIVIC"])
    info_dict["LINKED_DRUGS_RESPONSIVE"] = (
        get_therapy_information_CGI(row, biomarkers_df, "Responsive", ["40%", "20%", "20%", "20%", "20%", "20%"])
        if get_therapy_information_CGI(row, biomarkers_df, "Responsive", ["40%", "20%", "20%", "20%", "20%", "20%"])
        != ""
        else "None"
    )
    info_dict["LINKED_DRUGS_RESISTANT"] = (
        get_therapy_information_CGI(row, biomarkers_df, "Resistant", ["40%", "20%", "20%", "20%", "20%", "20%"])
        if get_therapy_information_CGI(row, biomarkers_df, "Resistant", ["40%", "20%", "20%", "20%", "20%", "20%"])
        != ""
        else "None"
    )
    info_dict["EVIDENCE_DESCRIPTION"] = check_if_nan(row["evidence_description_CIVIC"])
    info_dict["LITERATURE_CIVIC"] = check_if_nan(row["evidence_source_CIVIC"])

    return info_dict


def write_individual_report(row, template_html, report_path, biomarkers_df, metadata_path):
    """
    This function creates a report for a specific variant.

    :param row: row of the dataframe
    :type row: pandas.core.series.Series
    :param template_html: path to the template html file
    :type template_html: str
    :param report_path: path to the individual report directory
    :type report_path: str
    :param biomarkers_df: dataframe with the biomarkers linked to the therapies
    :type biomarkers_df: pandas.core.frame.DataFrame
    :param metadata_path: path to the metadata file
    :type metadata_path: str
    :return: None
    :rtype: None
    """
    info_dict = retrieve_info_from_row(row, biomarkers_df, metadata_path)

    report_html = "{}/{}.html".format(report_path, row["report_name"])

    with open(template_html, "r") as f:
        with open(report_html, "w") as r:
            for line in f:
                for i in info_dict:
                    if i in line:
                        line = line.replace(i, str(info_dict[i]))
                r.write(line)


# ============================================================================ #
#                               Build the Report
# ============================================================================ #


def create_report_htmls(outdir, basename, civic_path, logger):
    """
    Creates the overall report for all variants and
    the individual reports for each variant

    :param outdir: Path to report directory
    :type outdir: str
    :param basename: User given Project name
    :type basename: str
    :civic_path: Path to civic results
    :type civic_path: str
    :return: None
    :rtype: None
    """

    logger.info("Creating HTML Reports")

    # read in files
    vep_civic_cgi_merge = pd.read_csv(f"{outdir}/combined_files/civic_cgi_vep.tsv", sep="\t")
    biomarkers_df = pd.read_csv(f"{outdir}/combined_files/biomarkers_linked.tsv", sep="\t")
    metadata_civic = f"{civic_path}/metadata.txt"  # read reference genome from metadata file
    # get path to save individual reports
    report_path = f"{os.path.abspath(outdir)}/report/variant_reports"

    # add additional columns used in report
    vep_civic_cgi_merge["sources"] = vep_civic_cgi_merge.apply(lambda x: get_sources(x), axis=1)
    vep_civic_cgi_merge["report_name"] = add_variant_name_report(vep_civic_cgi_merge)
    vep_civic_cgi_merge["disease_names_report"] = vep_civic_cgi_merge.apply(
        lambda x: get_disease_names_CIViC(x), axis=1
    )
    vep_civic_cgi_merge["therapy_names_report"] = vep_civic_cgi_merge.apply(
        lambda x: get_therapy_names(x, civic_only=False), axis=1
    )
    vep_civic_cgi_merge["evidence_levels_comb"] = vep_civic_cgi_merge.apply(
        lambda x: assign_comb_evidence_labels(x), axis=1
    )
    vep_civic_cgi_merge["report_name"] = add_variant_name_report(vep_civic_cgi_merge)

    # create "pretty_html_tables" for each tier
    tier_list = ["tier_1", "tier_2", "tier_3"]
    tier_table_list = [create_tier_table(vep_civic_cgi_merge, i, report_path) for i in tier_list]

    # create the upset plots
    fig_kb, fig_tiers = create_upsetplots(vep_civic_cgi_merge, os.path.join(os.path.abspath(outdir), "report/plots"))

    # Overall Report
    write_overall_report(
        template_html=os.path.join(os.path.dirname(__file__), "templates/template_overall_upsetplots.html"),
        report_html=f"{os.path.abspath(outdir)}/report/{basename}_overall_report.html",
        fig_kb=encode_upsetplot(fig_kb),
        fig_tiers=encode_upsetplot(fig_tiers),
        tier_table_list=tier_table_list,
    )

    # Individual Reports

    # we only create individual reports for variants in tier 1,2 and 3
    vep_civic_cgi_merge = vep_civic_cgi_merge[vep_civic_cgi_merge["report_tier"].isin(["tier_1", "tier_2", "tier_3"])]

    vep_civic_cgi_merge.apply(
        lambda x: write_individual_report(
            x,
            template_html=os.path.join(os.path.dirname(__file__), "templates/template_individual.html"),
            report_path=report_path,
            biomarkers_df=biomarkers_df,
            metadata_path=metadata_civic,
        ),
        axis=1,
    )

    logger.info("Reports created")
