""" Functions that are used by multiple scripts bundled """

import gzip
import os
import shutil


def flatten(l_l):
    """
    flattens a list of lists consisting of strings

    :param l_l: list of lists
    :type values: list
    :return: flattened list
    :rtype: list
    """
    flattened_list = []
    for i in l_l:
        if isinstance(i, list):
            flattened_list.extend(flatten(i))
        else:
            flattened_list.append(i)
    return flattened_list


def gzipped(file_path):
    """
    Helper function to test if given vcf is gzipped.
    If so, the first 2 bytes are "1f 8b"

    :param file_path: Path to gzipped input file
    :type file_path: str
    """
    with open(file_path, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


def gunzip_compressed_files(file_path, logger):
    """
    gunzips gzipped vcf file

    :param file_path: Path to gzipped input file
    :type file_path: str
    """
    logger.info(f"Unzipping input file ({os.path.basename(os.path.normpath(file_path))})")

    if not file_path.endswith(".gz"):
        logger.error("Given file does not end with '.gz'")
        exit(1)
    else:
        with gzip.open(file_path, "rb") as f_in:
            with open(file_path[: -len(".gz")], "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return file_path[: -len(".gz")]


def get_num_from_chr(s):
    """
    extracts numerical value from chromosome (chr1 -> 1)
    :param s: chromosome of variant
    :type vcf_path: pyVCF3 CHROM object
    :return: numerical value from variant's chromosome
    :rtype: str
    """
    int_chr = s.split("chr")[1] if str(s).startswith("chr") else s
    return str(int_chr)
