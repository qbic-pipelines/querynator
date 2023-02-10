=====
Usage
=====

The ``querynator`` is used from the command line. Use the help function to display all options and choices.

.. code-block:: bash

    querynator --help
    querynator query-api-cgi --help
    querynator query-api-civic --help


Query the cancergenomeinterpeter - CGI
**************************************

A typical command to query the `cancergenomeinterpeter - CGI <https://www.cancergenomeinterpreter.org/home>`_:

.. code-block:: bash

    querynator query-api-cgi \
        -m mutations.[vcf,tsv,gtf] \
        -o sample_name \
        -g hg38 \
        -c 'Any cancer type' \
        --email your-cgi-account-mail@whatever.com \
        --token your-cgi-token

The command above downloads the following files from CGI.
For further information please refer to their `FAQ <https://www.cancergenomeinterpreter.org/faq#q18>`_.

.. code-block:: bash

    sample_name.cgi_results
    ├── drug_prescription.tsv
    ├── input01.tsv
    ├── metadata.txt
    └── mutation_analysis.tsv


Mutation, CNA & translocation analysis
======================================

If you run the command with all possible input files, you will obtain:

.. code-block:: bash

    # run command
    querynator query-api-cgi \
        --mutations mutations.[vcf,tsv,gtf] \
        --cnas cnas.tsv \
        --translocations translocations.tsv \
        -o sample_name \
        -g hg38 \
        -c 'Any cancer type' \
        --email your-cgi-account-mail@whatever.com \
        --token your-cgi-token

    # output
    sample_name.cgi_results
    ├── cna_analysis.tsv
    ├── drug_prescription.tsv
    ├── fusion_analysis.tsv
    ├── input01.tsv
    ├── input02.tsv
    ├── input03.tsv
    ├── metadata.txt
    └── mutation_analysis.tsv


Input file formats
==================

For detailed information please refer to `CGI formats <https://www.cancergenomeinterpreter.org/faq#q22>`_.
The genomic tabular format ``gtf`` is displayed below and contains partly the same columns as a ``vcf`` file (>v. 4.0) and is tab-separated.
The `sample column` is not mandatory, but recommended when more than one sample is contained in one file.

A mutations/variant file can have the extensions ``vcf``, ``tsv`` or ``gtf``. The column names can also be uppercase letters as in a ``vcf``.

.. list-table:: mutations.[vcf,tsv,gtf]
    :widths: 25 25 25 25 25
    :header-rows: 1

    *   - sample
        - #chrom/chr
        - pos
        - ref
        - alt
    *   - test1
        - chr4
        - 121369475
        - A
        - T
    *   - test2
        - chr10
        - 122630837
        - C
        - G


A copy number alterations file should be ``tsv`` and column names must be lowercase.

.. list-table:: cnas.tsv
    :widths: 25 25 25
    :header-rows: 1

    *   - sample
        - gene
        - cna
    *   - test1
        - ERBB2
        - amp
    *   - test2
        - TP53
        - del

A copy number alterations file should be ``tsv`` and column names must be lowercase.

.. list-table:: translocations.tsv
    :widths: 25 25
    :header-rows: 1

    *   - sample
        - fus
    *   - test1
        - BCR__ABL1
    *   - test2
        - PML__RARA


Genome build versions
=====================

.. note::
    The cancergenomeinterpeter will perform a liftover of the genomic coordinates to `hg38` if the parameter ``--genome hg19`` is used.


Query the Clinical Interpretations of Variants in Cancer - CIViC
****************************************************************

A typical command to query the `Clinical Interpretations of Variants in Cancer - CIViC <https://civicdb.org/welcome>`_:

.. code-block:: bash

    querynator query-api-civic \
        -v input_file.vcf,tsv,gtf \
        -o outdir \

The command above generates the following result files using `CIViCpy <https://docs.civicpy.org/>`_.

.. code-block:: bash

    outdir
    ├── civic_results.tsv
    └── metadata.txt


Input file format
==================

The querynator requires a ``vcf`` file in uncompressed or in `bgzipped format <http://www.htslib.org/doc/bgzip.html>`_ ``vcf.gz`` to query CIViC.

It is recommended (although not required) to provide an index-file (``vcf.gz.tbi``) with the input ``vcf`` file, e.g. using `tabix <http://www.htslib.org/doc/tabix.html>`_.
The index file must be stored in the same directory as the ``vcf`` file.
