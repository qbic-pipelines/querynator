=====
Usage
=====

The ``querynator`` is used from the command line. Use the help function to display all options and choices.

.. code-block:: bash

    querynator --help
    querynator query-api-cgi --help
    querynator query-api-civic --help

VCF Normalization
**************************************
We recommend to normalize the input vcf's using `bcftools norm <https://samtools.github.io/bcftools/bcftools.html>`_ before
running the querynator to unify the input:

The vcf file must be indexed to run ``bcftools norm``, e.g. using `tabix <http://www.htslib.org/doc/tabix.html>`_:

.. code-block:: bash

    tabix /path/to/vcf_file.vcf

Then run ``bcftools norm`` in the following way:

.. code-block:: bash

    bcftools norm \
        -a \
        --atom-overlaps . \
        -f /path/to/reference_file \
        /path/to/vcf




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

    sample_name
    ├── sample_name.cgi_results
    |   ├── drug_prescription.tsv
    |   ├── input01.tsv
    |   ├── metadata.txt
    |   └── mutation_analysis.tsv

.. note::
    The input variants for the CGI query must be sorted based on their coordinates.

.. note::
    Too many requests in too short amount of time can result in the email address used being blocked.

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
    sample_name
    ├── sample_name.cgi_results
    |   ├── cna_analysis.tsv
    |   ├── drug_prescription.tsv
    |   ├── fusion_analysis.tsv
    |   ├── input01.tsv
    |   ├── input02.tsv
    |   ├── input03.tsv
    |   ├── metadata.txt
    |   └── mutation_analysis.tsv
    └── sample_name.cgi_results.zip


Using the ``filter_vep`` `flag <https://querynator.readthedocs.io/en/latest/usage.html#filtering-benign-variants>`_, the querynator can filter out benign variants in ``vcf`` files before querying the knowledgebase (KB).


Input file formats
==================

For detailed information please refer to `CGI formats <https://www.cancergenomeinterpreter.org/faq#q22>`_.
The genomic tabular format ``gtf`` is displayed below and contains partly the same columns as a ``vcf`` file (>v. 4.0) and is tab-separated.
The `sample column` is not mandatory, but recommended when more than one sample is contained in one file.

A mutations/variant file can have the extensions ``vcf``, ``vcf.gz``, ``tsv`` or ``gtf``. The column names can also be uppercase letters as in a ``vcf``.

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

A translocation file should be ``tsv`` and column names must be lowercase.

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
        -v input_file.vcf \
        -o outdir \
        -g ref_genome [GRCh37, GRCh38, NCBI36]

The command above generates the following result files using `CIViCpy <https://docs.civicpy.org/>`_.

.. code-block:: bash

    sample_name
    ├── sample_name.civic_results.tsv
    └── metadata.txt

The querynator performs an ``exact`` search, meaning that variants in the KB must match the given coordinates, reference allele(s) and alternate allele(s) precisely.

Using the ``filter_vep`` `flag <https://querynator.readthedocs.io/en/latest/usage.html#filtering-benign-variants>`_, the querynator can filter out benign variants in ``vcf`` files before querying the KB.

Input file format
==================

The querynator requires a ``vcf`` file (>v. 4.0) in uncompressed or in `bgzipped format <http://www.htslib.org/doc/bgzip.html>`_ ``vcf.gz`` to query CIViC.

It is recommended (although not required) to provide an index-file (``vcf.gz.tbi``) with the input ``vcf`` file, e.g. using `tabix <http://www.htslib.org/doc/tabix.html>`_.
The index file must be stored in the same directory as the ``vcf`` file.


Filtering benign variants
****************************************************************

Variants that are classified as ``low Impact`` and ``synonymous variants`` will be filtered out based on their `Ensembl VEP
annotation <https://www.ensembl.org/info/docs/tools/vep/index.html>`_ if the additional flag ``filter_vep`` is set.
The filtering step can be applied before querying both KBs.
Currently filtering can only be applied on VEP annotated ``vcf`` files. In order to filter the file,
the querynator expects a ``vcf`` that was annotated using VEP's standard key (``CSQ``).

To filter, the following fields are required in the VEP info column:

- Consequence
- IMPACT

If ``filter_vep`` is set, the filtered and removed variants are given out as results in the ``vcf_files`` directory.

A typical command for a CIViC query:

.. code-block:: bash

    querynator query-api-civic \
        -v input_file.vcf,tsv,gtf \
        -o outdir \
        -g ref_genome [GRCh37, GRCh38, NCBI36] \
        --filter_vep

The command above generates the following result files using `CIViCpy <https://docs.civicpy.org/>`_.

.. code-block:: bash

    sample_name
    ├── vcf_files
    |   ├── sample_name.filtered_variants.vcf
    |   ├── sample_name.removed_variants.vcf
    ├── sample_name.civic_results.tsv
    └── metadata.txt

.. note::
    When the ``filter_vep`` flag is set a unique Querynator ID is added to the INFO column of each variant in the vcf file.
    The same ID is added to the ``sample_name.civic_results.tsv`` if CIViC is queried.


Create an HTML Report
**************************************

After querying the knowledgebases included in the querynator, it is possible to combine the results into one table
and to create an HTML report summarizing the most important features of each variant.


.. note::
    This functionality was specifically created to be included into the ``variantMTB`` nextflow pipeline
    which can be found `here`_. (Currently under development)

.. _here: https://github.com/qbic-pipelines/variantmtb


A typical command to create such a report:

.. code-block:: bash

    querynator create-report \
        --cgi_path path/to/cgi_results \
        --civic_path path/to/civic_results \
        --outdir path/to/save/results

The command above generates the following result directory:

.. code-block:: bash

    outdir
    ├── combined_files
    |   ├── alterations_vep.tsv
    |   ├── biomarkers_linked.tsv
    |   ├── civic_cgi_vep.tsv
    |   └── civic_vep.tsv
    ├── report
    |   |   ├── overall_report.html
    |   |   ├── variant_reports
    |   |   |    ├── chr1-1234-A-T.html
    |   |   |    ├── chr1-1245-C-G.html
    |   |   |    └── ...
    |   |   ├── plots
    |   |   |    ├── kb_upsetplot.png
    └── └── └──  └── tier_upsetplot.png



The command creates one overall report which includes some statistics and shows an overview of the most important variants in the project.
The ``Details`` column in the overall report links directly to a more detailed report on the variant in question.
