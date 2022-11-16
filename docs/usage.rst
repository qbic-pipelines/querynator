=====
Usage
=====

The ``querynator`` is used from the command line:

.. code-block:: bash

    querynator --help
    querynator query-api-cgi --help

A typical command to query the `cancergenomeinterpeter - CGI <https://www.cancergenomeinterpreter.org/home>`_:
    
.. code-block:: bash

    querynator query-api-cgi \
        -i variants.tsv \
        -o sample_name \
        -g hg38 \
        -c 'Liver & billiary tract' \
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

