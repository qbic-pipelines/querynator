==========
querynator
==========


.. image:: https://img.shields.io/pypi/v/querynator.svg
        :target: https://pypi.python.org/pypi/querynator

.. image:: https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fqbic-pipelines%2Fquerynator%2Fbadge%3Fref%3Dmaster&style=flat-square
        :target: https://actions-badge.atrox.dev/qbic-pipelines/querynator/goto?ref=master

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
        :target: https://github.com/psf/black

.. image:: https://img.shields.io/badge/code%20style-prettier-ff69b4.svg
        :target: https://github.com/prettier/prettier

.. image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
        :target: https://pycqa.github.io/isort

.. image:: https://readthedocs.org/projects/querynator/badge/?version=latest
        :target: https://querynator.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Python package to query cancer variant databases


* Free software: MIT license
* Documentation: https://querynator.readthedocs.io.


Features
--------

* Command-line tool to query the `cancergenomeinterpreter <https://www.cancergenomeinterpreter.org/home>`_ via its REST API
* Command-line tool to query the `Clinical Interpretation of Variants in Cancer (CIViC) Knowledgebase <https://civicdb.org/>`_ using `CIViCpy <https://docs.civicpy.org/en/latest/>`_

Credits
-------

This package uses the cancergenomeinterpreter.org REST API for data retrieval.

* Muiños, F., Martínez-Jiménez, F., Pich, O. et al. In silico saturation mutagenesis of cancer genes. Nature 596, 428–432 (2021). https://doi.org/10.1038/s41586-021-03771-1
* Tamborero, D. Rubio-Perez, C., Deu-Pons, J. et al., Cancer Genome Interpreter annotates the biological and clinical relevance of tumor alterations. Genome Medicine 10, (2018). doi: https://doi.org/10.1101/140475

This package uses the CIViCpy package for data retrieval from the CIViC database.

* Wagner, Alex H., et al. "CIViCpy: a python software development and analysis toolkit for the CIViC knowledgebase." JCO Clinical Cancer Informatics 4 (2020): 245-253. doi: https://doi.org/10.1200/CCI.19.00127
* Griffith, M., Spies, N., Krysiak, K. et al. CIViC is a community knowledgebase for expert crowdsourcing the clinical interpretation of variants in cancer. Nat Genet 49, 170–174 (2017). doi: https://doi.org/10.1038/ng.3774

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


