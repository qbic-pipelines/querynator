Changelog
============

0.6.0 - Keppler-452b Goldilocks (2024-10-10)
---------------------------------------------

**Added**

* option `--filter_evidence` to filter CIViC evidences based on type, direction, status, level and significance

**Fixed**

**Dependencies**

* updated `black` to version 24.3.0 due to moderate vulnerability CVE-2024-21503
* updated `requests` to version 2.31.0 due to moderate vulnerabilities CVE-2023-32681 and CVE-2024-35195

**Deprecated**

0.5.5 - Sulfur Io  (2024-09-09)
---------------------------------------------

**Added**

**Fixed**

* [#52](https://github.com/qbic-pipelines/querynator/issues/52): issue that lead to an inconsitent number of fields for the CIViC evidences
* [#51](https://github.com/qbic-pipelines/querynator/issues/51): CGI evidences are now filtered by the specified cancer type

**Dependencies**

**Deprecated**

0.5.4 - Sulfur Io  (2024-07-24)
---------------------------------------------

**Added**

**Fixed**

* [#43](https://github.com/qbic-pipelines/querynator/issues/43): issue with disease ontology file path, which resurfaced in docker container
* evidences without specified disease are now excluded if `--cancer` is specified

**Dependencies**

**Deprecated**


0.5.3 - Sulfur Io  (2024-07-20)
---------------------------------------------

**Added**

**Fixed**

* [#43](https://github.com/qbic-pipelines/querynator/issues/43): added disease ontology file to the package

**Dependencies**

**Deprecated**


0.5.2 - Sulfur Io  (2024-07-05)
---------------------------------------------

**Added**

**Fixed**

* [#37](https://github.com/qbic-pipelines/querynator/issues/37): fixed bug in `query-api-civic` that caused the containerized tool to fail when querying CIViC.

**Dependencies**

* pinned `pandas` to version 1.5.3 to fix deprecation issues with `pandas 2.x` in `querynator` and `civicpy`.

**Deprecated**


0.5.1 - Sulfur Io  (2024-07-03)
---------------------------------------------

**Added**

**Fixed**

* fixed a severe bug in `query-api-civic` that made the tool unusable, the bug was introduced in v0.5.0.

**Dependencies**

**Deprecated**


0.5.0 - Sulfur Io  (2024-06-23)
---------------------------------------------

**Added**

* CL option to specify cancer type when querying CIViC (`-c` `--cancer`)

**Fixed**

**Dependencies**

* pinned numpy to version 1.24.4 to fix pandas bug with numpy-v2

**Deprecated**


0.4.2 - Stormy Saturn  (2024-04-05)
---------------------------------------------

**Added**

**Fixed**

* Bug fixes to handle new CGI wildtype biomarkers
* civicpy cache will only be updated when civic is queried

**Dependencies**

**Deprecated**


0.4.1 - Stormy Saturn  (2023-06-13)
---------------------------------------------

**Added**

**Fixed**

* Bug fixes to include all evidence of CIViC

**Dependencies**

**Deprecated**

0.4.0 - Stormy Saturn  (2023-05-31)
---------------------------------------------

**Added**

**Fixed**

* Fixed functionality for new CGI file structure
* Fixed case when CIViC has no hits

**Dependencies**

**Deprecated**

0.3.3 - Iron Mercury  (2023-05-05)
---------------------------------------------

**Added**

**Fixed**

* Fixed API docs

**Dependencies**

**Deprecated**

0.3.2 - Iron Mercury  (2023-05-05)
---------------------------------------------

**Added**

**Fixed**

* Fixed version bump

**Dependencies**

**Deprecated**

0.3.1 - Iron Mercury  (2023-05-05)
---------------------------------------------

**Added**

**Fixed**

* Fixed import of site-packages in setup.py

**Dependencies**

**Deprecated**

0.3.0 - Iron Mercury  (2023-05-04)
---------------------------------------------

**Added**

* Added functionality to combine the results of the Knowledgebases in an HTML report
* Added possibility to have non-numerical chromosome columns in the input vcf
* Added deletion of CGI jobs from CGI Server after completion

**Fixed**

**Dependencies**

**Deprecated**

0.2.2 - Sour Venus  (2023-03-16)
---------------------------------------------

**Added**

* Optional VEP annotation based filtering
* Additional metadata
* Usage of pyVCF3 to read vcf files
* Querynator ID added for filtered vcf files
* All possible reference genomes for CIViC

**Fixed**

**Dependencies**

**Deprecated**

* Usage of pysam to read vcf files


0.2.1 - Sour Venus  (2023-02-16)
---------------------------------------------

**Added**

**Fixed**

* Rendering API docs

**Dependencies**

**Deprecated**

0.2.0 - Sour Venus  (2023-02-07)
---------------------------------------------

**Added**

* Added functionality to query the Clinical Interpretation of Variants in Cancer (CIViC) Knowledgebase
* Added possibility to query bgzipped files

**Fixed**

**Dependencies**

**Deprecated**

0.1.3 - Diamond Neptune  (2022-11-21)
---------------------------------------------

**Added**

**Fixed**

* Fix including module

**Dependencies**

**Deprecated**

0.1.2 - Diamond Neptune  (2022-11-18)
---------------------------------------------

**Added**

**Fixed**

* Fix installing requirements

**Dependencies**

**Deprecated**

0.1.1 -  Methane Titan (2022-11-18)
---------------------------------------------

**Added**

**Fixed**

* Github Actions publishing to PyPI
* Fix docs

**Dependencies**

**Deprecated**


0.1.0 - initial release (2022-11-18)
---------------------------------------------

**Added**

* First release on PyPI
* Created the package template with cookiecutter
* Functions to query the cancergenomeinterpreter REST API

**Fixed**

**Dependencies**

**Deprecated**
