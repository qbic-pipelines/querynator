.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/qbic-pipelines/querynator/issues.
Please check beforehand that there isn't `already an issue <https://github.com/qbic-pipelines/querynator/issues>`_.
about your idea to avoid duplicating work.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it. Please also check beforehand if there isn't already an `open PR <https://github.com/qbic-pipelines/querynator/pulls>`_.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

querynator could always use more documentation, whether as part of the
official querynator docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/qbic-pipelines/querynator/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up ``querynator`` for local development.

1. Fork the ``querynator`` repo on GitHub to your account.
2. Clone your fork locally:

.. code-block:: bash

    git clone git@github.com:your_name_here/querynator.git

3. Create a virtualenv or conda environment and install the required packages and the querynator for development:

.. code-block:: bash

    pip install --upgrade -r requirements-dev.txt
    pip install -e .

4. Create a branch for local development:

.. code-block:: bash

    git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.

5. When you're done making changes, check that your changes pass pytest:

.. code-block:: bash

    python3 -m pytest tests/ --color=yes

6. Docstrings should adhere to the `Sphinx docu <https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html>`_. These are used to automatically generate package documentation on the querynator website using Sphinx.
You can find this documentation `here <https://querynator.readthedocs.io/en/latest/index.html>`_.

If you would like to test the documentation, you can install Sphinx locally by following Sphinx's `installation instruction <https://www.sphinx-doc.org/en/master/usage/installation.html>`_.
Once done, you can run ``make -C docs html`` in the directory of ``qbic-pipelines/querynator``.
The HTML will then be generated in ``docs/api/_build/html`` and can be opened in any web browser.

7. All Python code in qbic-pipelines/querynator must be passed through the `Black Python code formatter <https://black.readthedocs.io/en/stable/>`_
through `isort <https://pycqa.github.io/isort/index.html>`_ and `prettier <https://prettier.io/>`_.
This ensures a harmonized code formatting style and harmonized imports throughout the package, from all contributors.
You can run Black, isort and prettier on the command-line for the entire repository.

.. code-block:: bash

  black .
  isort .
  prettier --write .

8. Commit your changes and push your branch to GitHub:

.. code-block:: bash

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

9. Submit a pull request against the ``dev`` branch through the GitHub website and wait for the code to be reviewed and merged.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a docstring, and add the feature to the list in README.rst.
3. The pull request should work for Python 3.8, 3.9, and 3.10 and for PyPy.


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in CHANGELOG.rst).
Then run:

.. code-block:: bash

  bump2version patch # possible: major / minor / patch
  git push
  git push --tags

Travis will then deploy to PyPI if tests pass.
