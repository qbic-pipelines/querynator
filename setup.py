#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("CHANGELOG.rst") as changelog_file:
    changelog = changelog_file.read()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

test_requirements = [
    "pytest>=3",
]

setup(
    author="Susanne Jodoin, Mark Polster",
    author_email="susanne.jodoin@qbic.uni-tuebingen.de, mark.polster@uni-tuebingen.de",
    python_requires=">=3.8,<3.11",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    description="Python package to query cancer variant databases",
    entry_points={
        "console_scripts": [
            "querynator=querynator.__main__:run_querynator",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + changelog,
    include_package_data=True,
    keywords="querynator",
    name="querynator",
    packages=find_packages(exclude=("docs")),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/qbic-pipelines/querynator",
    version="0.3.0",
    zip_safe=False,
)
