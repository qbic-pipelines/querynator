#!/usr/bin/env python

"""Tests for `querynator` package."""

import unittest
from click.testing import CliRunner

import querynator.__main__
from querynator.helper_functions import ontology


class testCLI(unittest.TestCase):

    def test_commandLineInterface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(querynator.__main__.querynator_cli)
        assert result.exit_code == 0

        """Test help message"""
        help_result = runner.invoke(querynator.__main__.querynator_cli, ["--help"])
        assert help_result.exit_code == 0
        assert "--help     Show this message and exit." in help_result.output

        """Test CGI help message"""
        cgi_help_result = runner.invoke(querynator.__main__.querynator_cli, ["query-api-cgi", "--help"])
        assert cgi_help_result.exit_code == 0
        assert "--help                          Show this message and exit." in cgi_help_result.output

        """Test CIVIC help message"""
        civic_help_result = runner.invoke(querynator.__main__.querynator_cli, ["query-api-civic", "--help"])
        assert civic_help_result.exit_code == 0
        assert "--help                          Show this message and exit." in civic_help_result.output

        """Test CREATE-REPORT help message"""
        civic_help_result = runner.invoke(querynator.__main__.querynator_cli, ["create-report", "--help"])
        assert civic_help_result.exit_code == 0
        assert "--help                 Show this message and exit." in civic_help_result.output

        """Test non-existing subcommand"""
        result = runner.invoke(querynator.__main__.querynator_cli, ["query-api-clinvar"])
        assert result.exit_code == 2
        assert "No such command" in result.output

        """Test non-existing option"""
        result = runner.invoke(querynator.__main__.querynator_cli, ["query-api-cgi", "--baz"])
        assert result.exit_code == 2
        assert "No such option" in result.output


class testOntology(unittest.TestCase):
    """Test the ontology classes"""

    def test_disease_ontology(self):
        """Test the Disease Ontology class"""
        doid = ontology.DO("querynator/helper_functions/doid.obo")
        print(doid)

        self.assertIn("DOID:4", doid.ids)
        self.assertIn("DOID:4", doid.terms)
        self.assertIn("DOID:4", doid)
        self.assertIn("DOID:4947", doid)

        term = doid.get("DOID:4947")
        self.assertEqual(term, doid.get(4947))
        self.assertEqual(term, doid.get("cholangiocarcinoma"))
        self.assertEqual(term, doid.get("Cholangiocarcinoma"))

        self.assertIn(doid.get("bile duct adenocarcinoma"), doid.get("cholangiocarcinoma").relationships)
        self.assertIn(doid.get("bile duct adenocarcinoma"), doid.get_all_ancestors("cholangiocarcinoma"))
        self.assertIn(doid.get("cholangiocarcinoma"), doid.get_all_ancestors("cholangiocarcinoma", includeSelf=True))

        self.assertIsNone(doid.get("DOID:123456"))
        self.assertIsNone(doid.get_from_name("foo"))


if __name__ == "__main__":
    unittest.main()
