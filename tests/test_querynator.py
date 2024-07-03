#!/usr/bin/env python

"""Tests for `querynator` package."""

import os
import shutil
import unittest

from click.testing import CliRunner

from querynator.__main__ import querynator_cli


class testCliHelp(unittest.TestCase):

    def setUp(self):
        self.runner = CliRunner()

    def test_commandLineInterface(self):
        """Test the CLI."""
        result = self.runner.invoke(querynator_cli)
        self.assertEqual(result.exit_code, 0)

    def test_helpMessage(self):
        """Test help message"""
        help_result = self.runner.invoke(querynator_cli, ["--help"])
        self.assertEqual(help_result.exit_code, 0)
        self.assertIn("--help     Show this message and exit.", help_result.output)
        self.assertNotIn("+++lalelu+++this is not part of the help message+++", help_result.output)

    def test_CgiHelp(self):
        """Test CGI help message"""
        cgi_help_result = self.runner.invoke(querynator_cli, ["query-api-cgi", "--help"])
        self.assertEqual(cgi_help_result.exit_code, 0)
        self.assertIn("--help                          Show this message and exit.", cgi_help_result.output)

    def test_CivicHelp(self):
        """Test CIVIC help message"""
        civic_help_result = self.runner.invoke(querynator_cli, ["query-api-civic", "--help"])
        self.assertEqual(civic_help_result.exit_code, 0)
        self.assertIn("--help                          Show this message and exit.", civic_help_result.output)

    def test_CreatereportHelp(self):
        """Test CREATE-REPORT help message"""
        report_help_result = self.runner.invoke(querynator_cli, ["create-report", "--help"])
        self.assertEqual(report_help_result.exit_code, 0)
        self.assertIn("--help                 Show this message and exit.", report_help_result.output)

    def test_nonExistingSubcommand(self):
        """Test non-existing subcommand"""
        result = self.runner.invoke(querynator_cli, ["query-api-clinvar"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("No such command", result.output)

    def test_nonExistingOption(self):
        """Test non-existing option"""
        result = self.runner.invoke(querynator_cli, ["query-api-cgi", "--baz"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("No such option", result.output)


class testCliRun(unittest.TestCase):
    unittest_dir = f"{os.getcwd()}/example_files/unittest_out"

    # helper functions
    def assertIsFile(self, path, msg=None):
        msg = f"{path} is not a file" if msg is None else msg
        self.assertTrue(os.path.isfile(path), msg)

    def get_testdir(self):
        """Get the test directory for the current test method."""
        method_name = str(self).split()[0]
        return f"{self.unittest_dir}/{method_name}"

    def clean_testdir(self):
        try:
            shutil.rmtree(self.get_testdir())
        except FileNotFoundError:
            pass

    # tests
    def setUp(self):
        self.runner = CliRunner()
        self.clean_testdir()

    def test_queryApiCivic(self):
        """Test query-api-civic"""
        outdir = self.get_testdir()
        result = self.runner.invoke(
            querynator_cli,
            [
                "query-api-civic",
                "--vcf",
                f"{os.getcwd()}/example_files/example.vcf",
                "--genome",
                "GRCh37",
                "--outdir",
                outdir,
                "--filter_vep",
            ],
        )
        self.assertEqual(result.exit_code, 0, "non-zero exit code")
        self.assertIsFile(f"{outdir}/{outdir.split('/')[-1]}.civic_results.tsv", "civic_results.tsv not created")

    def test_queryApiCivic_Cancer(self):
        """Test query-api-civic"""
        outdir = self.get_testdir()

        result = self.runner.invoke(
            querynator_cli,
            [
                "query-api-civic",
                "--vcf",
                "example_files/example.vcf",
                "--genome",
                "GRCh37",
                "--outdir",
                outdir,
                "--cancer",
                "Cholangiocarcinoma",
                "--filter_vep",
            ],
        )
        self.assertEqual(result.exit_code, 0, "non-zero exit code")
        self.assertIsFile(f"{outdir}/{outdir.split('/')[-1]}.civic_results.tsv", "civic_results.tsv not created")

    def test_queryApiCgi(self):
        """test querynator query-api-cgi with dummy credentials"""
        outdir = self.get_testdir()
        result = self.runner.invoke(
            querynator_cli,
            [
                "query-api-cgi",
                "--mutations",
                f"{os.getcwd()}/example_files/example.vcf",
                "--genome",
                "GRCh37",
                "--outdir",
                outdir,
                "--email",
                "dummy@internet.by",
                "--token",
                "thisisadummytoken",
                "--filter_vep",
            ],
        )
        self.assertEqual(result.exit_code, 1)

    def test_queryApiCgi_Cancer(self):
        """test querynator query-api-cgi with dummy credentials"""
        outdir = self.get_testdir()
        result = self.runner.invoke(
            querynator_cli,
            [
                "query-api-cgi",
                "--mutations",
                f"{os.getcwd()}/example_files/example.vcf",
                "--genome",
                "GRCh37",
                "--outdir",
                outdir,
                "--cancer",
                "Breast adenocarcinoma",
                "--email",
                "dummy@internet.by",
                "--token",
                "thisisadummytoken",
                "--filter_vep",
            ],
        )
        self.assertEqual(result.exit_code, 1)

    def test_queryApiCgi_invalidCancer(self):
        """test querynator query-api-cgi with dummy credentials"""
        outdir = self.get_testdir()
        result = self.runner.invoke(
            querynator_cli,
            [
                "query-api-cgi",
                "--mutations",
                f"{os.getcwd()}/example_files/example.vcf",
                "--genome",
                "GRCh37",
                "--outdir",
                outdir,
                "--cancer",
                "this is not cancer, 'tis but a flesh wound",
                "--email",
                "dummy@internet.by",
                "--token",
                "thisisadummytoken",
                "--filter_vep",
            ],
        )
        self.assertEqual(result.exit_code, 2)

    def test_createReport(self):
        """test querynator create-report on example files"""
        outdir = self.get_testdir()
        result = self.runner.invoke(
            querynator_cli,
            [
                "create-report",
                "--cgi_path",
                f"{os.getcwd()}/example_files/cgi_test_out",
                "--civic_path",
                f"{os.getcwd()}/example_files/civic_test_out",
                "--outdir",
                outdir,
            ],
        )
        self.assertEqual(result.exit_code, 0, "non-zero exit code")
        self.assertIsFile(f"{outdir}/combined_files/civic_cgi_vep.tsv", "combined file not created")
        self.assertIsFile(f"{outdir}/report/{outdir.split('/')[-1]}_overall_report.html", "report not created")

    def tearDown(self):
        """Clean up output directories after tests are run."""
        self.clean_testdir()


if __name__ == "__main__":
    unittest.main()
