#!/usr/bin/env python

"""Tests for `querynator` package."""

from unittest import mock

from click.testing import CliRunner

import querynator.__main__


@mock.patch("querynator.__main__.querynator_cli")
def test_header(mock_cli):
    """Just try to execute the header function"""
    querynator.__main__.run_querynator()


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(querynator.__main__.querynator_cli)
    assert result.exit_code == 0

    """Test help message"""
    help_result = runner.invoke(querynator.__main__.querynator_cli, ["--help"])
    assert help_result.exit_code == 0
    assert "--help     Show this message and exit." in help_result.output

    """Test non-existing subcommand"""
    result = runner.invoke(querynator.__main__.querynator_cli, ["query-api-clinvar"])
    assert result.exit_code == 2
    assert "No such command" in result.output

    """Test non-existing option"""
    result = runner.invoke(querynator.__main__.querynator_cli, ["query-api-cgi", "--baz"])
    assert result.exit_code == 2
    assert "No such option" in result.output
