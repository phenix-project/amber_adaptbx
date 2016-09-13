#!/usr/bin/env phenix.python
import sys
import pytest

args = sys.argv[1:] if sys.argv[1:]  else ['-vs', '.']

pytest.main(args)

__doc__ = """

This script is a thin wrapper of pytest (phenix.python -m pytest)

pytest document: http://doc.pytest.org/en/latest/

Usage
=====

(for those who are not familiar with pytest)

- Run full tests with verbose

    ./run_tests.py -vs

- Run a given test

    # specific test name
    ./run_tests.py -vs les/test_builder.py::test_run_sander_LES_min

    # specific folder
    ./run_tests.py -vs normal/

- Run a group of tests that you tag with a given name

    ./run_tests.py -vs -m slowtest
    # we tag the function with @pytest.makr.slowtest

- Many more

    ./run_tests.py --help

- I also made Makefile to run the most commont tasks.

Why using pytest?
=================

Many reasons but I can list some of them

- auto discover test files: You just need to make a filename starting with "test_"
and write your test (e.g. def test_abc)

    def test_abc():
        assert 1 + 1 == 2

- flexible in running a group of tests
"""
