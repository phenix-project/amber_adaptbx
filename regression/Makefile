# run some most common tests.
# Please see the doc in run_tests.py too

all:
	./run_tests.py -vs

# run with code coverage, require pytest-cov
all_coverage:
	./run_tests.py -vs --cov=amber_adaptbx --cov-report=html

# tests that are slow (e.g energy minimization with maxcyc=1000)
slow:
	./run_tests.py -vs -m slowtest

fast:
	./run_tests.py -vs -m "not slowtest"

les::
	./run_tests.py -vs ./les/

normal::
	./run_tests.py -vs ./normal/

# run tag: 
# @pytest.mark.minimization
min:
	./run_tests.py -vs -m minimization

amber_phenix_converter:
	./run_tests.py -vs ./les/test_amber_phenix_reorder.py
