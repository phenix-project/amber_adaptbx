# run some most common tests.
# Please see the doc in run_tests.py too

test.all: test.long

test.long:
	# about > 80 minutes
	# :D
	./run_tests.py -v

test.all.verbose:
	./run_tests.py -vs

# run with code coverage, require pytest-cov
test.all.coverage:
	./run_tests.py -vs --cov=amber_adaptbx --cov-report=html

# tests that are slow (e.g energy minimization with maxcyc=1000)
test.slow:
	./run_tests.py -vs -m slow

test.fast.verbose:
	./run_tests.py -vs -m "not slow and not medium"

test.fast:
	# about 2-3 minutes
	./run_tests.py -v -m "not slow and not medium"

test.short:
	# about 2-3 minutes
	./run_tests.py -v -m "not slow and not medium"

test.medium:
	# about 6 - 10 minutes
	./run_tests.py -v -m "medium"

test.medium.verbose:
	./run_tests.py -vs -m "medium"

test.collect:
	# only collect test names
	# TODO: make nicer display
	./run_tests.py -vs --collect-only

test.geometry_minimization:
	./run_tests.py -vs ./test_geometry_minimization/

test.les:
	./run_tests.py -vs ./les/

test.refinement:
	./run_tests.py -vs -k 'refinement'

test.les.refinement:
	./run_tests.py -vs les/test_refinement.py

test.non_les.refinement:
	./run_tests.py -vs normal/test_refinement.py

test.min:
	./run_tests.py -vs -k minimization

test.saved:
	./run_tests.py -vs -m saved

test.amber_phenix_converter:
	./run_tests.py -vs ./les/test_amber_phenix_reorder.py

test.non_les:
	./run_tests.py -vs ./normal/

test.non_les_fast:
	./run_tests.py -vs ./normal/ -m "not slow"

test.non_les_refine_phenix:
	./run_tests.py -vs ./normal/ -m "phenix"

test.vAla3:
	./run_tests.py -v -k 'vAla3'
