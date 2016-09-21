from amber_adaptbx.tests.utils import get_fn

saved_2igd_rst7_file = get_fn('2igdab.LES.rst7')
saved_2igd_prmtop_file = get_fn('2igdab.LES.prmtop')
saved_2igd_pdb_file = get_fn('2igdab_4phenix.LES.pdb')
saved_2igd_mtz_file = get_fn('2igd.mtz')

PDB_COLLECTION = [
    get_fn('2igd.pdb'),
    get_fn('4lzt/4lzt_no_BHOH.pdb'),
]

# must have the same order as PDB_COLLECTION
MTZ_COLLECTION = [
    get_fn('2igd.mtz'),
    get_fn('4lzt/4lzt.mtz'),
]

PDB_MTZ_COLLECTION = list(zip(PDB_COLLECTION, MTZ_COLLECTION))[::-1]
