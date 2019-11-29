"""
Copyright 2019 William Rochira at the University of York
Developed at York Structural Biology Laboratory - Cowtan group

 - Loads serialised rotamer library to calculate probability of
   a particular sidechain conformation from a (semi) continuous
   probability distribution
 - Included rotamer library data originally from the Richardson
   lab (https://github.com/rlabduke/reference_data)
"""
import os
import pickle

from .. import _defs


ROTAMER_LIBRARY = None


def load_rotamer_lib():
    global ROTAMER_LIBRARY
    with open(os.path.join(os.path.dirname(__file__), 'lib', 'rotamer_library.pkl' ), 'rb') as infile:
        ROTAMER_LIBRARY = pickle.load(infile)


def get_probability(code, chis):
    if ROTAMER_LIBRARY is None:
        load_rotamer_lib()
    if code not in ROTAMER_LIBRARY.keys() or _defs.SC_INCOMPLETE_STRING in chis:
        return None
