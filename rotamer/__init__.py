"""
Copyright 2019 William Rochira at the University of York
Developed at York Structural Biology Laboratory - Cowtan group

Included rotamer library data originally from the Richardson
lab (https://github.com/rlabduke/reference_data)
"""

import os
import pickle

from .. import _defs


ROTAMER_LIBRARY = None
ROTAMER_CENTRAL_VALUES = None


def _load_library():
    global ROTAMER_LIBRARY
    with open(os.path.join(os.path.dirname(__file__), 'data', 'library.pkl' ), 'rb') as infile:
        ROTAMER_LIBRARY = pickle.load(infile)


def _load_central_values():
    global ROTAMER_CENTRAL_VALUES
    central_values = { }
    with open(os.path.join(os.path.dirname(__file__), 'data', 'central_values.csv' ), 'r') as infile:
        infile.readline() # Skip header line
        for line in infile.readlines():
            splitline = line.strip().split(',')
            code = splitline[0]
            rot_name = splitline[1]
            chi_means = [ float(x) for x in splitline[2:6] if x != 'None' ]
            chi_sdevs = [ float(x) for x in splitline[6:10] if x != 'None' ]
            if code not in central_values:
                central_values[code] = [ ]
            central_values[code].append((rot_name, chi_means, chi_sdevs))
    ROTAMER_CENTRAL_VALUES = central_values


def _cv_sqdiff_scores(code, chis):
    rotamer_scores = { }
    chis = [ x for x in chis if x is not None ]
    for rot_name, chi_means, chi_sdevs in ROTAMER_CENTRAL_VALUES[code]:
        chis = chis[:len(chi_means)]
        sqdiffs = [ ((chis[i]-chi_means[i])/chi_sdevs[i])**2 for i in range(len(chis)) ]
        score = (sum(sqdiffs)/len(sqdiffs))**0.5
        rotamer_scores[rot_name] = score
    return rotamer_scores


def get_probability(code, chis):
    if ROTAMER_LIBRARY is None:
        _load_library()
    if code not in ROTAMER_LIBRARY.keys() or _defs.SC_INCOMPLETE_STRING in chis:
        return None


def get_cv_score(code, chis):
    if ROTAMER_CENTRAL_VALUES is None:
        _load_central_values()
    if code not in ROTAMER_CENTRAL_VALUES.keys() or _defs.SC_INCOMPLETE_STRING in chis:
        return None
    scores = _cv_sqdiff_scores(code, chis)
    best_score = min(scores.values())
    return best_score
