"""
Copyright 2020 William Rochira at the University of York
Developed at York Structural Biology Laboratory

Included rotamer data originally from the Richardson lab
(https://github.com/rlabduke/reference_data)
"""

import os
import pickle

import numpy as np

from .. import _defs


ROTAMER_LIBRARY_DATA = None
ROTAMER_CENTRAL_VALUES = None


def _load_library():
    global ROTAMER_LIBRARY_DATA
    try:
        import gzip
    except ImportError:
        print('WARNING: failed to import GZIP. Rotamer library cannot be loaded')
        return
    with gzip.open(os.path.join(os.path.dirname(__file__), 'data', 'library.gz'), 'rb') as infile:
        dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, compressed_byte_arrays = pickle.load(infile)
    classification_bytes = { }
    for code, compressed in compressed_byte_arrays.items():
        compressed = bytearray(compressed) # Need to cast from str in Python 2
        decompressed = bytearray()
        for combo_byte in compressed:
            b0 = combo_byte // 16
            b1 = combo_byte - 16 * b0
            decompressed.append(b0)
            decompressed.append(b1)
        classification_bytes[code] = decompressed
    ROTAMER_LIBRARY_DATA = (dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, classification_bytes)


def _load_central_values():
    global ROTAMER_CENTRAL_VALUES
    central_values = { }
    with open(os.path.join(os.path.dirname(__file__), 'data', 'central_values.csv'), 'r') as infile:
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
        sqdiffs = [ ]
        for i in range(len(chis)):
            deltas = (chis[i]-chi_means[i], chis[i]-chi_means[i]+360)
            best_delta = deltas[int(deltas[0]>deltas[1])]
            sqdiff = (best_delta/chi_sdevs[i])**2
            sqdiffs.append(sqdiff)
        score = (sum(sqdiffs)/len(sqdiffs))**0.5
        rotamer_scores[rot_name] = score
    return rotamer_scores


def get_classification(code, chis):
    if  _defs.SC_INCOMPLETE_STRING in chis:
        return None
    if ROTAMER_LIBRARY_DATA is None:
        _load_library()
    if ROTAMER_LIBRARY_DATA is None:
        return None
    dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, classification_bytes = ROTAMER_LIBRARY_DATA
    if code not in dim_offsets.keys():
        return None
    closest_values = [ ]
    chis = tuple([ x for x in chis if x is not None ][:len(dim_offsets[code])])
    for dimension, chi in enumerate(chis):
        dim_width = dim_bin_ranges[code][dimension][1] - dim_bin_ranges[code][dimension][0]
        if chi <= dim_bin_ranges[code][dimension][0]:
            chi += dim_width
        if chi >= dim_bin_ranges[code][dimension][1]:
            chi -= dim_width
        multiple = round((chi - dim_offsets[code][dimension]) / float(dim_bin_widths[code][dimension]))
        closest_value = dim_offsets[code][dimension] + multiple * dim_bin_widths[code][dimension]
        closest_values.append(closest_value)
    closest_values = tuple(closest_values)
    index = 0
    for dimension, chi in enumerate(closest_values):
        dim_offest = dim_offsets[code][dimension]
        dim_bin_width = dim_bin_widths[code][dimension]
        index += int((chi - dim_offest) / dim_bin_width * np.product(dim_num_options[code][dimension+1:]))
    return classification_bytes[code][index]


def get_cv_score(code, chis):
    if ROTAMER_CENTRAL_VALUES is None:
        _load_central_values()
    if code not in ROTAMER_CENTRAL_VALUES.keys() or _defs.SC_INCOMPLETE_STRING in chis:
        return None
    scores = _cv_sqdiff_scores(code, chis)
    best_score = min(scores.values())
    return best_score
