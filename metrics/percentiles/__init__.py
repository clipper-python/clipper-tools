"""
Copyright 2020 William Rochira at the University of York
Developed at York Structural Biology Laboratory
"""

import os

from .._defs import METRIC_NAMES, METRIC_POLARITIES, RESOLUTION_BIN_NAMES, RESOLUTION_PERCENTILES


INPUT_PATH = './data/percentiles_data.csv'
PERCENTILE_DATA = None


def _load_data():
    global PERCENTILE_DATA
    PERCENTILE_DATA = { }
    with open(os.path.join(os.path.dirname(__file__), INPUT_PATH), 'r') as infile:
        for i, line in enumerate(infile.readlines()):
            splitline = line.strip().split(',')
            if i == 0:
                metric_names = splitline[2:]
                for metric_name in metric_names:
                    PERCENTILE_DATA[metric_name] = { }
                    for bin_name in RESOLUTION_BIN_NAMES:
                        PERCENTILE_DATA[metric_name][bin_name] = { }
            else:
                bin_name = splitline[0]
                percentile = int(splitline[1])
                metric_values = [ float(x) for x in splitline[2:] ]
                for metric_name, metric_value in zip(metric_names, metric_values):
                    PERCENTILE_DATA[metric_name][bin_name][percentile] = metric_value
    for metric_name in METRIC_NAMES:
        if metric_name not in PERCENTILE_DATA.keys():
            raise Exception('no percentile data for metric: ' + metric_name)


def _get_bin_name(resolution):
    if resolution is None:
        return 'All'
    bin_id = 9
    for i, percentile in enumerate(sorted(RESOLUTION_PERCENTILES.keys())):
        percentile_resolution = RESOLUTION_PERCENTILES[percentile]
        if resolution < percentile_resolution:
            bin_id = i
            break
    bin_name = RESOLUTION_BIN_NAMES[bin_id]
    return bin_name


def get_percentile(metric_id, metric_value, resolution=None, normalise_polarity=False):
    if PERCENTILE_DATA is None:
        _load_data()

    if None in (metric_id, metric_value) or metric_id > len(METRIC_NAMES)-1:
        return None

    metric_name = METRIC_NAMES[metric_id]
    bin_name = _get_bin_name(resolution)
    determined_percentile = 100
    for percentile, percentile_value in PERCENTILE_DATA[metric_name][bin_name].items():
        if metric_value < percentile_value:
            determined_percentile = percentile
            break
    if normalise_polarity and METRIC_POLARITIES[metric_id] == -1:
        return 101 - determined_percentile
    else:
        return determined_percentile
