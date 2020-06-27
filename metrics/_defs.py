"""
Copyright 2020 William Rochira at the University of York
"""

SC_INCOMPLETE_STRING = 'INCOMPLETE SIDECHAIN'
NO_CONTEXT_MESSAGE = 'NOCONTEXT'
METRIC_NAMES = ('Ramachandran Score', 'Rotamer Score', 'Avg B-factor', 'Max B-factor', 'Mainchain Fit', 'Sidechain Fit')
METRIC_POLARITIES = (+1, -1, -1, -1, -1, -1)
RESOLUTION_BIN_NAMES = ('<10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '>90', 'All')
RESOLUTION_PERCENTILES = { 10 : 1.48,
                           20 : 1.67,
                           30 : 1.80,
                           40 : 1.92,
                           50 : 2.04,
                           60 : 2.20,
                           70 : 2.40,
                           80 : 2.60,
                           90 : 2.90 }
