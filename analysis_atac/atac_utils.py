"""
"""
import numpy as np
import pandas as pd

SAMPLE_TO_CONDCODE = {
    'P6a': 0,
    'P6b': 0,
    'P6c': 0,
    
    'P8a': 1,
    'P8b': 1,
    'P8c': 1,
    
    'P10a': 2,
    'P10b': 2,
    
    'P12a': 3,
    'P12b': 3,
    'P12c': 3,
        
    'P12DRa': 4,
    'P12DRb': 4,
    
    'P14a': 5,
    'P14b': 5,
        
    'P14DRa': 6,
    'P14DRb': 6,
    
    'P17a': 7,
    'P17b': 7,
        
    'P17DRa': 8,
    'P17DRb': 8,
        
    'P21a': 9,
    'P21b': 9,
        
    'P21DRa': 10,
    'P21DRb': 10,
}

CONDCODE_TO_COND = {
    0: 'P6',
    1: 'P8',
    2: 'P10',
    3: 'P12',
    4: 'P12DR',
    5: 'P14',
    6: 'P14DR',
    7: 'P17',
    8: 'P17DR',
    9: 'P21',
   10: 'P21DR',
}

# CONDCODE_TO_COND_v2 = {
#     0: 'P6NR',
#     1: 'P8NR',
#     2: 'P10NR',
#     3: 'P12NR',
#     4: 'P12DR',
#     5: 'P14NR',
#     6: 'P14DR',
#     7: 'P17NR',
#     8: 'P17DR',
#     9: 'P21NR',
#    10: 'P21DR',
# }

def merge_peaks(df, chrom, start, end):
    """df is a dataframe
    chrom, start, end are column numbers
    """
    return df[chrom].astype(str)+":"+df[start].astype(str)+"-"+df[end].astype(str)

def get_integer_bins(vec):
    """excluding 0; 
    bins = [0, 1, ..., max(vec)] + 0.5
    mids = [   1, ..., max(vec)]
    """
    bins = np.arange(0, np.max(vec)+1) + 0.5
    mids = np.arange(1, np.max(vec)+1)
    frq_sig_a, _ = np.histogram(vec, bins=bins)
    
    return frq_sig_a, mids 


def encode_sample(sample):
    """
    """
    return SAMPLE_TO_CONDCODE[sample]

def decode_cond(code):
    """
    """
    return CONDCODE_TO_COND[code]


def get_reverse_map(forward_map):
    """dictionary key -> item
    
    reverse 
        item (unique ones) -> key(s)
    """
    uniq_items = np.unique(pd.Series(forward_map).values) 
    
    rmap = {item: [] for item in uniq_items}

    for key, item in forward_map.items():
        rmap[item].append(key)
        
    return rmap

def get_organized_sample_dict():
    return get_reverse_map(SAMPLE_TO_CONDCODE)
