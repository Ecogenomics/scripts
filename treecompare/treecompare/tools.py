import os
from collections import defaultdict


def maybe_make_dirs(path):
    dir_name = os.path.dirname(path)
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)


def convert_tax_to_ranks(tax_str):
    out = defaultdict(set)
    for rank in tax_str.split(';'):
        out[rank[0]].add(rank)
    return out
