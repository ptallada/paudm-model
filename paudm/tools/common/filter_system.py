#!/usr/bin/env python
# encoding: UTF8

# Interface for loading the filter throughput into a Pandas
# dictionary.
import os
import glob
import pandas as pd
from pkg_resources import resource_filename

def _read_csv(path):
    """Read in the SED files. These are all stored in the same format."""

    band = os.path.basename(path).replace('.sed', '')
    data = pd.io.api.read_csv(path, sep=' ', comment='#', names=['lmb', 'resp'])
    data = data.set_index('lmb')

    return band, data

def filter_throughput():
    """The full filter throughput (filter, atmosphere, detector)."""

    d = resource_filename('paudm.resources', 'throughput')
    sed_glob = os.path.join(d, '*', '*.sed')
    D = dict(map(_read_csv, glob.glob(sed_glob)))
    resp = pd.concat(D)

    return resp
