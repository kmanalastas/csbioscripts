# This file is part of the csbioscripts package
#
# Copyright (c) 2024 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import urllib.request
import sys
import json
import os

def fetchpdbidresolution(pdbid):
    jsonfile = downloadpdbmetadata(pdbid)
    if os.path.exists(jsonfile):        
        with open(jsonfile, "r") as f:
            buf = json.load(f)
            if 'rcsb_entry_info' in buf:
                if 'resolution_combined' in buf['rcsb_entry_info']: 
                    cres = (buf['rcsb_entry_info']['resolution_combined'][0])
                    expmethod = buf['rcsb_entry_info']['experimental_method']
                    return cres, expmethod
            else:
                return None
    else:
        return None

def downloadpdbmetadata(pdbid, name=None):
    if name == None:
        name = f'{pdbid}.json'
    if not os.path.exists(name):
        urllib.request.urlretrieve(f'https://data.rcsb.org/rest/v1/core/entry/{pdbid}', name)
    return name


def downloadpage(baseurl, suffix, filename=None):
    if name == None:
        name = f'{suffix}.json'
    if not os.path.exists(name):
        urllib.request.urlretrieve(f'{baseurl}/{suffix}', name)
    return name
    


        
    