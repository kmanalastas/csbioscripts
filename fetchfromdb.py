# This file is part of the csbioscripts package
#
# Copyright (c) 2023 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import urllib.request
import sys

def uniprottopdb(uniprotid):
    jsonfile = downloadsiftsmapping(uniprotid)	

def downloadsiftsmapping(uniprotid, name=None):
    if name == None:
        name = f'{uniprotid}.json'
    urllib.request.urlretrieve(f'https://www.ebi.ac.uk/pdbe/api/mappings/{uniprotid}', name)
    return name

if __name__ == '__main__':
    upid = sys.argv[1]
    uniprottopdb(upid)
    