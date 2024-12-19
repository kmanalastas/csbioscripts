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

def uniprottopdb(uniprotid):
    results = []
    jsonfile = downloadsiftsmapping(uniprotid)
    with open(jsonfile, 'r') as f:
        buf = json.load(f)
    pdbentries = buf[uniprotid]['PDB']
    for pdbid in pdbentries:
        for instance in buf[uniprotid]['PDB'][pdbid]:
            match = {
                'pdbid': pdbid,
                'chain_id': instance['chain_id'],
                'unp_start': instance['unp_start'],
                'unp_end': instance['unp_end'],
                'pdb_start': instance['start']['residue_number'],
                'pdb_end': instance['end']['residue_number'],
            }
            results.append(match)
    return results

def downloadsiftsmapping(uniprotid, name=None):
    if name == None:
        name = f'{uniprotid}.json'
    urllib.request.urlretrieve(f'https://www.ebi.ac.uk/pdbe/api/mappings/{uniprotid}', name)
    return name

if __name__ == '__main__':
    upid = sys.argv[1]
    pdbrecs = uniprottopdb(upid)
    for i in pdbrecs:
        print (i)
    