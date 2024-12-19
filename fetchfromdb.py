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
    if not os.path.exists(name):
        urllib.request.urlretrieve(f'https://www.ebi.ac.uk/pdbe/api/mappings/{uniprotid}', name)
    return name

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


if __name__ == '__main__':
    upid = sys.argv[1]
    pdbrecs = uniprottopdb(upid)
    for i in pdbrecs:
        print (i)
        reso = fetchpdbidresolution(i['pdbid'])
        print ('resolution:', reso)
        
    