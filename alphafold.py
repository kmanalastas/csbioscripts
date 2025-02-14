# This file is part of the csbioscripts package
#
# Copyright (c) 2025 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import json
import os
from Bio import SeqIO

def af3localrunjson(runname, sequencelist, seed=[1], version=1, bondedpairlist=None, userccdstr=None):
    runjson = {}
    outname = f'{runname}.json'
    
    runjson['name'] = runname
    runjson['modelSeeds'] = seed
    runjson['sequences'] = sequencelist
    if bondedpairlist != None:
        runjson['bondedAtomPairs'] = bondedpairlist
    if userccdstr != None:
        runjson['userCCD'] = userccdstr
    runjson['dialect'] = 'alphafold3'
    runjson['version'] = version
    with open(outname, 'w') as f:
        json.dump(runjson, f)
        print (f'output file: {outname}')

def buildproteinsequence(chainid, sequence, modifications=None, unpairedmsa=None, unpairedmsapath=None, pairedmsa=None, pairedmsapath=None, templates=None):
    params = {}
    params['protein'] = {}
    params['protein']['id'] = chainid
    params['protein']['sequence'] = sequence
    if modifications != None:
        params['protein']['modifications'] = modifications
        
    # unpairedMsa and unpairedMsaPath are mutually exclusive
    if unpairedmsa != None:
        params['protein']['unpairedMsa'] = unpairedmsa
    elif unpairedmsapath != None:
        params['protein']['unpairedMsaPath'] = unpairedmsapath

    # pairedMsa and pairedMsaPath are mutually exclusive
    if pairedmsa != None:
        params['protein']['pairedMsa'] = pairedmsa
    elif pairedmsapath != None:
        params['protein']['pairedMsaPath'] = pairedmsapath
    
    if templates != None:
        params['protein']['templates'] = templates

    return params
    
def af3jobfromfasta(fastafile):
    basename = os.path.basename(fastafile)
    jobname, ext = os.path.splitext(basename)            
    records = list(SeqIO.parse(fastafile, 'fasta'))
    sequences = []
    chainid = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    for i,rec in enumerate(records):
        sequences.append(buildproteinsequence(chainid[i], str(rec._seq)))
    af3localrunjson(jobname, sequences)        
    
