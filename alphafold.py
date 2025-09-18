# This file is part of the csbioscripts package
#
# Copyright (c) 2025 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import json
import os
from Bio import SeqIO, PDB
from csbioscripts.crosslinks import parsecrosslinks_xlmstools_af3
from csbioscripts.fetchfromdb import downloadpage
from csbioscripts.parsers import findfieldinjson
import matplotlib.pyplot as plt
import numpy as np
import json

class AFdbmodel:
    def __init__(self, uniprotid, directory=None):
        self.uniprotid = uniprotid
        self.biopystructure = None
        self.pae = None
        self.fetchmodel(directory=directory)
    
    def fetchmodel(self, directory=None):
        # download PDB model
        pdbfile = f'{self.uniprotid}.pdb'
        pdbpath = downloadpage('https://alphafold.ebi.ac.uk/files', f'AF-{self.uniprotid}-F1-model_v4.pdb', directory=directory, filename=pdbfile)
        print (pdbpath)
        if pdbpath != None:
            parser = PDB.PDBParser()
            self.biopystructure = parser.get_structure(pdbpath, pdbpath)
            # download PAE plot
            paefile = f'{self.uniprotid}_pae.json'
            paepath = downloadpage('https://alphafold.ebi.ac.uk/files', f'AF-{self.uniprotid}-F1-predicted_aligned_error_v4.json', directory=directory, filename=paefile)
            self.pae = findfieldinjson(paefile, 'predicted_aligned_error')
        
def af3localrunjson(runname, sequencelist, seed=[1], version=1, bondedpairlist=None, userccdstr=None, crosslinks=None):
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
    if crosslinks != None:
        runjson['crosslinks'] = crosslinks
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
    
def af3jobfromfasta(fastafile, crosslinkfile=None, crosslinker=None, nseeds=None):
    basename = os.path.basename(fastafile)
    jobname, ext = os.path.splitext(basename)            
    records = list(SeqIO.parse(fastafile, 'fasta'))
    sequences = []
    chainid = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    for i,rec in enumerate(records):
        sequences.append(buildproteinsequence(chainid[i], str(rec._seq)))
    
    crosslinks = None
    if crosslinkfile != None:
        crosslinks = parsecrosslinks_xlmstools_af3(crosslinkfile, crosslinker)
    if nseeds != None:
        seeds = [i for i in range(1, nseeds+1)]
    else:
        seeds = [1]
    af3localrunjson(jobname, sequences, crosslinks=crosslinks, seed=seeds)
    
def plotpae(jsonfile):
    outfig = os.path.splitext(jsonfile)[0] + '_pae.pdf'
    pae = findfieldinjson(jsonfile, 'pae')
    seqlen = len(pae[0])
    if pae != None:
        token_chain_ids = findfieldinjson(jsonfile, 'token_chain_ids')
        breaks = [(i+i-1)/2 for i in range(1,seqlen) if token_chain_ids[i] != token_chain_ids[i-1]]
        plt.imshow(pae, cmap='Greens_r')
        for i in breaks:
            plt.plot([0, seqlen], [i,i], 'k--', alpha=0.5)
            plt.plot([i,i], [0, seqlen], 'k--', alpha=0.5)
        plt.colorbar(label='Expected position error (Ångströms)')
        plt.xlabel('Scored residue')
        plt.ylabel('Aligned residue')        
        plt.savefig(outfig)
        print (f'Saved PAE plot: {outfig}')

def highlightlowPAEinterchain(jsonfile, cxsession, maxpae=5):
    pae = findfieldinjson(jsonfile, 'pae')
    seqlen = len(pae[0])
    if pae == None:
        print (f'PAE values not found in {jsonfile}')
        return 0
    
    # parse low PAE values between chains
    mark = {}
    chain = findfieldinjson(jsonfile, 'token_chain_ids')
    resi = findfieldinjson(jsonfile, 'token_res_ids') 
    contactmap = np.where(np.array(pae) <= maxpae, 1, 0)
    for i in range(0, seqlen):
        for j in range(0, seqlen):
            if chain[i] != chain[j]:
                if contactmap[i,j] == 1:
                    mark = add_i(i, mark, chain, resi)
                    mark = add_i(j, mark, chain, resi)
    print (mark)
    
    # generate image with interchain low PAE colored in red
    outcxs = os.path.splitext(cxsession)[0] + f'_highlight_{str(maxpae)}.cxc'
    outfig = os.path.splitext(cxsession)[0] + f'_highlight_{str(maxpae)}.png'
    with open(outcxs, 'w') as f:
        f.write(f'open {cxsession}\n')
        for chid in mark.keys():
            f.write(f'color #1/{chid}:')
            f.write(f'{mark[chid][0]}')
            for rid in mark[chid][1:]:
                f.write(f', {rid}')
            f.write(' red cartoon\n')
        f.write(f'save {outfig}')
         
                    
def add_i(i, mark, chain, resi):
    if chain[i] in mark:
        if resi[i] not in mark[chain[i]]:
            mark[chain[i]].append(resi[i])
    else:
        mark[chain[i]] = [resi[i]]
    return mark
    
    
    
