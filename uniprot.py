# This file is part of the csbioscripts package
#
# Copyright (c) 2024 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import os
import json
from csbioscripts.fetchfromdb import downloadpage
from csbioscripts.pdb import PDBentry, tmscorewrapper
from csbioscripts.misc import lazycluster
import operator

class Uniprot:
    def __init__(self, upid):
        self.id = upid
        self.pdbentries = None
        self.nrchains = None
        self.conformations = None
        self.sequence = None
        
    def sortedPDBstructures(self, directory=None):
        jsonfile = downloadpage('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures', self.id, directory=directory)
        if os.path.exists(jsonfile):
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if self.id in buf:
                    pdbentries = buf[self.id]
                    self.pdbentries = pdbentries
                else:
                    self.pdbentries = []
    
    def representativepdbchains(self, directory=None):
        if self.nrchains != None:
            return self.nrchains
        else:
            if self.pdbentries == None:
                self.sortedPDBstructures(directory=directory)
            self.nrchains = []
            if len(self.pdbentries) > 0:
                for i in self.pdbentries:
                    discardi = False
                    for j in self.nrchains:
                        if self.checkoverlap(j, i):
                            discardi = True
                            break
                    if discardi == False:
                        self.nrchains.append(i)
            return self.nrchains
    
    def checkoverlap(self, pdba, pdbb, numadded=50):
        overlap = False
        s1 = pdba['unp_start']
        e1 = pdba['unp_end']
        s2 = pdbb['unp_start']
        e2 = pdbb['unp_end']
        if s2 >= s1 and s2 <= e1:
            extralen = e2 - e1
            if extralen <= numadded:
                overlap = True
            else:
                if extralen < 0.5*(e2-s2):
                    overlap = True
        elif e2 >= s1 and e2 <= e1:
            extralen = s1 - s2
            if extralen <= numadded:
                overlap = True
            else:
                if extralen < 0.5*(e2-s2):
                    overlap = True
        return overlap

    def getalternateconformations(self, removealtloc=False, removehetatm=False, keepligands=[], directory=None):
        if self.conformations != None:
            return self.conformations
        else:
            if self.nrchains == None:
                self.representativepdbchains(directory=directory)
            self.conformations = []
            for i in self.nrchains:
                allpdbreps = [j for j in self.pdbentries if self.checkoverlap(i,j)]
                chains = [] # paths to pdb files
                for j in allpdbreps:
                    #print ('j', j)
                    add = self.getchains(j, removealtloc=removealtloc, removehetatm=removehetatm, keepligands=keepligands, directory=directory)
                    if len(add) > 0:
                        chains += add
                if len(chains) > 0: 
                    clustered = lazycluster(chains, tmscorewrapper, operator.ge, 0.9)
                    self.conformations.append(clustered)
            return self.conformations
    
    def getchains(self, pdbrec, removealtloc=False, removehetatm=False, keepligands=[], directory=None):
        chains = []
        pdbent = PDBentry(pdbrec['pdb_id'])
        pdbent.fetchbiopythonstructure(directory=directory)
        if pdbent.biopystruct != None:
            if removealtloc:
                pdbent.removealtloc
            if removehetatm:
                pdbent.biopystruct = pdbent.removeheteroatoms(exceptions=keepligands)
            chains += pdbent.printchainaspdb(pdbrec['chain_id'], separate=True)
        return chains
    
    def getsequence(self, directory=None):
        if self.sequence != None:
            return self.sequence
        else:
            seq = ''
            fastafile = f'{self.id}.fasta'
            fastafile = downloadpage('https://rest.uniprot.org/uniprotkb/', fastafile, directory=directory, filename=fastafile)
            with open(fastafile, 'r') as f:
                for line in f:
                    if line[0] != '>':
                        seq += line.strip()
            self.sequence = seq
            return self.sequence
            
                
        
        
        
        
            
                
        
        
