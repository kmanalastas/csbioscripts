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
        
    def sortedPDBstructures(self):
        jsonfile = downloadpage('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures', self.id)
        if os.path.exists(jsonfile):
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if self.id in buf:
                    pdbentries = buf[self.id]
                    self.pdbentries = pdbentries
                else:
                    self.pdbentries = []
    
    def representativepdbchains(self):
        if self.nrchains != None:
            return self.nrchains
        else:
            if self.pdbentries == None:
                self.sortedPDBstructures()
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

    def getalternateconformations(self):
        if self.conformations != None:
            return self.conformations
        else:
            if self.nrchains == None:
                self.representativepdbchains()
            self.conformations = []
            for i in self.nrchains:
                allpdbreps = [j for j in self.pdbentries if self.checkoverlap(i,j)]
                chains = [] # paths to pdb files
                for j in allpdbreps:
                    #print ('j', j)
                    chains += self.getchains(j)
                clustered = lazycluster(chains, tmscorewrapper, operator.ge, 0.5)
                self.conformations.append(clustered)
            return self.conformations
    
    def getchains(self, pdbrec):
        chains = []
        pdbent = PDBentry(pdbrec['pdb_id'])
        pdbent.fetchbiopythonstructure()
        chains += pdbent.printchainaspdb(pdbrec['chain_id'], separate=True)
        return chains
            
        
        
        
        
            
                
        
        
