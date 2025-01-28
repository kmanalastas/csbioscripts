# This file is part of the csbioscripts package
#
# Copyright (c) 2024 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import os
import json
from fetchfromdb import downloadpage

class Uniprot:
    def __init__(self, upid):
        self.id = upid
        self.pdbentries = None
        self.nrchains = None
        
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
        
        
        
        
            
                
        
        
