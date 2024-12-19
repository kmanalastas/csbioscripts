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
        
    def sortedPDBstructures(self):
        jsonfile = downloadpage('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures', self.id)
        if os.path.exists(jsonfile):
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                pdbentries = buf[self.id]
                self.pdbentries = pdbentries
                
        
        
