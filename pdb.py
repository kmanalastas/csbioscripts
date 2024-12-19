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

class PDBentry:
    def __init__(self, pdbid):
        self.id = pdbid
        self.resolution = None
        self.expmethod = None
        self.filepath = None
        
    def fetchpdbidresolution(self):
        jsonfile = downloadpage('https://data.rcsb.org/rest/v1/core/entry', self.id)
        if os.path.exists(jsonfile):        
            with open(jsonfile, "r") as f:
                buf = json.load(f)
                if 'rcsb_entry_info' in buf:
                    if 'resolution_combined' in buf['rcsb_entry_info']: 
                        cres = (buf['rcsb_entry_info']['resolution_combined'][0])
                        expmethod = buf['rcsb_entry_info']['experimental_method']
                        self.resolution, self.expmethod = cres, expmethod
    
    

