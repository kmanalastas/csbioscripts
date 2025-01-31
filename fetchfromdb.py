# This file is part of the csbioscripts package
#
# Copyright (c) 2024 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import urllib.request
import os


def downloadpage(baseurl, suffix, filename=None):
    if filename == None:
        filename = f'{suffix}.json'
    if not os.path.exists(filename):
        try:
            response = urllib.request.urlopen(f'{baseurl}/{suffix}')
            webContent = response.read().decode('UTF-8')
            
            with open(filename, 'w') as f:
                f.write(webContent)
        except:
            print (f'No entries for {suffix}')
        
    return filename
    


        
    