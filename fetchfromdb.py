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
        urllib.request.urlretrieve(f'{baseurl}/{suffix}', filename)
    return filename
    


        
    