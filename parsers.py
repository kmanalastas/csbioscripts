# This file is part of the csbioscripts package
#
# Copyright (c) 2023 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import json

def findfieldinjson(infile, field):
    with open(infile, "r") as f:
        buf = json.load(f)
        ans = recursivesearch(buf, field)
    return ans

def recursivesearch(level, field):
    if field in level:
        return level[field]
    else:
        for i in level:
            return recursivesearch(i, field)
