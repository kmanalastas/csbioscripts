# This file is part of the csbioscripts package
#
# Copyright (c) 2025 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

def parsecrosslinks_xlmstools_af3(crosslinkfile, crosslinker):
    outlist = []
    with open(crosslinkfile, 'r') as f:
        entry = {}
        entry['name'] = crosslinker
        entry['residue_pairs'] = []
        for line in f:
            tmp = line.strip().split('|')
            if len(tmp) == 4:   # is a crosslink
                ch1 = tmp[1]
                res1 = int(tmp[0])
                ch2 = tmp[3]
                res2 = int(tmp[2])
                curxl = [[ch1, res1], [ch2, res2]]
                entry['residue_pairs'].append(curxl)
        outlist.append(entry)
    return outlist    