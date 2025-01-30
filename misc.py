# This file is part of the csbioscripts package
#
# Copyright (c) 2025 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

def lazycluster(listofthings, func, comparison, threshold):
    clusters = [[listofthings[0]]]
    for i in listofthings[1:]:
        merged = False
        for cj in clusters:
            for jmem in cj:
                if comparison(func(jmem, i), threshold):
                    cj.append(i)
                    merged = True
                    break
            if merged == True:
                break
        if merged == False:
            clusters.append([i])
    return clusters
                
    
    
    