# This file is part of the csbioscripts package
#
# Copyright (c) 2025 - Topf Lab, Leibniz-Institut für Virologie
# Hamburg, Germany.
#
# This module was developed by:
#   Karen Manalastas-Cantos    <karen.manalastas-cantos AT cssb-hamburg.de>

import numpy as np

def lazycluster(listofthings, func, comparison, threshold):
    clusters = [[listofthings[0]]]
    for i in listofthings[1:]:
        merged = False
        icompared = False
        for cj in clusters:
            for jmem in cj:
                if func(jmem, i) != None:   # function worked
                    icompared = True
                    if comparison(func(jmem, i), threshold):
                        cj.append(i)
                        merged = True
                        break
            if merged == True:
                break
        if merged == False and icompared == True:   # only add to list if function did not fail
            clusters.append([i])        
    return clusters

def distancematrix(a, b):
    distmat = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return distmat

                
    
    
    