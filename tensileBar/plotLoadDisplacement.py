#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:42:48 2017

@author: c8441146
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import itertools


if __name__ == "__main__":    
    marker = itertools.cycle((',', '+', '.', 'o', '*')) 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for suffix in ['_5el.csv', '_25el.csv']:
        print(suffix)
        time = np.loadtxt('time'+suffix)
        u = np.loadtxt('nodeDisplacements'+suffix)
        rf = np.loadtxt('reactionForces'+suffix)
        ax.plot(u[:,0], rf[:,0], label=suffix)#, marker=marker.next())
        
    plt.legend()
    plt.grid()
    plt.show()

#else:
#    xFile = sys.argv[1]
#    xIdx = int(sys.argv[2])
#    yFile = sys.argv[3]
#    yIdx = int(sys.argv[4])
#    
#    xData = np.loadtxt(xFile)[:,xIdx] if np.loadtxt(xFile).ndim>1 else np.loadtxt(xFile)
#    yData = np.loadtxt(yFile)[:,yIdx] if np.loadtxt(yFile).ndim>1 else np.loadtxt(yFile)
#    
#    plt.plot(xData, yData, label=yFile)
#    plt.legend()
#    plt.grid()
#    plt.show()