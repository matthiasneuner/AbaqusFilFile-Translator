# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 14:54:55 2016

@author: c8441146
"""

import sys
import numpy as np

if __name__ == "__main__":
    
    localStateVars = 6
    globalStateVars = 2*6
    noIntPoint = 4
    alphaPIdx = 1
    alphaDIdx = 3
    
    totalStatPerInt = localStateVars + globalStateVars

with open('modLeon/exportFil.inp','w') as f:
    f.write("*defineElementType, element=U3, shape=quad4\n")
    f.write("*ensightPerNodeVariable, set=mainPart, exportName=nodeDisplacements, source=U, dimensions=3,\n")
    f.write("0,1\n")
    
    for i in range(noIntPoint):
        f.write("*ensightPerElementVariable, set=ALL, exportName=ALPHAP_IP"+str(i+1)+", source=SDV\n")
        f.write(str(totalStatPerInt*(i)+alphaPIdx)+'\n')
        f.write("*ensightPerElementVariable, set=ALL, exportName=ALPHAD_IP"+str(i+1)+", source=SDV\n")
        f.write(str(totalStatPerInt*(i)+alphaDIdx)+'\n')
        f.write("*ensightPerElementVariable, set=ALL, exportName=STRESS_IP"+str(i+1)+", source=SDV\n")
        formattedString = ''
        for j in range(6):
            formattedString += str(totalStatPerInt*(i)+localStateVars+j)
            if j<5:
                formattedString += ','
        f.write(formattedString+'\n')
    for i in range(noIntPoint):
        f.write("*ensightPerElementVariable, set=ALL, exportName=STRAIN_IP"+str(i+1)+", source=SDV\n")
        formattedString = ''
        for j in range(6):
            formattedString += str(totalStatPerInt*(i)+localStateVars+6+j)
            if j<5:
                formattedString += ','
        f.write(formattedString+'\n')
        
