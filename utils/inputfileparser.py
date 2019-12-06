"""
Created on Tue Sep 22 18:56:03 2015

Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""  

import numpy as np
from os.path import dirname, join
import textwrap
import shlex

dTypes = {int : "integer",
          float: "float",
          str: "string",}
    
typeMappings = {    '*defineElementType' :                  ('assign an ensight Shape to an Abaqus Element',
                        {'element' :                        (str, 'Abaqus (User) Element'),
                            'shape' :                       (str, 'Ensight Shape, can be any of: point bar2 bar3 tria3 tria6 quad4 quad8 tetra4 tetra10 pyramid5 pyramid13 penta6 penta15 hexa8 hexa20 nsided nfaced '), }),

                    '*ensightCaseOptions' :                 ('modify Ensight export options',
                        {'discardTime' :                    (str, 'discard Time values and replace by enumeration of time steps'), }),

                    '*ensightPerNodeVariable':              ("define an Ensight per node variable for export",
                        {'set':                             (str, "Abaqus element set") ,
                         'exportName':                      (str, "export name of the variable"),
                         'source':                          (str, 'Abaqus variable identifier'),
                         'dimensions':                      (int, "(optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing components will be zero filled"),
                         'timeSet':                         (int, "(optional), define a timeset, for 'different' timelines"),
                         'values':                          (str, "(optional), define a index/slice to extract a subarray from the total result array (per Element)"),
                         'f(x)':                            (str, "(optional), apply a mathematical/array expression on the result array (per Element, slow!)"),
                         'fillMissingValues':               (float, "(optional), fill missing nodal values with a constant values, requires specified dimensions (slow!)"),
                         }),
                                                       
                    '*ensightPerElementVariable':           ("define an Ensight per element variable for export",
                        {'set':                             (str, "Abaqus element set") ,
                         'exportName':                      (str, "export name of the variable"),
                         'source':                          (str, 'Abaqus variable identifier'),
                         'dimensions':                      (int, "(optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing components will be zero filled"),
                         'nIntegrationPoints':              (int, "(optional), define a periodical pattern for a repeatet extraction (e.g. for results @ GaussPts)"),
                         'integrationPointDataDistance':    (int, "(optional), define a periodical pattern: initial constant offset"),
                         'integrationPointDataOffset':      (int, "(optional), define a periodical pattern: offset between extraction points"),
                         'timeSet':                         (int, "(optional), define a timeset, for 'different' timelines"),
                         'values':                          (str, "(optional), define a index/slice to extract a subarray from the total result array (per Element)"),
                         'f(x)':                            (str, "(optional), apply a mathematical/array expression on the result array (per Element, slow!)")}),
                               
                    '*include':                             ("(optional) load extra .inp file (fragment), use relative path to current .inp",
                        {'input':                           (str, "filename")}),
                }
                
def getMapType(kw, varName):
    kwDict = typeMappings.get(kw, (None,{}) )[1]    
    mType = kwDict.get(varName, [str])[0]
    return mType

def parseInputFile(fileName, currentKeyword = None, existingFileDict = None):
    """ parse an Abaqus like input file to generate an dictionary with its content"""
    if not existingFileDict:
        fileDict = { key : [] for key in typeMappings.keys()}
    else:
        fileDict = existingFileDict
    keyword = currentKeyword
    with open(fileName) as f:
        for l in f:
            lexer = shlex.shlex(l.strip(), posix=True)
            lexer.whitespace_split = True
            lexer.whitespace = ','
            lineElements = [x.strip() for x in lexer]
            if not lineElements or lineElements[0].startswith("**"):
                # line is a comment
                pass
            elif lineElements[0].startswith("*"):
                # line is keywordline
                lastkeyword = keyword
                keyword = lineElements[0]
                optionAssignments = lineElements[1:]
                
                objectentry = {}
                objectentry['inputFile'] = fileName # save also the filename of the original inputfile!
                    
                for ass in optionAssignments:
                    opts = ass.split("=")
                    optKey = opts[0].strip()
                    val = opts[1].strip()
                    mType = getMapType(keyword, optKey)
                    if mType is not None:
                        objectentry[optKey] = mType(val)
                    else:
                        objectentry[optKey] = val
                        
                #special treatment for *include:
                if keyword == "*include":
                    includeFile = objectentry['input']
                    parseInputFile(join(dirname(fileName), 
                                                          includeFile), 
                                                     currentKeyword=lastkeyword,
                                                     existingFileDict=fileDict)
                    keyword = lastkeyword
                else:
                    fileDict[keyword].append(objectentry)
    return fileDict

def printKeywords():
    """ print the input file language set"""
    kwString = "    {:}    "
    kwDataString = "        {:30}{:14}"    
    
    wrapper = textwrap.TextWrapper(width=120,replace_whitespace=False)
    for kw, (kwDoc,optiondict) in sorted(typeMappings.items()):
        wrapper.initial_indent = kwString.format(str(kw))
        wrapper.subsequent_indent = " "*len(wrapper.initial_indent)
        print(wrapper.fill(kwDoc))
        print('')
        
        for key in sorted(optiondict.keys()):
            optionName = key
            dType, description = optiondict[key]
            wrapper.initial_indent = kwDataString.format(str(optionName),dTypes[dType])
            wrapper.subsequent_indent = " "*len(wrapper.initial_indent)
            print(wrapper.fill(description))
        print("\n")
