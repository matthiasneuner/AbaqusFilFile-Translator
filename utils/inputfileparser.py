# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 18:56:03 2015

@author: matthias
"""  

import numpy as np
from os.path import dirname, join
import textwrap

dTypes = {int : "integer",
          float: "float",
          str: "string",          
          "npInt": "numpy int array",
          }
    
typeMappings = {    '*defineElementType' : ('assign an ensight Shape to an Abaqus Element',
                                            {'element' : (str, 'Abaqus (User) Element'),
                                             'shape' :  (str, 'Ensight Shape'),
                                             }),
                    '*ensightCaseOptions' :  ('modify Ensight export options',
                                            {'discardTime' : (str, ''),
                                             }),

                    '*ensightPerNodeVariable':        ("define an Ensight per node variable for export",

                        {'set':             (str, "Abaqus element set") ,
                         'exportName':      (str, "export name of the variable"),
                         'source':          (str, 'Abaqus variable identifier'),
                         'dimensions':      (int, "(optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing components will be zero filled"),
                         'data':            (int, "indices of the Abaqus source Variable")}),
                                                       
                    '*ensightPerElementVariable':        ("define an Ensight per element variable for export",

                        {'set':             (str, "Abaqus element set") ,
                         'exportName':      (str, "export name of the variable"),
                         'source':          (str, 'Abaqus variable identifier'),
                         'dimensions':      (int, "(optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing components will be zero filled"),
                         'periodicalPattern':(int, "(optional), define a periodical pattern for extraction (e.g. for GaussPts)"),
                         'periodicalShift': (int, "(optional), define a periodical pattern for extraction (e.g. for GaussPts)"),
                         'data':            ("npInt", "indices of the Abaqus source Variable")}),
                                                          
                    '*csvPerElementOutput':        ("export element results to a csv file",
                        {'set':             (str, "Abaqus element set") ,
                         'exportName':      (str, "export name of the variable"),
                         'source':          (str, 'Abaqus variable identifier'),
                         'fmt':             (str, '(optional) formatter for csv output by np.savetxt'),
#                         'dimensions':      (int, "(optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing components will be zero filled"),
                         'periodicalPattern':(int, "(optional), define a periodical pattern for extraction (e.g. for GaussPts)"),
                         'periodicalShift': (int, "(optional), define a periodical pattern for extraction (e.g. for GaussPts)"),
                         'data':            ("npInt", "indices of the Abaqus source Variable")}),
                                                    
                    '*csvPerNodeOutput':        ("export node results to a csv file",
                        {'node':             (int, "(alternative) specify an Abaqus node label") ,
                         'nSet':            (str, '(alternative) specify an Abaqus node set'),
                         'exportName':      (str, "export name of the variable"),
                         'source':          (str, 'Abaqus variable identifier'),
                         'fmt':             (str, '(optional) formatter for .csv output by np.savetxt'),
                         'data':            ("npInt", "indices of the Abaqus source Variable")}),
                                                 
                      '*exportTimeHistory':        ("export the time history as csv file",
                            {'exportName':      (str, "export Name of the variable"),
                             }),
                               
                    '*include': ("(optional) load extra .inp file (fragment), use relative path to current .inp",
                        {'input':           (str, "filename")}),
                        

                }
                
def getMapType(kw, varName):
    kwDict = typeMappings.get(kw, (None,{}) )[1]    
    mType = kwDict.get(varName, [str])[0]
    return mType
    
def npPrepare(data):
    """ replaces 'x' in datalines with np.infs"""
    return data.replace("x", "inf")

def parseInputFile(fileName, currentKeyword = None, existingFileDict = None):
    """ parse an Abaqus like input file to generate an dictionary with its content"""
    if not existingFileDict:
        fileDict = { key : [] for key in typeMappings.keys()}
    else:
        fileDict = existingFileDict
    keyword = currentKeyword
    with open(fileName) as f:
        for l in f:
            lineElements = [x.strip() for x in l.split(",")]
            lineElements=list(filter(None,lineElements))
            if not lineElements or lineElements[0].startswith("**"):
                # line is a comment
                pass
            elif lineElements[0].startswith("*"):
                # line is keywordline
                lastkeyword = keyword
                keyword = lineElements[0]
                optionAssignments = lineElements[1:]
                
                objectentry = {}
                objectentry['data'] = []
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
                
            else:
                # line is a dataline
                data = lineElements
                mType = getMapType(keyword, "data")
                if mType is not None:
                    if mType == "npInt":
                        data = np.array([npPrepare(x) for x in data], dtype = np.int)
                    else:    
                        data = [mType(d) for d in data]
                fileDict[keyword][-1]['data'].append(data)
    
    return fileDict

def printKeywords():
    """ print the input file language set"""
    kwString = "    {:}    "
    kwDataString = "        {:22}{:20}"    
    
    wrapper = textwrap.TextWrapper(width=80,replace_whitespace=False)
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
    

def mergeNumpyDataLines(multiLineData):
    flattenedMatProps = [p for row in multiLineData for p in row]
    return np.array(flattenedMatProps, dtype=np.float)
    
def mergeDictDataLines(multiLineData):
    d = {key:val for (key, val) in multiLineData}
    return d
    
