# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 09:18:51 2015

@author: c8441141
"""
import numpy as np
import os

def writeCFloat(f, ndarray):
    np.asarray(ndarray, dtype=np.float32).tofile(f)
def writeCInt(f, ndarray):
    np.asarray(ndarray, dtype=np.int32).tofile(f)
def writeC80(f, string):
    np.asarray(string, dtype='a80').tofile(f)

variableTypes = {"scalar" : 1,
                 "vector" : 3,
                 "tensor": 9}

ensightPerNodeVariableTypes = {
                               1 : 'scalar per node',
                               3 : 'vector per node',
                               6 : 'tensor per node',
                               9 : 'tensor9 per node'}  
                                   
ensightPerElementVariableTypes = {
                                  1 : 'scalar per element',
                                  3 : 'vector per element',
                                  6 : 'tensor per element',
                                  9 : 'tensor9 per element'}
                 
class EnsightUnstructuredPart:
    """ define an unstructured part, by a list of nodes and a dictionary of elements.
    Each dictionary entry consists of a list of tuples of elementlabel and nodelist: 
    {strElementType : [ ( intLabel = None, [nodeList] ) ]}"""
    def __init__(self, description, partNumber, elements, nodes, nodeLabels):
        self.structureType = "coordinates"
        self.nodes = nodes    # 2D numpyArray
        self.nodeLabels = nodeLabels    #[integer]
        self.elements = elements             #elementSet entries: {strElementType : [ ( intLabel = None, [nodeList] ) ]}
        self.description = description  #string, describing the part; max. 80 characters
        self.partNumber = partNumber 
        
    def writeToFile(self, binaryFileHandle, printNodeLabels=True, printElementLabels=True):
        if len(self.nodes.shape) < 2:
            pass
        elif self.nodes.shape[1] < 3:
            print('Ensight Unstructered Part: Extending node coordinate array to 3 dims')
            self.nodes = np.hstack( (self.nodes, np.zeros((self.nodes.shape[0],3-self.nodes.shape[1])) ))
        nNodes = self.nodes.shape[0]
        f = binaryFileHandle  #shortcut to functions
        
        writeC80(f, 'part')
        writeCInt(f, self.partNumber)
        writeC80(f, self.description)
        writeC80(f, 'coordinates')
        writeCInt(f, nNodes)
        
        # nodes
        if printNodeLabels:             
            writeCInt(f, self.nodeLabels)
        writeCFloat(f, self.nodes.T)
        
        # elements
        for elemType, elemList in self.elements.items():
            writeC80(f, elemType)
            writeCInt(f, len(elemList))
            if printElementLabels:
                for elemID, nodes in elemList:             
                    writeCInt(f, elemID)
            for elemID, nodes in elemList:             
                writeCInt(f, np.asarray(nodes, np.int32)[:]+1 )

                    
class EnsightTimeSet:
    """ defines a set which may be used by EnsightGeometry, EnsightStructuredPart, EnsightUnstructuredPart and is written into the case file"""
    def __init__(self,number=1, description="timeStepDesc",  fileNameStartNumber=0, fileNameNumberIncrement=1, timeValues = None):
        self.number = number
        self.description = description
        self.fileNameStartNumber = fileNameStartNumber
        self.fileNameNumberIncrement = fileNameNumberIncrement
        self.timeValues = timeValues if timeValues is not None else []
        
class EnsightGeometry:
    """ container class for one or more EnsightParts at a certain time state, handles also the file writing operation"""
    def __init__(self, name="geometry", descriptionLine1 ="", descriptionLine2 ="", ensightPartList = None , nodeIdOption="given", elementIdOption="given"):
        self.name = name
        self.descLine1 = descriptionLine1
        self.descLine2 = descriptionLine2
        self.partList = ensightPartList if ensightPartList is not None else []
        self.nodeIdOption = nodeIdOption
        self.elementIdOption = elementIdOption
    
    def writeToFile(self, fileHandle):
        f = fileHandle
        writeC80(f, self.descLine1)
        writeC80(f, self.descLine2)
        writeC80(f, "node id "+self.nodeIdOption)
        writeC80(f, "element id "+self.elementIdOption)
        
        if self.nodeIdOption == "given" or self.nodeIdOption == "ignore":
            printNodeLabels = True
        else:# assign or off
            printNodeLabels = False
            
        if self.elementIdOption == "given" or self.nodeIdOption == "ignore":
            printElementLabels = True
        else: # assign or off 
            printElementLabels = False
            
        for part in self.partList:
            part.writeToFile(f, printNodeLabels, printElementLabels)

class EnsightVariableTrend:
    """ container class for the time dependent evolution of one variable,
        establishes the connection between EnsightVariable entities and a EnsighTimeSet"""
    def __init__(self, ensightTimeSet, variableName,  ensightVariableList = None, variableType = "scalar per node", description = "variableTrendDescription"):
        self.timeSet = ensightTimeSet
        self.variableName = variableName
        self.variableList = ensightVariableList if ensightVariableList is not None else []
        self.variableType = variableType
        self.description = description
        
class EnsightPerNodeVariable:
    """ container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }"""
    def __init__(self, name, variableDimension, ensightPartsDict = None) :
        self.name = name
        self.description = name
        self.partsDict = ensightPartsDict or {} # { EnsightPart: np.array(variableValues) }
        self.variableDimension = variableDimension
        self.varType = ensightPerNodeVariableTypes[variableDimension]
        
    def writeToFile(self, fileHandle, ):
        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, (structureType,  values) in self.partsDict.items():          
            writeC80(f, 'part')
            writeCInt(f, ensightPartID)
            writeC80(f, structureType)
            writeCFloat(f, values.T)
            if values.shape[1] < self.variableDimension:
                writeCFloat(f, np.zeros( (values.shape[0],self.variableDimension - values.shape[1])))
            
class EnsightPerElementVariable:
    """ container class for data for one certain variable, defined for one or more parts (classification by partID), at a certain time state.
    For each part the structuretype ("coordinate" or "block") has to be defined.
    Each part-variable assignment is defined by a dictionary entry of type: { EnsightPart: np.array(variableValues) }"""
    def __init__(self, name, variableDimension, ensightPartsDict = None, ) :
        self.name = name
        self.description = name
        self.partsDict = ensightPartsDict or {} # { EnsightPart: np.array(variableValues) }
        self.varType = ensightPerElementVariableTypes[variableDimension]
        self.variableDimension = variableDimension
        
    def writeToFile(self, fileHandle):
        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, elTypeDict in self.partsDict.items():          
            writeC80(f, 'part')
            writeCInt(f, ensightPartID)
            for elType, values in elTypeDict.items():
                writeC80(f, elType)
                writeCFloat(f, values.T)
            if values.shape[1] < self.variableDimension:
                writeCFloat(f, np.zeros( (values.shape[0],self.variableDimension - values.shape[1])))
                     
class EnsightChunkWiseCase:
    
    def __init__(self, directory, caseName, writeTransientSingleFiles = True):
        self.directory = directory
        self.caseName = caseName
#        self.caseFileNamePrefix = os.path.join(directory, caseName)  
        self.caseFileNamePrefix = caseName
        self.writeTransientSingleFiles = writeTransientSingleFiles
        self.timeAndFileSets = {}
        self.geometryTrends = {}
        self.variableTrends = {}
        self.fileHandles = {}
        self.currentTime = 0.0
        
    def setCurrentTime(self,  timeValue):
        self.currentTime = timeValue
        
    def writeGeometryTrendChunk(self, ensightGeometry, timeAndFileSetNumber=1):
        
        if timeAndFileSetNumber != None:
        
            if not timeAndFileSetNumber in self.timeAndFileSets:
                self.timeAndFileSets[timeAndFileSetNumber] = EnsightTimeSet(timeAndFileSetNumber, 'timeset', 0,1)
                self.timeAndFileSets[timeAndFileSetNumber].timeValues.append ( self.currentTime )
                
            elif self.currentTime > self.timeAndFileSets[timeAndFileSetNumber].timeValues[-1]:
                self.timeAndFileSets[timeAndFileSetNumber].timeValues.append(self.currentTime)
        
        if ensightGeometry.name not in self.fileHandles:
            fileName = ('{:}'*3).format(self.caseFileNamePrefix,
                                        ensightGeometry.name,
                                        ".geo",
                                        )
            self.fileHandles[ensightGeometry.name] = open(fileName, mode='wb')
        
        f = self.fileHandles[ensightGeometry.name]
        
        if not ensightGeometry.name in self.geometryTrends:
            self.geometryTrends[ensightGeometry.name] = timeAndFileSetNumber
            writeC80(f, 'C Binary')
                
        if self.writeTransientSingleFiles:
            writeC80(f, 'BEGIN TIME STEP')
            ensightGeometry.writeToFile(f)
            writeC80(f, 'END TIME STEP')
        
    def writeVariableTrendChunk(self, ensightVariable, timeAndFileSetNumber=2):
        
        if not timeAndFileSetNumber in self.timeAndFileSets:
            self.timeAndFileSets[timeAndFileSetNumber] = EnsightTimeSet(timeAndFileSetNumber, 'timeset', 0,1)
            self.timeAndFileSets[timeAndFileSetNumber].timeValues.append ( self.currentTime )
            
        elif self.currentTime > self.timeAndFileSets[timeAndFileSetNumber].timeValues[-1]:
            self.timeAndFileSets[timeAndFileSetNumber].timeValues.append(self.currentTime)
        
        if ensightVariable.name not in self.fileHandles:
            fileName = ('{:}'*3).format(self.caseFileNamePrefix,
                            ensightVariable.name,
                            ".var")
            self.fileHandles[ensightVariable.name] = open(fileName, mode='wb')
        
        f = self.fileHandles[ensightVariable.name]
        
        if not ensightVariable.name in self.variableTrends:
            self.variableTrends[ensightVariable.name] = timeAndFileSetNumber, ensightVariable.varType
            writeC80(f, 'C Binary')
                
        if self.writeTransientSingleFiles:
            writeC80(f, 'BEGIN TIME STEP')
            ensightVariable.writeToFile(f)
            writeC80(f, 'END TIME STEP')
                    
    def finalize(self, discardTimeMarks=False, closeFileHandles=True):
        
        if closeFileHandles:
            for f in self.fileHandles.values():
                f.close()
        
        caseFName = self.caseFileNamePrefix+'.case'
        with open(caseFName ,mode='w') as cf:
            cf.write("FORMAT\n")
            cf.write("type: ensight gold\n")
        
            cf.write("TIME\n")
            for setNum, timeSet in self.timeAndFileSets.items():
                cf.write("time set: "+str(setNum)+" noDesc\n")
                cf.write("number of steps: "+str(len (timeSet.timeValues))  + "\n")
                cf.write("filename start number: " + str(timeSet.fileNameStartNumber) +"\n")
                cf.write("filename increment: " +str(timeSet.fileNameNumberIncrement) +"\n")
                cf.write("time values: ")
                for i, timeVal in enumerate(timeSet.timeValues):
                    if discardTimeMarks:
                        cf.write('{:}'.format(i) +"\n")
                    else:
                        cf.write('{:1.10f}'.format(timeVal) +"\n")
                
            if self.writeTransientSingleFiles:
                cf.write("FILE\n")
                for timeSet in self.timeAndFileSets.values():
                    cf.write("file set: {:} \n".format(timeSet.number))
                    cf.write("number of steps: {:} \n".format(len(timeSet.timeValues)))
            
            cf.write("GEOMETRY\n")
            for geometryName, tAndFSetNum in self.geometryTrends.items():
                cf.write("model: {:} \n".format(self.caseFileNamePrefix+geometryName+".geo"))
                
            cf.write("VARIABLE\n")
            for variableName, (tAndFSetNum, variableType) in self.variableTrends.items():
                cf.write("{:}: {:} {:} {:} {:}.var\n".format(
                    variableType,
                    tAndFSetNum,
                    tAndFSetNum,
                    variableName,
                   self.caseFileNamePrefix+ variableName))

                    
