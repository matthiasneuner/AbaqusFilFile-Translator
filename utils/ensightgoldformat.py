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
    def __init__(self,description="partDescription",partNumber=1, elements = None, nodes=None, nodeLabels=None):
        self.structureType = "coordinates"
        self.nodes = nodes    # 2D numpyArray
        self.nodeLabels = nodeLabels or []    #[integer]
        self.elements = elements or {}            #elementSet entries: {strElementType : [ ( intLabel = None, [nodeList] ) ]}
        self.description = description  #string, describing the part; max. 80 characters
        self.partNumber = partNumber 
        return
        
    def writeToFile(self, binaryFileHandle, printNodeLabels=True, printElementLabels=True):
#        print(self.nodes)
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
        
        if printNodeLabels:                 
            writeCInt(f, self.nodeLabels)
        writeCFloat(f, self.nodes.T)
                
        for elemType, elemList in self.elements.items():
#            if len(elemList) is not 0:
                
            writeC80(f, elemType)
            writeCInt(f, len(elemList))
            if printElementLabels:
                writeCInt(f, elemID)
            for elemID, nodes in elemList:             
                writeCInt(f, np.asarray(nodes, np.int32)[:]+1 )

class EnsightPointPart(EnsightUnstructuredPart):
    """derived UnstructuredPart to represent a single point"""
    def __init__(self, partNumber, coordinates, description="single_point"):
        super().__init__(description, partNumber, {"point" : [(1, [0])] }, coordinates, nodeLabels=[1])

class EnsightStructuredPart:
    """ define a structured part using a 3d array for the blockdimensions (for shell elements: one dimension may be 1),
        TODO: IMPLEMENT NODELABELS AND ELEMENTLABELS PRINTING
    """
    def __init__(self,description="partDescription",partNumber=1, blockDimensions=None, nodeCoordinates=None, nodeLabels=None):
        self.structureType = "block"
        self.nodeLabels = nodeLabels if nodeLabels is not None else []    #[integer]
        self.description = description  #string, describing the part; max. 80 characters
        self.partNumber = partNumber    #int of partnumber
        self.blockDimensions = blockDimensions if blockDimensions is not None else []
        self.nodeCoordinates = nodeCoordinates if nodeCoordinates is not None else []
        self.variables = {}
        return
        
    def writeToFile(self, binaryFileHandle, printNodeLabels=False, printElementLabels=False):
        f = binaryFileHandle  #shortcut to functions
        writeC80(f, 'part')
        writeCInt(f, self.partNumber)
        writeC80(f, self.description)
        writeC80(f, 'block')
        
        writeCInt(f, self.blockDimensions)
        writeCFloat(f, self.nodeCoordinates.T)
                
        #todo: implement nodeLabels and elementLabels printing(!)
                    
class EnsightTimeSet:
    """ defines a set which may be used by EnsightGeometry, EnsightStructuredPart, EnsightUnstructuredPart and is written into the case file"""
    def __init__(self,number=1, description="timeStepDesc",  fileNameStartNumber=0, fileNameNumberIncrement=1, timeValues = None):
        self.number = number
        self.description = description
        self.fileNameStartNumber = fileNameStartNumber
        self.fileNameNumberIncrement = fileNameNumberIncrement
        self.timeValues = timeValues if timeValues is not None else []
        
class EnsightGeometryTrend:
    """ container class for the time dependent evolution of the geometry,
        establishes the connection between the geometry entities and a EnsightTimeSet"""
    def __init__(self, ensightTimeSet, ensightGeometryList = None):
        self.timeSet = ensightTimeSet
        self.geometryList = ensightGeometryList if ensightGeometryList is not None else []
        
class EnsightGeometry:
    """ container class for one or more EnsightParts at a certain time state, handles also the file writing operation"""
    def __init__(self, name="geometry", descriptionLine1 ="", descriptionLine2 ="", ensightPartList = None , nodeIdOption="assign", elementIdOption="assign"):
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
        
    def writeToFile(self, fileHandle):
        f = fileHandle
        writeC80(f, self.description)
        for ensightPartID, elTypeDict in self.partsDict.items():          
            writeC80(f, 'part')
            writeCInt(f, ensightPartID)
            for elType, values in elTypeDict.items():
                writeC80(f, elType)
                writeCFloat(f, values.T)
            
class EnsightCase:
    """ master class for visualization export. 
    Minimum information is a list of timeSets, containing a minimum of 1 timeSet,
    and a geometryTrend instance, containing a minimum of one geometryInstance associated with the timeSet.
    This class also handles the file writing operations for all associated entities (geometry, parts, variables, casefile)
    """
    def __init__(self, ensightTimeSetList, ensightGeometryTrend , ensightVariableTrendsList = None):
        self.timeSetList = ensightTimeSetList
        self.geometryTrend = ensightGeometryTrend
        self.variableTrendList = ensightVariableTrendsList if ensightVariableTrendsList is not None else []
    
    def writeCaseToFiles(self, directory, casename, writeTransientSingleFiles=True):
        caseFileName = os.path.join(directory, casename + ".case")
        
        with open(caseFileName, mode="w") as cf:
            
            cf.write("FORMAT\n")
            cf.write("type: ensight gold\n")
        
            cf.write("TIME\n")
            for timeSet in self.timeSetList:
                cf.write("time set: "+str(timeSet.number)+" "+timeSet.description+"\n")
                cf.write("number of steps: "+str(len ( timeSet.timeValues))  + "\n")
                cf.write("filename start number: " + str(timeSet.fileNameStartNumber) +"\n")
                cf.write("filename increment: " +str(timeSet.fileNameNumberIncrement) +"\n")
                cf.write("time values: ")
                for timeVal in timeSet.timeValues:
                    cf.write(str(timeVal) +"\n")
                
            if writeTransientSingleFiles:
                fs = "1"
                cf.write("FILE\n")
                for timeSet in self.timeSetList:
                    cf.write("file set: {:} \n".format(timeSet.number))
                    cf.write("number of steps: {:} \n".format(len(timeSet.timeValues)))
            else:
                fs = " "
        
            
            cf.write("GEOMETRY\n")
            print(directory)
            print(casename)
            geometryFileName = os.path.join(directory, casename)       
            geometryTSn = self.geometryTrend.timeSet.number               
            geometryFSn = self.geometryTrend.timeSet.number               
            cf.write("model: {:} {:} {:} \n".format(geometryTSn, geometryFSn, casename+".geo***"))
            
            initFileNumber =  self.geometryTrend.timeSet.fileNameStartNumber
            numberIncrement = self.geometryTrend.timeSet.fileNameNumberIncrement
            
            currentFileNumber = initFileNumber
            
            if writeTransientSingleFiles:
                with open(geometryFileName + ".geo" + str(currentFileNumber).zfill(3),mode='wb') as f:
                    writeC80(f, 'C Binary')
                    for geom in self.geometryTrend.geometryList:
                        writeC80(f, 'BEGIN TIME STEP')
                        geom.writeToFile(f)
                        writeC80(f, 'END TIME STEP')
            else:
                for geom in self.geometryTrend.geometryList:
                    with open(geometryFileName + ".geo" + str(currentFileNumber).zfill(3),mode='wb') as f:
                        writeC80(f, 'C Binary')
                        geom.writeToFile(f)
                        currentFileNumber += numberIncrement

            if self.variableTrendList:
                cf.write("VARIABLE\n")
                for variableTrend in self.variableTrendList:
                    varFileName = os.path.join(directory, variableTrend.variableName)     
                    cf.write("{:}: {:} {:} {:} {:}.var***\n".format(
                        variableTrend.variableType,
                        variableTrend.timeSet.number,
                        variableTrend.timeSet.number,
                        variableTrend.description,
                        variableTrend.variableName))

                    initFileNumber = variableTrend.timeSet.fileNameStartNumber
                    numberIncrement= variableTrend.timeSet.fileNameNumberIncrement
                    currentFileNumber =  initFileNumber
                
                    if writeTransientSingleFiles:             
                        fileName = varFileName + ".var"+str(currentFileNumber).zfill(3)                    
                        with open(fileName, mode="wb") as f:
                            for variable in variableTrend.variableList:
                                writeC80(f, 'BEGIN TIME STEP')
                                variable.writeToFile(f)
                                writeC80(f, 'END TIME STEP')
                    else:
                        for variable in variableTrend.variableList:
                            fileName = varFileName + ".var"+str(currentFileNumber).zfill(3)
                            with open(fileName, mode="wb") as f:
                                variable.writeToFile(f)
                                currentFileNumber += numberIncrement                        
                     
class EnsightChunkWiseCase:
    
    def __init__(self, directory, caseName, writeTransientSingleFiles = True):
        self.directory = directory
        self.caseName = caseName
        self.caseFileNamePrefix = os.path.join(directory, caseName)  
        print(directory)
        print(caseName)
        print(self.caseFileNamePrefix)
        self.writeTransientSingleFiles = writeTransientSingleFiles
        self.timeAndFileSets = {}
        self.geometryTrends = {}
        self.variableTrends = {}

    def setCurrentTime(self, timeAndFileSetNumber, timeValue):
        if not timeAndFileSetNumber in self.timeAndFileSets:
            self.timeAndFileSets[timeAndFileSetNumber] = EnsightTimeSet(timeAndFileSetNumber, 'noDesc', 0,1)
        tfSet = self.timeAndFileSets[timeAndFileSetNumber]
        tfSet.timeValues.append(timeValue)
        
    def writeGeometryTrendChunk(self, ensightGeometry, timeAndFileSetNumber=1):
        
        fileName = ('{:}'*3).format(self.caseFileNamePrefix,
                                    ensightGeometry.name,
                                    ".geo",
#                                    str(0).zfill(3)
                                    )
        print(fileName) 
        print(self.caseFileNamePrefix)
        if not ensightGeometry.name in self.geometryTrends:
            self.geometryTrends[ensightGeometry.name] = timeAndFileSetNumber
            with open(fileName ,mode='wb') as f:
                writeC80(f, 'C Binary')
                
        if self.writeTransientSingleFiles:
                with open(fileName, mode='ab') as f:
                    writeC80(f, 'BEGIN TIME STEP')
                    ensightGeometry.writeToFile(f)
                    writeC80(f, 'END TIME STEP')
        
    def writeVariableTrendChunk(self, ensightVariable, timeAndFileSetNumber=2):
        
        fileName = ('{:}'*3).format(self.caseFileNamePrefix,
                                    ensightVariable.name,
                                    ".var",
#                                    str(0).zfill(3)
                                    )
        if not ensightVariable.name in self.variableTrends:
            self.variableTrends[ensightVariable.name] = timeAndFileSetNumber, ensightVariable.varType
            with open(fileName ,mode='wb') as f:
                writeC80(f, 'C Binary')
                
        if self.writeTransientSingleFiles:
                with open(fileName, mode='ab') as f:
                    writeC80(f, 'BEGIN TIME STEP')
                    ensightVariable.writeToFile(f)
                    writeC80(f, 'END TIME STEP')
                    
    def finalize(self):
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
                for timeVal in timeSet.timeValues:
                    cf.write(str(timeVal) +"\n")
                
            if self.writeTransientSingleFiles:
                cf.write("FILE\n")
                for timeSet in self.timeAndFileSets.values():
                    cf.write("file set: {:} \n".format(timeSet.number))
                    cf.write("number of steps: {:} \n".format(len(timeSet.timeValues)))
            
            cf.write("GEOMETRY\n")
            for geometryName, tAndFSetNum in self.geometryTrends.items():
#                geometryTSn = tAndFSetNum              
#                geometryFSn = tAndFSetNum              
                geometryTSn = ''
                geometryFSn = ''              
                cf.write("model: {:} {:} {:} \n".format(geometryTSn, geometryFSn, self.caseFileNamePrefix+geometryName+".geo"))
                
            cf.write("VARIABLE\n")
            for variableName, (tAndFSetNum, variableType) in self.variableTrends.items():
                cf.write("{:}: {:} {:} {:} {:}.var\n".format(
                    variableType,
                    tAndFSetNum,
                    tAndFSetNum,
                    variableName,
                   self.caseFileNamePrefix+ variableName))

                    
