"""
Created on Fri Oct 28 13:43:53 2016

Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np
from collections import defaultdict
import utils.ensightgoldformat as es

def filInt(word):
    return word.view('<i8').ravel()
def filString(word):
    return word.view('a8')
def filStrippedString(word):
    return filString(word).tostring().decode('utf-8').strip()
def filDouble(word):
    return word.view('<d').ravel()
def filFlag(word):
    return word[0:4].view('<i')[0]

def makeExtractionFunction(expression, symbol='x'):
    """ make a simple f(x) expression from string"""
    return lambda x : eval ( expression, globals(), { symbol : x} )

def sliceFromString(string, shift=0):
    """ generate a slice from a string, which can represent a slice or an index"""
    if ':' in string:
        a, b = string.split(':')
        return slice ( int(a) + shift, int(b)  +shift )
    else:
        return slice(int(string) + shift, int(string)+1 +shift )
    
class Node:
    def __init__(self, label=None, coords=None):
        self.label = label
        self.coords = coords

class Element:
    def __init__(self, label, shape, nodes):
        self.label = label
        self.shape = shape
        self.nodes = nodes

class ElSet:
    def __init__(self, name, elements, ensightPartID = None):
        self.name  = name
        self.ensightPartID = ensightPartID
        
        self.elements = defaultdict(list)
        self.appendElements(elements)
            
        self._elementLabels = None
        self._reducedNodes = None
        self._reducedElements = None
        self._reducedNodeCoords3D = None
        
    def appendElements(self, elements):
        for element in elements:
            self.elements[element.shape].append(element)
            
    def getEnsightCompatibleReducedNodes(self, ):
        self._reducedNodes = dict([ (node.label, node) for elementsByShape in self.elements.values() 
                                                            for element in elementsByShape 
                                                                for node in element.nodes ])
        return self._reducedNodes
    
    def getEnsightCompatibleReducedNodeCoords(self,):
        self._reducedNodeCoords3D = np.asarray([node.coords for node in self.getEnsightCompatibleReducedNodes().values()])
        return self._reducedNodeCoords3D
    
    def getEnsightCompatibleElementNodeIndices(self,):
        self._reducedNodeIndices = {node : i for (i, node) in enumerate(self.getEnsightCompatibleReducedNodes().keys()) }
        self._reducedElements = dict()
        for eShape, elements in self.elements.items():
            self._reducedElements[eShape] = [ (e.label , [ self._reducedNodeIndices[n.label] for n in e.nodes] ) for e in elements  ]
        return self._reducedElements
    
class NSet:
    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes
        
    def appendNodes(self, nodes):
        self.nodes += nodes

class ExportJob:
    def __init__(self, exportName, dimensions, timeSetID, perSetJobs = None, writeEmptyTimeSteps=True ):
        self.exportName = exportName
        self.dimensions = dimensions
        self.timeSetID = timeSetID
        self.perSetJobs = perSetJobs or {}
        self.writeEmptyTimeSteps = writeEmptyTimeSteps
    
class PerSetJob:
    def __init__(self, setName, variableSource, extractionSlice=False, extractionFunction=False, offset=False, fillMissingValues=False):
        self.setName = setName
        self.extractionSlice = extractionSlice
        self.extractionFunction = extractionFunction
        self.offset = offset
        self.source = variableSource
        self.fillMissingValues = fillMissingValues

class ExportEngine:
                        
    def __init__(self, exportJobs, caseName):
        
        self.perElementJobs = self.collectPerElementJobs( exportJobs.get('*ensightPerElementVariable', []) )
        self.perNodeJobs    = self.collectPerNodeJobs (  exportJobs.get('*ensightPerNodeVariable', []) )
        self.ensightElementTypeMappings = { x['element'] : x['shape'] for x in exportJobs.get('*defineElementType', [])}
            
        self.allNodes = defaultdict(Node)
        
        self.allElements = {}
        self.nSets = {}
        self.elSets = {}
        
        self.currentState = 'model setup'
        self.currentIncrement = {}
        self.currentAbqElSet = None
        self.currentSetName = 'ALL'
        self.currentIpt = 1
        self.nIncrements = 0
        self.timeHistory = []
        self.labelCrossReferences = {}
        self.ensightCase = es.EnsightChunkWiseCase('.', caseName)
        self.ensightCaseDiscardTimeMarks = False
        
        self.knownRecords = { 
            1: ('Element header record', self.elementHeaderRecord),
            5: ('SDV output', lambda x: self.handlePerElementOutput(x, 'SDV')),
            11: ('S output', lambda x : self.handlePerElementOutput(x, 'S')),
            21: ('E output', lambda x : self.handlePerElementOutput(x, 'E')),
            101: ('U output', lambda x : self.handlePerNodeOutput(x, 'U')),
            102: ('V output', lambda x : self.handlePerNodeOutput(x, 'V')),
            103: ('A output', lambda x : self.handlePerNodeOutput(x, 'A')),
            104: ('RF output', lambda x : self.handlePerNodeOutput(x, 'RF')),
            201: ('NT output', lambda x : self.handlePerNodeOutput(x, 'NT')),
            1501: ('Surface definition header', self.surfaceDefHeader),
            1502: ('Surface facet',  lambda x : None),
            1900: ('element definition', self.addElement),
            1901: ('node definition', self.addNode),
            1902: ('active dof', lambda x : None),
            1911: ('output request definition', self.outputDefinition),
            1921: ('heading', lambda x : None),
            1922: ('heading', lambda x : None),
            1931: ('node set definition', self.addNodeset),
            1932: ('node set definition cont.', self.contAddNodeset),
            1933: ('element set definition', self.addElset),
            1934: ('element set definition cont.', self.contAddElset),
            1940: ('label cross reference', self.addLabelCrossReference),
            2000: ('start increment', self.addIncrement),
            2001: ('end increment', self.finishAndParseIncrement),}
    
    def computeRecord(self, recordLength, recordType, recordContent):
        if recordType in self.knownRecords:
            doc, action = self.knownRecords[recordType]     
            action(recordContent)
            return True
        else:
            print("{:<20}{:>6}{:>10}{:>4}".format('unknown record:',recordType, ' of length', recordLength))
            return False
            
    def finishAndParseIncrement(self, recordContent):
        if self.currentState == 'model setup':
            
            # replace all label references by the respective labels            
            for key, label in self.labelCrossReferences.items():
                strKey =  key #str(key)
                if strKey in self.nSets:
                    self.nSets[label] = self.nSets[strKey]
                    self.nSets[label].name = label
                    del self.nSets[strKey]
                if strKey in self.elSets:
                    self.elSets[label] = self.elSets[strKey]
                    self.elSets[label].name = label
                    del self.elSets[strKey]            
            
            geometryTimesetNumber = None 
            geometry = self.createEnsightGeometryFromModel()
            self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)
            
        elif self.currentState == 'surface definition':
            pass
            
        elif self.currentState == 'increment parsing':
            self.nIncrements +=1
            self.ensightCase.setCurrentTime(self.currentIncrement['tTotal'])
            self.timeHistory.append(self.currentIncrement['tTotal'])
            
            print('parsing increment {:>5}  tTotal:{:>16.5f}'.format(self.nIncrements,self.currentIncrement['tTotal']))
            
            for exportJob in self.perNodeJobs.values():
                    enSightVar = self.createEnsightPerNodeVariableFromPerNodeJob(exportJob)
                    if enSightVar:
                        self.ensightCase.writeVariableTrendChunk(enSightVar, exportJob.timeSetID)
                        del enSightVar
                                        
            for exportJob in self.perElementJobs.values():
                
                    enSightVar = self.createEnsightPerElementVariableFromPerElementJob(exportJob)
                    if enSightVar:
                        self.ensightCase.writeVariableTrendChunk(enSightVar, exportJob.timeSetID)
                        del enSightVar
            
            if self.nIncrements % 10 == 0:
                # intermediate saving ...
                self.ensightCase.finalize(discardTimeMarks = self.ensightCaseDiscardTimeMarks, closeFileHandles=False)
                
            del self.currentIncrement
            self.currentIncrement = {}
        
    def finalize(self):
        self.ensightCase.finalize(discardTimeMarks = self.ensightCaseDiscardTimeMarks)
    
    def createEnsightGeometryFromModel(self):
        partList = []
        partNumber = 1
        
        ALLSet = ElSet('ALL', self.allElements.values())
        self.elSets['ALL'] = ALLSet
        
        for elSet in self.elSets.values():

            elSetPart = es.EnsightUnstructuredPart(elSet.name, 
                                                   partNumber, 
                                                   elSet.getEnsightCompatibleElementNodeIndices(), 
                                                   elSet.getEnsightCompatibleReducedNodeCoords(), 
                                                   list(elSet.getEnsightCompatibleReducedNodes().keys()))
            elSet.ensightPartID = partNumber
            partList.append(elSetPart)
            partNumber +=1
            
        geometry = es.EnsightGeometry('geometry', '-', '-', partList, 'given', 'given')
        return geometry

    def collectPerElementJobs(self, entries):
        """ Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets """
        jobs = {}
        for entry in entries:
            source = entry['source']
            setName = entry['set']
            nPattern = entry.get('nIntegrationPoints', False )
            for i in range(nPattern):
                exportName = '{:}_{:}'.format(entry['exportName'], i)
                if source == 'SDVUEL':
                    nPatternInitialOffset = entry.get('integrationPointDataOffset', 0)
                    nPatternDistance = entry.get('integrationPointDataDistance', 0)
                    perSetJob = PerSetJob ( setName,
                                        variableSource = 'SDV@1',
                                        extractionSlice = sliceFromString ( entry['values']) if 'values' in entry else False,
                                        extractionFunction = makeExtractionFunction ( entry['f(x)']) if 'f(x)' in entry else False,
                                        offset = nPatternDistance * i + nPatternInitialOffset if nPattern else False,
                                        )
                else:
                    perSetJob = PerSetJob ( setName,
                                        variableSource = '{:}@{:}'.format(entry['source'], i+1 ),
                                        extractionSlice = sliceFromString ( entry['values']) if 'values' in entry else False,
                                        extractionFunction = makeExtractionFunction (entry['f(x)']) if 'f(x)' in entry else False,
                                        offset=False,)
                
                dimensions = entry.get('dimensions', False)
                timeSet = entry.get('timeSet', 1)
                
                if exportName not in jobs:
                    jobs[exportName] = ExportJob ( exportName, dimensions, timeSet, None, True)
                    
                jobs[exportName].perSetJobs[setName] = perSetJob

        return jobs

    def collectPerNodeJobs(self, entries):
        """ Collect all defined per node jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets """
        jobs = {}
        for entry in entries:

            setName = entry['set']
            exportName = entry['exportName'] 
            
            perSetJob = PerSetJob ( setName,
                                    entry['source'],
                                    extractionSlice = sliceFromString ( entry['values']) if 'values' in entry else False,
                                    extractionFunction = makeExtractionFunction ( entry['f(x)']) if 'f(x)' in entry else False,
                                    offset = False,
                                    fillMissingValues = entry.get('fillMissingValues', False),
                                    )
            
            dimensions = entry.get('dimensions', False)
            timeSet = entry.get('timeSet', 1)
            
            if exportName not in jobs:
                jobs[exportName] = ExportJob ( exportName, dimensions, timeSet, None, True )
                
            jobs[exportName].perSetJobs[setName] = perSetJob

        return jobs 
    
    def createEnsightPerNodeVariableFromPerNodeJob(self, exportJob ):
        
        partsDict = {}
        for setName, perSetJob in exportJob.perSetJobs.items():
            elSet = self.elSets[setName]

            if type(perSetJob.fillMissingValues) == float:
                results = np.full ( ( len ( elSet.getEnsightCompatibleReducedNodes()), exportJob.dimensions ) , perSetJob.fillMissingValues) 
                for i, node in enumerate( elSet.getEnsightCompatibleReducedNodes().keys() ):
                    if node in self.currentIncrement['nodeResults'][perSetJob.source]:
                        vals = self.currentIncrement['nodeResults'][perSetJob.source][node]
                        results[i, 0: vals.shape[0]] = vals
            
            else:
                results = np.asarray([ self.currentIncrement['nodeResults']
                                                                    [perSetJob.source]
                                                                    [node]
                                                                     for node in elSet.getEnsightCompatibleReducedNodes().keys() ] )

            if perSetJob.extractionFunction:
                results = np.apply_along_axis( perSetJob.extractionFunction, axis = 1, arr = results[:, perSetJob.offset :]  )
                results = np.reshape ( results, ( results.shape[0], -1))

            if perSetJob.extractionSlice:
                results = results[:, perSetJob.extractionSlice]

            variableLength = results.shape[1]

            partsDict[elSet.ensightPartID] =  ('coordinates', results)

        # determine the dimension of the variable .. either it is specified or the length of the variable is taken
        variableDimension = exportJob.dimensions or variableLength
            
        if partsDict or  exportJob.writeEmptyTimeSteps:
            return  es.EnsightPerNodeVariable(exportJob.exportName, variableDimension, partsDict)
        else:
            return None
        
        
    def createEnsightPerElementVariableFromPerElementJob(self, exportJob):
        
        partsDict = {}
        for setName, perSetJob in exportJob.perSetJobs.items():
            elSet = self.elSets[setName]
            resultLocation = perSetJob.source
            
            incrementVariableResults = self.currentIncrement['elementResults'][resultLocation][setName]
            incrementVariableResultsArrays = {}
            for ensElType, elDict in incrementVariableResults.items():

                results = np.asarray([ elDict[el.label] for el in elSet.elements[ensElType] ])
                
                if perSetJob.offset:
                    results = results[:, perSetJob.offset :]

                if perSetJob.extractionFunction:
                    results = np.apply_along_axis( perSetJob.extractionFunction, axis = 1, arr = results )
                    results = np.reshape ( results, ( results.shape[0], -1)) #ensure that dimensions are kept

                if perSetJob.extractionSlice:
                    results = results[:, perSetJob.extractionSlice]

                incrementVariableResultsArrays[ensElType] = results 
                variableLength = results.shape[1]
            
            partsDict[elSet.ensightPartID] = incrementVariableResultsArrays 
            
        if partsDict or  exportJob.writeEmptyTimeSteps :
            variableDimension = exportJob.dimensions or variableLength
            return es.EnsightPerElementVariable(exportJob.exportName, variableDimension, partsDict, )
        else:
            return None
        
        
    def outputDefinition(self, recordContent):
        """
        Attributes:  	
        1  –  Flag for element-based output (0), nodal output (1), modal output (2), or element set energy output (3).
        2  –  Set name (node or element set) used in the request (A8 format). This attribute is blank if no set was specified.
        3  –  Element type (only for element output, A8 format)."""
        
        flag = filFlag(recordContent[0])
        if flag == 0:
            setName = filStrippedString(recordContent[1])
            elType = filStrippedString(recordContent[2])
            self.currentEnsightElementType = self.ensightElementTypeMappings[elType]
                        
        elif flag == 1:
            setName = filStrippedString(recordContent[1])
            
        if not setName:
            setName = 'ALL'
        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]
            
        self.currentSetName = setName
        
    def elementHeaderRecord(self, rec):
        elNum = filFlag(rec[0]) 
        self.currentElementNum = elNum 
        self.currentIpt = filFlag(rec[1])
           
    def handlePerElementOutput(self, rec, location):
        res = filDouble(rec)
        currentIncrement = self.currentIncrement    
        currentSetName = self.currentSetName
        currentEnsightElementType =self.currentEnsightElementType
        currentElementNum = self.currentElementNum
        
        location += '@'+str(self.currentIpt)

        if location not in currentIncrement['elementResults']:
            currentIncrement['elementResults'][location] = {}

        if currentSetName not in currentIncrement['elementResults'][location]:
            currentIncrement['elementResults'][location][currentSetName] = {}

        if currentEnsightElementType not in currentIncrement['elementResults'][location][currentSetName]:
            currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType] =  {}
        
        targetLocation = currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType]
        if currentElementNum not in targetLocation:
            targetLocation[currentElementNum] = res
        else:
            # continuation of an existing record
            targetLocation[currentElementNum] = np.concatenate( (targetLocation[currentElementNum], res))

    def handlePerNodeOutput(self, rec, location):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])

        if location not in self.currentIncrement['nodeResults']:
            self.currentIncrement['nodeResults'][location] = {}
        self.currentIncrement['nodeResults'][location][node] = vals
    
    def addNode(self, recordContent):
        label = filInt(recordContent[0])[0]   
        coords = filDouble(recordContent[1:4])

        #make coords 3D, always!
        if ( coords.shape[0] < 3 ):
            coords = np.pad(coords, (0, 3-coords.shape[0]), mode='constant')

        node = self.allNodes[label]
        node.label = label

        node.coords = coords
    
    def addElement(self, recordContent):
        elNum = filInt(recordContent[0])[0]
        elType = filStrippedString(recordContent[1])
        elabqNodes = filInt(recordContent[2:])
        
        self.allElements[elNum]= Element(elNum, self.ensightElementTypeMappings[elType], [ self.allNodes[n] for n in elabqNodes ])
    
    def addElset(self, recordContent):
        setName = filStrippedString(recordContent[0])
        if not setName:
            setName = 'ALL'
        elif setName in  self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]
        self.currentSetName = setName
        abqElements = filInt(recordContent[1:])
        
        self.elSets[setName] = ElSet(setName, [self.allElements[e] for e in abqElements])

    def contAddElset(self, recordContent):
        abqElements = filInt(recordContent)
        self.elSets[self.currentSetName].appendElements([self.allElements[e] for e in abqElements] )
    
    def addNodeset(self, recordContent):
        setName = filStrippedString(recordContent[0])
        
        if not setName:
            setName = 'ALL'
        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]
            
        self.currentSetName = setName
        abqNodes = filInt(recordContent[1:])
        
        self.nSets[setName] = NSet(setName, [self.allNodes[n] for n in abqNodes ])
    
    def contAddNodeset(self, recordContent):
        abqNodes = filInt(recordContent)
        self.nSets[self.currentSetName].appendNodes( [self.allNodes[n] for n in abqNodes ] )
    
    def addIncrement(self, recordContent):
        self.currentState = 'increment parsing'
        r = recordContent
        tTotal, tStep = filDouble(r[0:2])
        
        nStep, nInc = filInt(r[5:7])
        timeInc = filDouble(r[10])[0]
        currentIncrement = self.currentIncrement
        currentIncrement['tTotal'] = tTotal
        currentIncrement['nInc'] = nInc
        currentIncrement['tStep'] = tStep
        currentIncrement['nStep'] = nStep
        currentIncrement['timeInc'] = timeInc
        currentIncrement['elementResults'] = {} 
        currentIncrement['nodeResults'] = {} 
        
    def addLabelCrossReference(self, recordContent):
        r = recordContent
        intKey = filFlag( r[0] )
        label = filStrippedString ( r[1:] )
        self.labelCrossReferences[ str(intKey) ] = label
        
    def surfaceDefHeader(self, recordContent):
        self.currentState = 'surface definition'

