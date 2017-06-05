#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:43:53 2016

@author: matthias
"""


"""

TODO: Add Input parameter for dimension ( necessary as abaqus provides no info about the used dim!)


"""

import numpy as np
from collections import OrderedDict, defaultdict
import utils.ensightgoldformat as es

from utils.OrderedDefaultDict import OrderedDefaultDict

mathFunctions = {'sum' : np.sum,
                 'mean': np.mean
                 }

def filInt(word):
#     return word[0:4].view( '<4i')
    return word.view('<i8').ravel()
def filString(word):
    return word.view('a8')
def filStrippedString(word):
    return filString(word).tostring().decode('utf-8').strip()
def filDouble(word):
    return word.view('<d').ravel()
def filFlag(word):
    return word[0:4].view('<i')[0]


def sliceFromString(string, shift=0):
    if ':' in string:
        a, b = string.split(':'),
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
        
        self.elements = OrderedDefaultDict(list)
        self.appendElements(elements)
            
        self._elementLabels = None
        self._reducedNodes = None
        self._reducedElements = None
        self._reducedNodeCoords3D = None
        
    def appendElements(self, elements):
        for element in elements:
            self.elements[element.shape].append(element)
            
    def getEnsightCompatibleReducedNodes(self, ):
        if self._reducedNodes is None:
            self._reducedNodes = OrderedDict([ (node.label, node) for elementsByShape in self.elements.values() 
                                                                for element in elementsByShape 
                                                                    for node in element.nodes ])
        return self._reducedNodes
    
    def getEnsightCompatibleReducedNodeCoords(self,):
        if self._reducedNodeCoords3D is None:
            self._reducedNodeCoords3D = np.asarray([node.coords for node in self.getEnsightCompatibleReducedNodes().values()])
            if self._reducedNodeCoords3D.shape[1] < 3:
                self._reducedNodeCoords3D = np.hstack( (self._reducedNodeCoords3D, np.zeros((self._reducedNodeCoords3D.shape[0],3 -self._reducedNodeCoords3D.shape[1])) ))   
        return self._reducedNodeCoords3D
    
    def getEnsightCompatibleElementNodeIndices(self,):
        if self._reducedElements is None:
            self._reducedNodeIndices = {node : i for (i, node) in enumerate(self.getEnsightCompatibleReducedNodes().keys()) }
            self._reducedElements = OrderedDict()
            for eShape, elements in self.elements.items():
                self._reducedElements[eShape] = [ (e.label , [ self._reducedNodeIndices[n.label] for n in e.nodes] ) for e in elements  ]
        return self._reducedElements
    
class NSet:
    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes
        
    def appendNodes(self, nodes):
        self.nodes += nodes

class ExportEngine:
    ensightElementTypeMappings = {
                                  'CPS4': 'quad4',
                                  'B21': 'bar2',
               }
                        
    def __init__(self, exportJobs, exportName):
        self.perElementJobs = exportJobs.get('*ensightPerElementVariable', [])
        for entry in self.perElementJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
#            entry['dimensions'] = entry.get('dimensions', len(entry['data']))
            entry['periodicalJobs'] = []
#            entry['dimension'] = 
            for i in range(entry.get('periodicalPattern', 1)):
                entry['periodicalJobs'].append ( sliceFromString( entry['data'][0][0], shift=entry.get('periodicalShift', 0) * i  ) )
            
        self.perNodeJobs = exportJobs.get('*ensightPerNodeVariable', [])
        for entry in self.perNodeJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['dimensions'] = entry.get('dimensions', len(entry['data']))
            entry['slice'] = sliceFromString (entry['data'][0][0])
            
        for entry in exportJobs.get('*defineElementType', []):
            self.ensightElementTypeMappings[entry['element']] = entry['shape']

        self.csvPerElementSetJobs = exportJobs.get('*csvPerElementOutput', [])
        for entry in self.csvPerElementSetJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['csvData'] = []
        
        self.csvPerNodeJobs = exportJobs.get('*csvPerNodeOutput', [])
        for entry in self.csvPerNodeJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['csvData'] = []
            entry['slice'] = sliceFromString (entry['data'][0][0])
            if entry.get('math', False):
                entry['math'] = mathFunctions.get(entry['math'], False)
            
        self.allNodes = OrderedDefaultDict(Node)
        
        self.allElements = OrderedDict()
        self.nSets = {}
        self.elSets = {}
        
        self.currentState = 'model setup'
        self.currentIncrement = {}
        self.currentAbqElSet = None
        self.currentSetName = 'mainPart'
        self.currentIpt = 1
        self.nIncrements = 0
        self.timeHistory = []
        self.labelCrossReferences = {}

        if exportJobs['*exportTimeHistory']:
            self.exportTimeHistory = exportJobs['*exportTimeHistory'][0].get('exportName', 'timeHistory')
        else:
            self.exportTimeHistory = False
        
        self.ensightCase = es.EnsightChunkWiseCase('.', exportName)
        self.ensightCaseDiscardTimeMarks = False# exportJobs.get('discardTime', False)
        
        self.knownRecords = { 
            1: ('Element header record', self.elementHeaderRecord),
            5: ('SDV output', lambda x: self.handlePerElementOutput(x, 'SDV', appendGaussPt=False)),
            11: ('S output', lambda x : self.handlePerElementOutput(x, 'S')),
            21: ('E output', lambda x : self.handlePerElementOutput(x, 'E')),
            101: ('U output', self.handleUOutput),
            104: ('RF output', self.handleRFOutput),
            201: ('NT output', self.handleNTOutput),
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
            timeSetID = 1
            self.ensightCase.setCurrentTime(timeSetID, self.currentIncrement['tTotal'])
            print('parsing increment {:>5}  tTotal:{:>16.5f}'.format(self.nIncrements,self.currentIncrement['tTotal']))
            
            for entry in self.perNodeJobs:
                jobElSetPartName = entry['set']
                resultLocation = entry['source']
#                if not 'slice' in entry:
#                    entry['slice'] = sliceFromDataline(entry['data'][0])
#                resultIndices = entry['data'][0]
                resultIndices = entry['slice']
                resultTypeLength = entry['dimensions'] 
                jobName = entry['exportName']
                enSightVar = self.createEnsightPerNodeVariableFromJob(jobElSetPartName, jobName, resultLocation, resultIndices, resultTypeLength)
                self.ensightCase.writeVariableTrendChunk(enSightVar, timeSetID)
                del enSightVar
                                    
            for entry in self.perElementJobs:
                jobElSetPartName = entry['set']
                resultLocation = entry['source']
#                resultIndices = entry['data'][0]
#                if not 'slice' in entry:
#                    entry['slices'] = sliceFromDataline(entry['data'][0])
                jobName = entry['exportName']
#                dimension = entry['dimensions']
                for i, periodicalJobSlice in enumerate(entry['periodicalJobs']):
                    
#                nCount = entry.get('periodicalPattern', 1)
#                nShift = entry.get('periodicalShift', 0)
#                for i in range(nCount):
                    enSightVar = self.createEnsightPerElementVariableFromJob(jobElSetPartName, 
                            jobName + (str(i+1) if len (entry['periodicalJobs']) > 1 else ''), 
                            resultLocation, periodicalJobSlice)
                    self.ensightCase.writeVariableTrendChunk(enSightVar, timeSetID)
                del enSightVar
                
#            for entry in self.csvPerElementSetJobs:
#                currentRow = []
#                jobElSetPartName = entry['set']
#                resultLocation = entry['source']
##                resultIndices = entry['data'][0]
##                nCount = entry.get('periodicalPattern', 1)
##                nShift = entry.get('periodicalShift', 0)
#                for elDict in self.currentIncrement['elementResults'][resultLocation][jobElSetPartName].values():
##                    
##                    for elResult in elDict.values():
##                        for i in range(nCount):
#                            currentRow.append( elResult[resultIndices + i*nShift].ravel() )
#                entry['csvData'].append(currentRow)
                            
            for entry in self.csvPerNodeJobs:
                if 'nodes' not in entry:
                    # first run (=first increment)
                    if 'node' in entry:
                        entry['nodes'] = self.allNodes[entry['node']]
                    elif 'nSet' in entry:
                        entry['nodes'] = self.nSets[entry['nSet']].nodes
                        
                resultLocation = entry['source']
                resultIndices = entry['slice']
                currentRow = [self.currentIncrement['nodeResults'][resultLocation][node.label][resultIndices] for node in entry['nodes']]
                if entry.get('math', False):
                    currentRow = entry['math'](currentRow)
                entry['csvData'].append(currentRow)
            
            if self.nIncrements % 10 == 0:
                # intermediate saving ...
                self.ensightCase.finalize(discardTimeMarks = self.ensightCaseDiscardTimeMarks, closeFileHandles=False)
                
            del self.currentIncrement
            self.currentIncrement = {}
        
    def finalize(self):
        self.ensightCase.finalize(discardTimeMarks = self.ensightCaseDiscardTimeMarks)
        
        for entry in self.csvPerElementSetJobs:
            jobName = entry['exportName']
            table = np.asarray(entry['csvData'] )
            np.savetxt('{:}.csv'.format(jobName), table, fmt=entry.get('fmt', '%.6e'), )
                
        for entry in self.csvPerNodeJobs:
            jobName = entry['exportName']
            table = np.asarray(entry['csvData']  )
            np.savetxt('{:}.csv'.format(jobName), table, fmt=entry.get('fmt', '%.6e'), )
                
        if self.exportTimeHistory:
            completeTimeHistory = np.asarray(self.timeHistory)
            np.savetxt('{:}.csv'.format(self.exportTimeHistory), completeTimeHistory, fmt='%.6e', )
    
    def createEnsightGeometryFromModel(self):
        
        partList = []
        partNumber = 1
        
        mainPartSet = ElSet('mainPart', self.allElements.values())
        self.elSets['mainPart'] = mainPartSet
        
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
    
    def createEnsightPerNodeVariableFromJob(self, jobElSetPartName, jobName, resultLocation, resultIndices, resultTypeLength):
        elSet = self.elSets[jobElSetPartName]

        nodalVarTable = np.asarray([ self.currentIncrement['nodeResults'][resultLocation][node] for node in elSet.getEnsightCompatibleReducedNodes().keys() ] )
        
        partsDict = {elSet.ensightPartID : ('coordinates', nodalVarTable)}
        enSightVar = es.EnsightPerNodeVariable(jobName, resultTypeLength, partsDict)
        
        return enSightVar
        
    def createEnsightPerElementVariableFromJob(self, jobElSetPartName, jobName, resultLocation, resultIndices):
        
        elSet = self.elSets[jobElSetPartName]
        incrementVariableResults = self.currentIncrement['elementResults'][resultLocation][jobElSetPartName]
        
        incrementVariableResultsArrays = {}
        for ensElType, elDict in incrementVariableResults.items():
            incrementVariableResultsArrays[ensElType] = np.asarray([ elDict[el.label][resultIndices] for el in elSet.elements[ensElType] ])

        
        varDim = list(incrementVariableResultsArrays.values())[0].shape[1]
#        print(type(varDim))
        enSightVar = es.EnsightPerElementVariable(jobName, varDim, None, )
        enSightVar.partsDict[elSet.ensightPartID] = incrementVariableResultsArrays 
        return enSightVar
        
    def outputDefinition(self, recordContent):
        """
        Attributes:  	
        1  –  Flag for element-based output (0), nodal output (1), modal output (2), or element set energy output (3).
        2  –  Set name (node or element set) used in the request (A8 format). This attribute is blank if no set was specified.
        3  –  Element type (only for element output, A8 format)."""
        
        flag = filFlag(recordContent[0])
        if flag == 0:
            setName = filStrippedString(recordContent[1])
            self.currentEnsightElementType = self.ensightElementTypeMappings[filStrippedString(recordContent[2])]
            
        elif flag == 1:
            setName = filStrippedString(recordContent[1])
            
        if not setName:
            setName = 'mainPart'
        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]
            
        self.currentSetName = setName
        
    def elementHeaderRecord(self, rec):
#        elNum = filInt(rec[0])[0]
        elNum = filFlag(rec[0]) # march 2017: filInt does not worki in tensilebar?
        self.currentElementNum = elNum 
        self.currentIpt = filFlag(rec[1])
           
    def handlePerElementOutput(self, rec, location, appendGaussPt=True):
        res = filDouble(rec)
        currentIncrement = self.currentIncrement    
        currentSetName = self.currentSetName
        currentEnsightElementType =self.currentEnsightElementType
        currentElementNum = self.currentElementNum
        
        if appendGaussPt:
            location += '@'+str(self.currentIpt)
        
        if currentEnsightElementType not in currentIncrement['elementResults'][location][currentSetName]:
            currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType] =  {}
        
        targetLocation = currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType]
        if currentElementNum not in targetLocation:
            targetLocation[currentElementNum] = res
        else:
            targetLocation[currentElementNum] = np.concatenate( (targetLocation[currentElementNum], res))
    
    def handleUOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:3])
        self.currentIncrement['nodeResults']['U'][node] = vals

    def handleRFOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])
        self.currentIncrement['nodeResults']['RF'][node] = vals

    def handleNTOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1])
        self.currentIncrement['nodeResults']['NT'][node] = vals
    
    def addNode(self, recordContent):
        label = filInt(recordContent[0])[0]   
        coords = filDouble(recordContent[1:4])
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
            setName = 'mainPart'
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
            setName = 'defaultNodeSet'
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
        currentIncrement['elementResults'] = defaultdict(lambda: defaultdict(lambda: dict() )) 
        currentIncrement['nodeResults'] = defaultdict(OrderedDict)
        self.timeHistory.append(tTotal)
        
    def addLabelCrossReference(self, recordContent):
        r = recordContent
        intKey = filFlag( r[0] )
        label = filStrippedString ( r[1:] )
        self.labelCrossReferences[ str(intKey) ] = label
        
    def surfaceDefHeader(self, recordContent):
        self.currentState = 'surface definition'

