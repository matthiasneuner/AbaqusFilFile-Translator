#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:43:53 2016

@author: matthias
"""

import numpy as np
from collections import OrderedDict, defaultdict
from utils.OrderedDefaultDict import OrderedDefaultDict
import utils.ensightgoldformat as es

def filInt(word):
    # return word[0:4].view( '<4i')
    return word.view('<i8').ravel()
def filString(word):
    return word.view('a8')
def filStrippedString(word):
    return filString(word).tostring().decode('utf-8').strip()
def filDouble(word):
    return word.view('<d').ravel()
def filFlag(word):
    return word[0:4].view('<i')[0]

class ExportEngine:
    ensightElementTypeMappings = {
                                  'CPS4': 'quad4',
               }
                        
    def __init__(self, exportJobs, exportName):
        
        self.perElementJobs = exportJobs.get('*ensightPerElementVariable', [])
        for entry in self.perElementJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['dimensions'] = entry.get('dimensions', len(entry['data']))
        self.perNodeJobs = exportJobs.get('*ensightPerNodeVariable', [])
        for entry in self.perNodeJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['dimensions'] = entry.get('dimensions', len(entry['data']))
            
        for entry in exportJobs.get('*defineElementType', []):
            self.ensightElementTypeMappings[entry['element']] = entry['shape']
            
        self.abqNodes = OrderedDict()
        self.abqElements = {}
        self.abqElSets = {}
        self.abqElSetToEnsightPartMappings = {}
        
        self.currentState = 'model setup'
        self.currentIncrement = {}
        self.currentAbqElSet = None
        self.currentSetName = 'mainPart'
        self.currentIpt = 1
        self.nIncrements = 0
        
        self.ensightCase = es.EnsightChunkWiseCase('.', exportName)
        
        self.knownRecords = { 1 : ('Element header record', self.elementHeaderRecord),
            5: ('SDV output', lambda x: self.handlePerElementOutput(x, 'SDV', appendGaussPt=False)),
            11 :('S output', lambda x : self.handlePerElementOutput(x, 'S')),
            21 :('E output', lambda x : self.handlePerElementOutput(x, 'E')),
            101: ('U output', self.handleUOutput),
            1901 : ('node definition', self.addNode),
            1900: ('element definition', self.addElement),
            1911: ('output request definition', self.outputDefinition),
            1921: ('heading ????', lambda x : None),
            1931: ('node set definition', self.addabqNodeset),
            1933: ('element set definition', self.addElset),
            1934: ('element set definition cont.', self.contAddElset),
            2000: ('start increment', self.addIncrement),
            2001: ('end increment', self.finishAndParseIncrement),
}
    
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
            geometryTimesetNumber = None 
            geometry = self.createEnsightGeometryFromModel()
            self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)
            self.currentState = 'increment parsing'
            
        elif self.currentState == 'increment parsing':
            self.nIncrements +=1
            timeSetID = 1
            self.ensightCase.setCurrentTime(timeSetID, self.currentIncrement['tTotal'])
            print('{:<25}{:>15.5f}'.format('parsing increment tTotal:',self.currentIncrement['tTotal']))
            
            for entry in self.perNodeJobs:
                jobElSetPartName = entry['set']
                resultLocation = entry['source']
                resultIndices = entry['data']
                resultTypeLength = entry['dimensions'] 
                jobName = entry['exportName']
                enSightVar = self.createEnsightPerNodeVariableFromJob(jobElSetPartName, jobName, resultLocation, resultIndices, resultTypeLength)
                self.ensightCase.writeVariableTrendChunk(enSightVar, timeSetID)
                del enSightVar
                                    
            for entry in self.perElementJobs:
                jobElSetPartName = entry['set']
                resultLocation = entry['source']
                resultIndices = entry['data'][0]
                jobName = entry['exportName']
                nCount = entry.get('periodicalPattern', 1)
                nShift = entry.get('periodicalShift', 0)
                for i in range(nCount):
                    enSightVar = self.createEnsightPerElementVariableFromJob(jobElSetPartName, 
                            jobName + (str(i+1) if nCount > 1 else ''), 
                            resultLocation, resultIndices + i*nShift)
                    self.ensightCase.writeVariableTrendChunk(enSightVar, timeSetID)
                del enSightVar
                
            del self.currentIncrement
            self.currentIncrement = {}
        
    def finalize(self):
        self.ensightCase.finalize()
    
    def createEnsightGeometryFromModel(self):
        geometry = es.EnsightGeometry("geometry", "descGEO", "desc2GEO")
        allNodeLabelList = list(self.abqNodes.keys())
        allNodes = np.asarray( list(self.abqNodes.values()))
        if allNodes.shape[1] < 3:
            allNodes = np.hstack( (allNodes, np.zeros((allNodes.shape[0],3 -allNodes.shape[1])) ))       
        allElements = np.asarray(  list(self.abqElements.keys()))
        self.abqElSets['mainPart'] = allElements
        
        for ensPartNumber ,( elSetName,  elSet) in enumerate(self.abqElSets.items()):
            ensPartNumber+=1 # ensight is not able to begin with #0
            self.abqElSetToEnsightPartMappings[elSetName] = ensPartNumber
            elSetPart = es.EnsightUnstructuredPart(elSetName, ensPartNumber)
            elSetPart.nodes = allNodes#mainPart.nodes # np.asarray ( list(abqNodes.values()) )           # all abqNodes
            elSetNodeLabels = allNodeLabelList#ist( abqNodes.keys () )
            for ensightElType in self.ensightElementTypeMappings.values():
                elSetPart.elements[ensightElType] = []
            for elLabel in elSet:
                element = self.abqElements[elLabel]
                elType, elNodeLabels = element
                elNodeIndices = [ elSetNodeLabels.index(i) for i in elNodeLabels]
                elSetPart.elements[self.ensightElementTypeMappings[elType]].append( (elLabel, elNodeIndices)   )
            geometry.partList.append(elSetPart)
        return geometry
        
    def createEnsightPerNodeVariableFromJob(self, jobElSetPartName, jobName, resultLocation, resultIndices, resultTypeLength):
        varDict = self.currentIncrement['nodeResults'][resultLocation]
        nodalVarTable = np.asarray([nodeVals[resultIndices] for nodeNum, nodeVals in varDict.items() ])    
        partsDict = {self.abqElSetToEnsightPartMappings[jobElSetPartName] : ('coordinates', nodalVarTable)}
        enSightVar = es.EnsightPerNodeVariable(jobName, resultTypeLength, partsDict)
        return enSightVar
        
    def createEnsightPerElementVariableFromJob(self, jobElSetPartName, jobName, resultLocation, resultIndices):
        enSightVar = es.EnsightPerElementVariable(jobName, len(resultIndices),)
        varDict = self.currentIncrement['elementResults'][resultLocation][jobElSetPartName]

        varDict = { ensElType : np.asarray([ np.concatenate(chunks) for chunks in elResults.values() ])[:,resultIndices]
                                           for ensElType, elResults in varDict.items()}

        enSightVar.partsDict[self.abqElSetToEnsightPartMappings[jobElSetPartName]] = varDict
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
        self.currentSetName = setName
        
    def elementHeaderRecord(self, rec):
        elNum = filInt(rec[0])[0]
        self.currentElementNum = elNum #abqElements[elNum]
        self.currentIpt = filFlag(rec[1])
           
    def handlePerElementOutput(self, rec, location, appendGaussPt=True):
        res = filDouble(rec)
        currentIncrement = self.currentIncrement    
        currentSetName = self.currentSetName
        currentEnsightElementType =self.currentEnsightElementType
        currentElementNum = self.currentElementNum
        
        if appendGaussPt:
            location += '@'+str(self.currentIpt)

        insertLocation = currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType][currentElementNum]
        insertLocation.append(res)
    
    def handleUOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])
        self.currentIncrement['nodeResults']['U'][node] = vals
    
    def addNode(self, recordContent):
        node = filInt(recordContent[0])[0]   
        coords = filDouble(recordContent[1:])
        self.abqNodes[node] = coords
    
    def addElement(self, recordContent):
        elNum = filInt(recordContent[0])[0]
        elType = filStrippedString(recordContent[1])
        elabqNodes = filInt(recordContent[2:])
        self.abqElements[elNum] = (elType, elabqNodes)
    
    def addElset(self, recordContent):
        setName = filStrippedString(recordContent[0])
        if not setName:
            setName = 'mainPart'
        self.currentSetName = setName
        abqElements = filInt(recordContent[1:])
        self.abqElSets[setName] = abqElements

    def contAddElset(self, recordContent):
        abqElements = filInt(recordContent)
        self.abqElSets[self.currentSetName] = np.concatenate( [self.abqElSets[self.currentSetName], abqElements])
    
    def addabqNodeset(self, recordContent):
        pass
    
    def addIncrement(self, recordContent):
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
        currentIncrement['elementResults'] = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: OrderedDefaultDict(list))))
        currentIncrement['nodeResults'] = defaultdict(OrderedDict)

