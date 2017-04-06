#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:43:53 2016

@author: matthias
"""

import numpy as np
from collections import OrderedDict, defaultdict
import utils.ensightgoldformat as es

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

class ExportEngine:
    ensightElementTypeMappings = {
                                  'CPS4': 'quad4',
               }
                        
    def __init__(self, exportJobs, exportName):
        
        # TODO : CHECK 'INDICES' OF perNodeJobs and perElementJobs
        
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

        self.csvPerElementSetJobs = exportJobs.get('*csvPerElementOutput', [])
        for entry in self.csvPerElementSetJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['csvData'] = []
        
        self.csvPerNodeJobs = exportJobs.get('*csvPerNodeOutput', [])
        for entry in self.csvPerNodeJobs:
            entry['exportName'] = entry.get('exportName', entry['source']) 
            entry['csvData'] = []
#            if entry['math']:
#                entry['math'] = mathFunctions.get(entry['math'], False)
            
        self.abqNodes = OrderedDict()
        self.abqElements = {}
        self.abqNodeSets = {}
        self.abqElSets = {}
        self.abqElSetToEnsightPartMappings = {}
        
        self.currentState = 'model setup'
        self.currentIncrement = {}
        self.currentAbqElSet = None
        self.currentSetName = 'mainPart'
        self.currentIpt = 1
        self.nIncrements = 0
        self.timeHistory = []

        if exportJobs['*exportTimeHistory']:
            self.exportTimeHistory = exportJobs['*exportTimeHistory'][0].get('exportName', 'timeHistory')
        else:
            self.exportTimeHistory = False
        
        self.ensightCase = es.EnsightChunkWiseCase('.', exportName)
        self.ensightCaseDiscardTimeMarks = exportJobs.get('discardTime', False)
        
        self.knownRecords = { 
            1: ('Element header record', self.elementHeaderRecord),
            5: ('SDV output', lambda x: self.handlePerElementOutput(x, 'SDV', appendGaussPt=False)),
            11: ('S output', lambda x : self.handlePerElementOutput(x, 'S')),
            21: ('E output', lambda x : self.handlePerElementOutput(x, 'E')),
            101: ('U output', self.handleUOutput),
            104: ('RF output', self.handleRFOutput),
            201: ('NT output', self.handleNTOutput),
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
            1940: ('label cross reference', lambda x : None),
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
                
            for entry in self.csvPerElementSetJobs:
                currentRow = []
                jobElSetPartName = entry['set']
                resultLocation = entry['source']
                resultIndices = entry['data'][0]
                nCount = entry.get('periodicalPattern', 1)
                nShift = entry.get('periodicalShift', 0)
                for elDict in self.currentIncrement['elementResults'][resultLocation][jobElSetPartName].values():
                    for elResult in elDict.values():
                        for i in range(nCount):
                            currentRow.append( elResult[resultIndices + i*nShift].ravel() )
                entry['csvData'].append(currentRow)
                            
            for entry in self.csvPerNodeJobs:
                if 'nodes' not in entry:
                    # first run (=first increment)
                    if 'node' in entry:
                        entry['nodes'] = [entry['node']]
                    elif 'nSet' in entry:
                        entry['nodes'] = self.abqNodeSets[entry['nSet']]
                        
                resultLocation = entry['source']
                resultIndices = entry['data'][0]
                currentRow = [self.currentIncrement['nodeResults'][resultLocation][node][resultIndices] for node in entry['nodes']]
#                if entry['math']:
#                    currentRow = entry['math'](currentRow)
                entry['csvData'].append(currentRow)
                
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
        allNodeLabelList = list(self.abqNodes.keys())
        allNodes = np.asarray( list(self.abqNodes.values()))
        if allNodes.shape[1] < 3:
            allNodes = np.hstack( (allNodes, np.zeros((allNodes.shape[0],3 -allNodes.shape[1])) ))       
        allElements = np.asarray(  list(self.abqElements.keys()))
        self.abqElSets['mainPart'] = allElements
        
        partList = []
        for ensPartNumber ,( elSetName,  elSet) in enumerate(self.abqElSets.items()):
            ensPartNumber+=1 # ensight is not able to begin with #0
            self.abqElSetToEnsightPartMappings[elSetName] = ensPartNumber
            elSetNodeLabels = allNodeLabelList#ist( abqNodes.keys () )
            setElements = { ensightElType : [] for ensightElType in self.ensightElementTypeMappings.values()}
            for elLabel in elSet:
                element = self.abqElements[elLabel]
                elType, elNodeLabels = element
                elNodeIndices = [ elSetNodeLabels.index(i) for i in elNodeLabels]
                setElements[self.ensightElementTypeMappings[elType]].append( (elLabel, elNodeIndices)   )
            elSetPart = es.EnsightUnstructuredPart(elSetName, ensPartNumber, setElements, allNodes, elSetNodeLabels)
            partList.append(elSetPart)
        
        geometry = es.EnsightGeometry('geometry', '-', '-', partList, 'given', 'given')
        return geometry
        
    def createEnsightPerNodeVariableFromJob(self, jobElSetPartName, jobName, resultLocation, resultIndices, resultTypeLength):
        varDict = self.currentIncrement['nodeResults'][resultLocation]
        nodalVarTable = np.asarray([nodeVals[resultIndices] for nodeNum, nodeVals in varDict.items() ])    
        partsDict = {self.abqElSetToEnsightPartMappings[jobElSetPartName] : ('coordinates', nodalVarTable)}
        enSightVar = es.EnsightPerNodeVariable(jobName, resultTypeLength, partsDict)
        return enSightVar
        
    def createEnsightPerElementVariableFromJob(self, jobElSetPartName, jobName, resultLocation, resultIndices):
        enSightVar = es.EnsightPerElementVariable(jobName, len(resultIndices),)
        incrementVariableResults = self.currentIncrement['elementResults'][resultLocation][jobElSetPartName]
        
        incrementVariableResultsArrays = {}
        for ensElType, elDict in incrementVariableResults.items():
            incrementVariableResultsArrays[ensElType] = np.asarray([row[resultIndices] for row in  elDict.values() ])

        enSightVar.partsDict[self.abqElSetToEnsightPartMappings[jobElSetPartName]] = incrementVariableResultsArrays 
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
#        elNum = filInt(rec[0])[0]
        elNum = filFlag(rec[0]) # march 2017: filInt does not worki in tensilebar?
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
        
        if currentEnsightElementType not in currentIncrement['elementResults'][location][currentSetName]:
            newDict = OrderedDict(zip(self.abqElSets[currentSetName], [None] * len(self.abqElSets[currentSetName]) ))
            currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType] =  newDict
        
        targetLocation = currentIncrement['elementResults'][location][currentSetName][currentEnsightElementType]
        if targetLocation[currentElementNum] is None:
            targetLocation[currentElementNum] = res
        else:
            targetLocation[currentElementNum] = np.concatenate( (targetLocation[currentElementNum], res))
    
    def handleUOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])
        self.currentIncrement['nodeResults']['U'][node] = vals

    def handleRFOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])
        self.currentIncrement['nodeResults']['RF'][node] = vals

    def handleNTOutput(self, rec):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])
        self.currentIncrement['nodeResults']['NT'][node] = vals
    
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
    
    def addNodeset(self, recordContent):
        setName = filStrippedString(recordContent[0])
        if not setName:
            setName = 'defaultNodeSet'
        self.currentSetName = setName
        abqNodes = filInt(recordContent[1:])
        self.abqNodeSets[setName] = abqNodes
    
    def contAddNodeset(self, recordContent):
        abqNodes = filInt(recordContent)
        self.abqNodeSets[self.currentSetName] = np.concatenate( [self.abqNodeSets[self.currentSetName], abqNodes])
    
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
        currentIncrement['elementResults'] = defaultdict(lambda: defaultdict(lambda: dict() )) 
        currentIncrement['nodeResults'] = defaultdict(OrderedDict)
        self.timeHistory.append(tTotal)

