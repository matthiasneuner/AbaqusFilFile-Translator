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
    return word.view("<i8").ravel()


def filString(word):
    return word.view("a8")


def filStrippedString(word):
    return filString(word).tostring().decode("utf-8").strip()


def filDouble(word):
    return word.view("<d").ravel()


def filFlag(word):
    return word[0:4].view("<i")[0]


def makeExtractionFunction(expression, symbol="x"):
    """make a simple f(x) expression from string"""
    return lambda x: eval(expression, globals(), {symbol: x})


def sliceFromString(string, shift=0):
    """generate a slice from a string, which can represent a slice or an index"""
    if ":" in string:
        a, b = string.split(":")
        return slice(int(a) + shift, int(b) + shift)
    else:
        return slice(int(string) + shift, int(string) + 1 + shift)


def RecursiveDefaultDict():
    return defaultdict(RecursiveDefaultDict)
    # return dict()


class ElementDefinition:
    def __init__(self, label, elementType, nodeLabels):
        self.label = label
        self.shape = elementType
        self.nodeLabels = nodeLabels


class ElSetDefinition:
    def __init__(self, name, elementLabels):
        self.name = name
        self.elementLabels = elementLabels
        self.appendElementLabels(elementLabels)

    def appendElementLabels(self, elementLabels):
        self.elementLabels += elementLabels


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
    def __init__(self, name, elements, ensightPartID=None):
        self.name = name
        self.ensightPartID = ensightPartID

        self.elementsByShape = defaultdict(list)

        for element in elements:
            self.elementsByShape[element.shape].append(element)

        self.reducedNodes = self._getEnsightCompatibleReducedNodes()
        self.reducedNodeIndices = self._getEnsightCompatibleElementNodeIndices()
        self.reducedElements = self._getEnsightCompatibleElements()
        self.reducedNodeCoords3D = self._getEnsightCompatibleReducedNodeCoords()

    def _getEnsightCompatibleReducedNodes(
        self,
    ):
        reducedNodes = dict(
            [
                (node.label, node)  # (node, self.allNodes[node])
                for elementsByShape in self.elementsByShape.values()
                for element in elementsByShape
                for node in element.nodes
            ]
        )

        return reducedNodes

    def _getEnsightCompatibleReducedNodeCoords(
        self,
    ):
        reducedNodeCoords3D = np.asarray(
            [node.coords for node in self.reducedNodes.values()]
        )

        return reducedNodeCoords3D

    def _getEnsightCompatibleElementNodeIndices(
        self,
    ):
        reducedNodeIndices = {
            node: i for (i, node) in enumerate(self.reducedNodes.keys())
        }
        return reducedNodeIndices

    def _getEnsightCompatibleElements(
        self,
    ):
        reducedElements = dict()

        for eShape, elements in self.elementsByShape.items():
            reducedElements[eShape] = [
                (e.label, [self.reducedNodeIndices[n.label] for n in e.nodes])
                for e in elements
            ]
        return reducedElements


class NSetDefinition:
    def __init__(self, name, nodeLabels):
        self.name = name
        self.nodeLabels = nodeLabels

    def appendNodeLabels(self, nodeLabels):
        self.nodeLabels += nodeLabels


class NSet:
    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes


class EnsightExportJob:
    def __init__(self, name, dimensions, timeSetID, writeEmptyTimeSteps=True):
        self.exportName = name
        self.dimensions = dimensions
        self.timeSetID = timeSetID
        self.entries = {}
        self.writeEmptyTimeSteps = writeEmptyTimeSteps


class EnsightPerSetJobEntry:
    def __init__(
        self,
        job,
        setName,
        result,
        location,
        which,
        extractionSlice=None,
        extractionFunction=None,
        offset=None,
        fillMissingValuesTo=None,
    ):
        self.job = job
        self.setName = setName
        self.extractionSlice = extractionSlice
        self.extractionFunction = extractionFunction
        self.offset = offset
        self.result = result
        self.location = location
        self.which = which
        self.fillMissingValuesTo = (
            fillMissingValuesTo
            if fillMissingValuesTo is None
            else float(fillMissingValuesTo)
        )


class ExportEngine:
    def __init__(self, exportJobs, caseName):
        self.uelSdvToQpJobs = self.collectUelSDVToQpJobs(
            exportJobs["*UELSDVToQuadraturePoints"]
        )
        self.qpAverageJobs = self.collectQpAverageJobs(
            exportJobs["*computeAverageOverQuadraturePoints"]
        )

        self.perElementJobs = self.collectExportJobs(
            exportJobs["*ensightPerElementVariableJob"]
        )
        self.perNodeJobs = self.collectExportJobs(
            exportJobs["*ensightPerNodeVariableJob"]
        )

        self.perElementJobs = self.collectPerElementJobEntries(
            exportJobs["*ensightPerElementVariableJobEntry"],
            self.perElementJobs,
        )
        self.perNodeJobs = self.collectPerNodeJobEntries(
            exportJobs["*ensightPerNodeVariableJobEntry"],
            self.perNodeJobs,
        )

        self.ensightElementTypeMappings = {
            x["element"]: x["shape"] for x in exportJobs["*defineElementType"]
        }
        self.ignoreLastNodesForElType = {
            x["element"]: x["number"]
            for x in exportJobs["*ignoreLastNodesForElementType"]
        }

        self.allNodes = {}
        # add a default node, to which abaqus falls back if it creates node in place (e.g, for hex27 elements in contact)
        self.allNodes[0] = Node(0, np.array([0.0, 0.0, 0.0]))

        self.elementDefinitions = {}
        self.elSetDefinitions = {}
        self.nSetDefinitions = {}

        self.allElements = {}
        self.nSets = {}
        self.elSets = {}

        self.currentState = "model setup"
        self.currentIncrement = {}
        self.currentAbqElSet = None
        self.currentSetName = "ALL"
        self.currentIpt = 1
        self.nIncrements = 0
        self.timeHistory = []
        self.labelCrossReferences = {}
        self.ensightCase = es.EnsightChunkWiseCase(".", caseName)
        self.ensightCaseDiscardTimeMarks = False

        self.knownRecords = {
            1: ("Element header record", self.elementHeaderRecord),
            5: ("SDV output", lambda x: self.handlePerElementOutput(x, "SDV")),
            11: ("S output", lambda x: self.handlePerElementOutput(x, "S")),
            21: ("E output", lambda x: self.handlePerElementOutput(x, "E")),
            22: ("PE output", lambda x: self.handlePerElementOutput(x, "PE")),
            101: ("U output", lambda x: self.handlePerNodeOutput(x, "U")),
            102: ("V output", lambda x: self.handlePerNodeOutput(x, "V")),
            103: ("A output", lambda x: self.handlePerNodeOutput(x, "A")),
            104: ("RF output", lambda x: self.handlePerNodeOutput(x, "RF")),
            201: ("NT output", lambda x: self.handlePerNodeOutput(x, "NT")),
            1501: ("Surface definition header", self.surfaceDefHeader),
            1502: ("Surface facet", lambda x: None),
            1900: ("element definition", self.addElement),
            1901: ("node definition", self.addNode),
            1902: ("active dof", lambda x: None),
            1911: ("output request definition", self.outputDefinition),
            1921: ("heading", lambda x: None),
            1922: ("heading", lambda x: None),
            1931: ("node set definition", self.createNodeSetDefinition),
            1932: ("node set definition cont.", self.contNodeSetDefinition),
            1933: ("element set definition", self.addElset),
            1934: ("element set definition cont.", self.contAddElset),
            1940: ("label cross reference", self.addLabelCrossReference),
            2000: ("start increment", self.addIncrement),
            2001: ("end increment", self.finishAndParseIncrement),
        }

    def computeRecord(self, recordLength, recordType, recordContent):
        if recordType in self.knownRecords:
            doc, action = self.knownRecords[recordType]
            action(recordContent)
            return True
        else:
            print(
                "{:<20}{:>6}{:>10}{:>4}".format(
                    "unknown record:", recordType, " of length", recordLength
                )
            )
            return False

    def finishAndParseIncrement(self, recordContent):
        if self.currentState == "model setup":
            ALLSet = ElSetDefinition("ALL", list(self.elementDefinitions.keys()))
            self.elSetDefinitions["ALL"] = ALLSet

            # Time to create Elements, ElSets and NodeSets from the definitions!
            for elDef in self.elementDefinitions.values():
                self.allElements[elDef.label] = Element(
                    elDef.label,
                    elDef.shape,
                    [self.allNodes[label] for label in elDef.nodeLabels],
                )
                # el.nodes = [self.allNodes[label] for label in el.nodes]

            for elSetDef in self.elSetDefinitions.values():
                self.elSets[elSetDef.name] = ElSet(
                    elSetDef.name,
                    [self.allElements[label] for label in elSetDef.elementLabels],
                )

            for nSetDef in self.nSetDefinitions.values():
                self.nSets[nSetDef.name] = NSet(
                    nSetDef.name, [self.allNodes[n] for n in nSetDef.nodeLabels]
                )

            # replace all label references by the respective labels
            for key, label in self.labelCrossReferences.items():
                strKey = key  # str(key)
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

        elif self.currentState == "surface definition":
            pass

        elif self.currentState == "increment parsing":
            self.nIncrements += 1
            self.ensightCase.setCurrentTime(self.currentIncrement["tTotal"])
            self.timeHistory.append(self.currentIncrement["tTotal"])

            print("*" * 80)
            print("increment contained element results for")
            print(
                "\n".join(
                    [
                        " {:5} [{:}]".format(
                            resName, ", ".join([s for s in resEntries])
                        )
                        for resName, resEntries in self.currentIncrement[
                            "elementResults"
                        ].items()
                    ]
                )
            )
            print("")
            print("increment contained node results for")
            print(
                "\n".join(
                    [
                        " {:5} [{:10} nodes]".format(resName, len(resEntries))
                        for resName, resEntries in self.currentIncrement[
                            "nodeResults"
                        ].items()
                    ]
                )
            )
            print("")

            print("exporting...")

            for exportJob in self.perNodeJobs.values():
                enSightVar = self.createEnsightPerNodeVariableFromPerNodeJob(exportJob)
                if enSightVar:
                    self.ensightCase.writeVariableTrendChunk(
                        enSightVar, exportJob.timeSetID
                    )
                    del enSightVar

            # operate on elemetal results (e.g. compute average over quadraturePoint )
            for uelSdvToQpJob in self.uelSdvToQpJobs:
                self.computeUelSdvToQp(uelSdvToQpJob)

            for qpAverageJob in self.qpAverageJobs:
                self.computeQpAverage(qpAverageJob)

            for exportJob in self.perElementJobs.values():
                enSightVar = self.createEnsightPerElementVariableFromPerElementJob(
                    exportJob
                )
                if enSightVar:
                    self.ensightCase.writeVariableTrendChunk(
                        enSightVar, exportJob.timeSetID
                    )
                    del enSightVar

            if self.nIncrements % 10 == 0:
                # intermediate saving ...
                self.ensightCase.finalize(
                    discardTimeMarks=self.ensightCaseDiscardTimeMarks,
                    closeFileHandles=False,
                )

            # data might consume a lot of memory, so we delete it (explicitly)
            del self.currentIncrement
            # and create a new dict for a new increment; it is filled with the next start increment entry in the fil file
            self.currentIncrement = dict()

    def finalize(self):
        self.ensightCase.finalize(discardTimeMarks=self.ensightCaseDiscardTimeMarks)

    def createEnsightGeometryFromModel(self):
        partList = []
        partNumber = 1

        for elSet in self.elSets.values():
            elSetPart = es.EnsightUnstructuredPart(
                elSet.name,
                partNumber,
                elSet._getEnsightCompatibleElements(),
                elSet._getEnsightCompatibleReducedNodeCoords(),
                list(elSet._getEnsightCompatibleReducedNodes().keys()),
            )
            elSet.ensightPartID = partNumber
            partList.append(elSetPart)
            partNumber += 1

        geometry = es.EnsightGeometry("geometry", "-", "-", partList, "given", "given")
        return geometry

    def collectQpAverageJobs(self, entries):
        jobs = []
        for entry in entries:
            jobs.append(entry)
        return jobs

    def computeQpAverage(self, job):
        result = job["result"]
        setName = job["set"]

        setResults = self.currentIncrement["elementResults"][result][setName]

        for elTypeResults in setResults.values():
            for elResults in elTypeResults.values():
                elResults["computed"]["average"] = np.mean(
                    [qpRes for qpRes in elResults["qps"].values()], axis=0
                )

    def collectUelSDVToQpJobs(self, entries):
        jobs = []

        for entry in entries:
            offset = entry["qpInitialOffset"]
            nQps = entry["qpCount"]
            qpDistance = entry["qpDistance"]
            entry["qpSlices"] = [
                slice(offset + i * qpDistance, offset + (i + 1) * qpDistance)
                for i in range(nQps)
            ]

            jobs.append(entry)

        return jobs

    def computeUelSdvToQp(self, job):
        setName = job["set"]
        destination = job["destination"]
        qpSlices = job["qpSlices"]

        source = self.currentIncrement["elementResults"]["SDV"][setName]

        destination = self.currentIncrement["elementResults"][job["destination"]][
            setName
        ]

        for ensElType, elements in source.items():
            for elLabel, uelResults in elements.items():
                uelSdv = uelResults["qps"][1]
                qpsData = [uelSdv[qpSlice] for qpSlice in qpSlices]

                destination[ensElType][elLabel]["qps"] = {
                    (i + 1): qpData for i, qpData in enumerate(qpsData)
                }

    def collectExportJobs(self, jobDefinitions):
        """Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        jobs = {}
        for jobDef in jobDefinitions:
            jobName = jobDef["name"]
            dimensions = int(jobDef["dimensions"])
            timeSet = jobDef.get("timeSet", 1)

            jobs[jobName] = EnsightExportJob(
                jobName, dimensions, timeSet, writeEmptyTimeSteps=True
            )

        return jobs

    def collectPerElementJobEntries(self, entryDefinitions, perElementVariableJobs):
        """Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        # jobs = {}
        for entry in entryDefinitions:
            job = perElementVariableJobs[entry["job"]]

            result = entry["result"]
            setName = entry["set"]
            loc = entry["location"]

            which = entry["which"]

            if loc == "qps":
                which = int(w)

            perSetJob = EnsightPerSetJobEntry(
                job,
                setName,
                result=result,
                location=loc,
                which=which,
                extractionSlice=sliceFromString(entry["values"])
                if "values" in entry
                else None,
                extractionFunction=makeExtractionFunction(entry["f(x)"])
                if "f(x)" in entry
                else None,
                offset=None,
            )

            job.entries[setName] = perSetJob

        return perElementVariableJobs

    def collectPerNodeJobEntries(self, entryDefinitions, perNodeVariableJobs):
        """Collect all defined per node jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        # jobs = {}
        for entry in entryDefinitions:
            job = perNodeVariableJobs[entry["job"]]
            setName = entry["set"]

            jobEntry = EnsightPerSetJobEntry(
                job,
                setName,
                result=entry["result"],
                location=None,
                which=None,
                extractionSlice=sliceFromString(entry["values"])
                if "values" in entry
                else None,
                extractionFunction=makeExtractionFunction(entry["f(x)"])
                if "f(x)" in entry
                else None,
                offset=None,
                fillMissingValuesTo=entry.get("fillMissingValuesTo", None),
            )

            job.entries[setName] = jobEntry

        return perNodeVariableJobs

    def createEnsightPerNodeVariableFromPerNodeJob(self, exportJob):
        partsDict = {}
        for setName, jobEntry in exportJob.entries.items():
            elSet = self.elSets[setName]

            print(" {:<20} / {:<28}".format(exportJob.exportName, setName))

            # collect all result, do not yet make a numpy array, as the results array might be ragged, or not present for all nodes
            setNodeIndices = elSet.reducedNodes.keys()

            results = [
                self.currentIncrement["nodeResults"][jobEntry.result].get(node, None)
                for node in setNodeIndices
            ]

            if jobEntry.extractionSlice is not None:
                results = [
                    r[jobEntry.extractionSlice] if r is not None else r for r in results
                ]

            if jobEntry.extractionFunction is not None:
                results = [
                    perSetJobEntry.extractionFunction(r) if r is not None else r
                    for r in results
                ]

            if jobEntry.fillMissingValuesTo is not None:
                defaultResult = np.full(
                    (exportJob.dimensions,), jobEntry.fillMissingValuesTo
                )
                d = exportJob.dimensions

                # first fill up all the results we have
                results = [
                    np.append(
                        r,
                        (jobEntry.fillMissingValuesTo,) * (d - r.shape[0]),
                    )
                    if r is not None and r.shape[0] != d
                    else r
                    for r in results
                ]

                # then set all those we don't have for certain nodes
                results = [defaultResult if r is None else r for r in results]

            try:
                results = np.asarray(results, dtype=float)
            except:
                raise Exception(
                    "Failed to set up all results {:} for all nodes in {:}. Try using fillMissingValuesTo= option?".format(
                        jobEntry.result, setName
                    )
                )

            setVariableDimensions = results.shape[1]

            if setVariableDimensions != exportJob.dimensions:
                raise Exception(
                    "Variable dimension {:} in set {:} does not match the defined job dimension of {:} in job '{:}'. Consider using the 'fillMissingValuesTo' option for the export entry.".format(
                        setVariableDimensions,
                        setName,
                        exportJob.dimensions,
                        exportJob.exportName,
                    )
                )

            partsDict[elSet.ensightPartID] = ("coordinates", results)

        if partsDict or exportJob.writeEmptyTimeSteps:
            return es.EnsightPerNodeVariable(
                exportJob.exportName, exportJob.dimensions, partsDict
            )
        else:
            return None

    def createEnsightPerElementVariableFromPerElementJob(self, exportJob):
        partsDict = {}
        for setName, perSetJobEntry in exportJob.entries.items():
            elSet = self.elSets[setName]
            result = perSetJobEntry.result
            location = perSetJobEntry.location
            which = perSetJobEntry.which

            incrementVariableResults = self.currentIncrement["elementResults"][result][
                setName
            ]
            incrementVariableResultsArrays = {}

            print(" {:<20} / {:<28}".format(exportJob.exportName, setName))

            for ensElType, elDict in incrementVariableResults.items():
                try:
                    results = np.asarray(
                        [
                            elDict[el.label][location][which]
                            for el in elSet.elementsByShape[ensElType]
                        ],
                        dtype=float,
                    )
                except:
                    raise Exception(
                        "Failed to retrieve result '{:}' in '{:}/{:}' for set {:}. Does it exist?".format(
                            result, location, which, setName
                        )
                    )

                if perSetJobEntry.offset:
                    results = results[:, perSetJobEntry.offset :]

                if perSetJobEntry.extractionFunction:
                    results = np.apply_along_axis(
                        perSetJobEntry.extractionFunction, axis=1, arr=results
                    )
                    results = np.reshape(
                        results, (results.shape[0], -1)
                    )  # ensure that dimensions are kept

                if perSetJobEntry.extractionSlice:
                    results = results[:, perSetJobEntry.extractionSlice]

                incrementVariableResultsArrays[ensElType] = results
                setVariableDimensions = results.shape[1]

            if setVariableDimensions != exportJob.dimensions:
                raise Exception(
                    "Variable dimension {:} in set {:} does not match the defined job dimension of {:} in job '{:}'. Consider using the 'fillMissingValuesTo' option for the export entry.".format(
                        setVariableDimensions,
                        setName,
                        exportJob.dimensions,
                        exportJob.exportName,
                    )
                )

            partsDict[elSet.ensightPartID] = incrementVariableResultsArrays

        if partsDict or exportJob.writeEmptyTimeSteps:
            # variableDimension = exportJob.dimensions or variableLength
            return es.EnsightPerElementVariable(
                exportJob.exportName,
                exportJob.dimensions,
                partsDict,
            )
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
            setName = "ALL"

        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]

        self.currentSetName = setName

    def elementHeaderRecord(self, rec):
        elNum = filFlag(rec[0])
        self.currentElementNum = elNum
        self.currentIpt = filFlag(rec[1])

    def handlePerElementOutput(self, rec, result):
        res = filDouble(rec)
        currentIncrement = self.currentIncrement
        currentSetName = self.currentSetName
        currentEnsightElementType = self.currentEnsightElementType
        qp = self.currentIpt
        currentElementNum = self.currentElementNum

        targetLocation = currentIncrement["elementResults"][result][currentSetName][
            currentEnsightElementType
        ]

        if qp not in targetLocation[currentElementNum]["qps"]:
            targetLocation[currentElementNum]["qps"][qp] = res
        else:
            # continuation of an existing record
            targetLocation[currentElementNum]["qps"][qp] = np.concatenate(
                (targetLocation[currentElementNum]["qps"][qp], res)
            )

    def handlePerNodeOutput(self, rec, location):
        node = filInt(rec[0])[0]
        vals = filDouble(rec[1:])

        if location not in self.currentIncrement["nodeResults"]:
            self.currentIncrement["nodeResults"][location] = {}
        self.currentIncrement["nodeResults"][location][node] = vals

    def addNode(self, recordContent):
        label = filInt(recordContent[0])[0]
        coords = filDouble(recordContent[1:4])

        # make coords 3D, always!
        if coords.shape[0] < 3:
            coords = np.pad(coords, (0, 3 - coords.shape[0]), mode="constant")

        if label in self.allNodes:
            print("Node {:} already exists!".format(label))
            exit(0)

        self.allNodes[label] = Node(label, coords)

    def addElement(self, recordContent):
        elNum = filInt(recordContent[0])[0]
        elType = filStrippedString(recordContent[1])
        elabqNodes = filInt(recordContent[2:])

        nodes = [n for n in elabqNodes]

        if elType in self.ignoreLastNodesForElType:
            nodes = nodes[0 : -self.ignoreLastNodesForElType[elType]]

        self.elementDefinitions[elNum] = ElementDefinition(
            elNum, self.ensightElementTypeMappings[elType], nodes
        )

    def addElset(self, recordContent):
        setName = filStrippedString(recordContent[0])
        if not setName:
            setName = "ALL"
        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]
        self.currentSetName = setName
        abqElements = filInt(recordContent[1:])

        self.elSetDefinitions[setName] = ElSetDefinition(
            setName, [e for e in abqElements]
        )

    def contAddElset(self, recordContent):
        abqElements = filInt(recordContent)
        self.elSetDefinitions[self.currentSetName].appendElementLabels(
            [e for e in abqElements]
        )

    def createNodeSetDefinition(self, recordContent):
        setName = filStrippedString(recordContent[0])

        if not setName:
            setName = "ALL"
        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]

        self.currentSetName = setName
        abqNodes = filInt(recordContent[1:])

        self.nSetDefinitions[setName] = NSetDefinition(setName, [n for n in abqNodes])

    def contNodeSetDefinition(self, recordContent):
        abqNodes = filInt(recordContent)
        # self.nSets[self.currentSetName].appendNodes([self.allNodes[n] for n in abqNodes])
        self.nSetDefinitions[self.currentSetName].appendNodeLabels(
            [n for n in abqNodes]
        )

    def addIncrement(self, recordContent):
        self.currentState = "increment parsing"
        r = recordContent
        tTotal, tStep = filDouble(r[0:2])

        nStep, nInc = filInt(r[5:7])
        timeInc = filDouble(r[10])[0]
        currentIncrement = self.currentIncrement
        currentIncrement["tTotal"] = tTotal
        currentIncrement["nInc"] = nInc
        currentIncrement["tStep"] = tStep
        currentIncrement["nStep"] = nStep
        currentIncrement["timeInc"] = timeInc
        currentIncrement["elementResults"] = RecursiveDefaultDict()
        currentIncrement["nodeResults"] = RecursiveDefaultDict()
        print("*" * 80)
        print(
            "processing increment {:>5}  tTotal:{:>16.5f}".format(
                self.nIncrements, self.currentIncrement["tTotal"]
            )
        )

    def addLabelCrossReference(self, recordContent):
        r = recordContent
        intKey = filFlag(r[0])
        label = filStrippedString(r[1:])
        self.labelCrossReferences[str(intKey)] = label

    def surfaceDefHeader(self, recordContent):
        self.currentState = "surface definition"
