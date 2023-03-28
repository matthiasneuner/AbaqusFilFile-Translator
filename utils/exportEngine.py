"""
Created on Fri Oct 28 13:43:53 2016

Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np
from collections import defaultdict
from utils.ensightExporter import EnsightExporter


def filInt(word: np.ndarray):
    """Convert a. fil word to a 8 byte integer.

    Parameters
    ----------
    word
        The fil word.

    Returns
    -------
    type
        The integer.
    """

    return word.view("<i8").ravel()


def filString(word):
    """Convert a. fil word to string.

    Parameters
    ----------
    word
        The fil word.

    Returns
    -------
    type
        The string.
    """

    return word.view("a8")


def filStrippedString(word):
    return filString(word).tostring().decode("utf-8").strip()


def filDouble(word):
    """Convert a. fil word to double precision float.

    Parameters
    ----------
    word
        The fil word.

    Returns
    -------
    type
        The float.
    """

    return word.view("<d").ravel()


def filFlag(word):
    """Convert a. fil word to a flag.

    Parameters
    ----------
    word
        The fil word.

    Returns
    -------
    type
        The flag.
    """

    return word[0:4].view("<i")[0]


def RecursiveDefaultDict():
    return defaultdict(RecursiveDefaultDict)


class ElementDefinition:
    def __init__(self, label: int, elementType: str, nodeLabels: list[int]):
        """A description of an element, not a discrete instance.

        Parameters
        ----------
        label
            The label of the element.
        elementType
            The shape of the element.
        nodeLabels
            The list of node labels for this element.
        """

        self.label = label
        self.shape = elementType
        self.nodeLabels = nodeLabels


class ElSetDefinition:
    def __init__(self, name: str, elementLabels: list[int]):
        """A description of an element set, not a discrete instance.

        Parameters
        ----------
        name
            The name of the set.
        elementLabels
            The list of element labels for this set.
        """

        self.name = name
        self.elementLabels = []
        self.appendElementLabels(elementLabels)

    def appendElementLabels(self, elementLabels: list[int]):
        """Add more element labels to the set description.

        Parameters
        ----------
        elementLabels
            The list of element labels to be added.
        """

        self.elementLabels += elementLabels


class Node:
    def __init__(self, label: int, coords: np.array):
        """A spatial node.

        Parameters
        ----------
        label
            The label of this node.
        coords
            The coordinates of this node.
        """

        self.label = label
        self.coords = coords


class Element:
    def __init__(self, label: int, shape: str, nodes: list[Node]):
        """A discrete element instance.

        Parameters
        ----------
        label
            The label of this element.
        shape
            The shape of this element.
        nodes
            The list of nodes.
        """

        self.label = label
        self.shape = shape
        self.nodes = nodes


class ElSet:
    def __init__(self, name: str, elements: list[Element]):
        """A discrete element set instance.

        Parameters
        ----------
        name
            The name of this set.
        elements
            The list of elements in this set.
        ensightPartID
            The ensight ID.
        """

        self.name = name
        # self.ensightPartID = ensightPartID

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
    def __init__(self, name: str, nodeLabels: list[int]):
        """A description of a node set, not a discrete instance.

        Parameters
        ----------
        name
            The name of the set.
        nodeLabels
            The list of node labels for this set.
        """
        self.name = name
        self.nodeLabels = nodeLabels

    def appendNodeLabels(self, nodeLabels: list[int]):
        """Add more node labels to this set

        Parameters
        ----------
        nodeLabels
            The list of node labels to be added.
        """
        self.nodeLabels += nodeLabels


class NSet:
    def __init__(self, name, nodes):
        """A discrete node set instance.

        Parameters
        ----------
        name
            The name of this set.
        nodes
            The list of nodes in this set.
        """
        self.name = name
        self.nodes = nodes


class ExportEngine:
    def __init__(self, inputFile, caseName):
        self.uelSdvToQpJobs = self.collectUelSDVToQpJobs(
            inputFile["*UELSDVToQuadraturePoints"]
        )
        self.qpAverageJobs = self.collectQpAverageJobs(
            inputFile["*computeAverageOverQuadraturePoints"]
        )

        self.ignoreLastNodesForElType = {
            x["element"]: x["number"]
            for x in inputFile["*ignoreLastNodesForElementType"]
        }

        self.ensightElementTypeMappings = {
            x["element"]: x["shape"] for x in inputFile["*defineElementType"]
        }

        self.ensightExporter = EnsightExporter(caseName, inputFile)

        self.nodes = {}
        # add a default node, to which abaqus falls back if it creates node in place (e.g, for hex27 elements in contact)
        self.nodes[0] = Node(0, np.array([0.0, 0.0, 0.0]))

        self.elementDefinitions = {}
        self.elSetDefinitions = {}
        self.nSetDefinitions = {}

        self.elements = {}
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
            1900: ("element definition", self.addElementDefinition),
            1901: ("node definition", self.addNode),
            1902: ("active dof", lambda x: None),
            1911: ("output request definition", self.outputDefinition),
            1921: ("heading", lambda x: None),
            1922: ("heading", lambda x: None),
            1931: ("node set definition", self.createNodeSetDefinition),
            1932: ("node set definition cont.", self.contNodeSetDefinition),
            1933: ("element set definition", self.addElsetDefinition),
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

    def finishAndParseIncrement(self, recordContent: np.ndarray):
        """A .fil increment (or the model setup) is finished. Time to write results and geometry!

        Parameters
        ----------
        recordContent
            The fil record. It is empty.
        """

        if self.currentState == "model setup":
            ALLSet = ElSetDefinition("ALL", list(self.elementDefinitions.keys()))
            self.elSetDefinitions["ALL"] = ALLSet

            # Time to create Elements, ElSets and NodeSets from the definitions!
            for elDef in self.elementDefinitions.values():
                self.elements[elDef.label] = Element(
                    elDef.label,
                    elDef.shape,
                    [self.nodes[label] for label in elDef.nodeLabels],
                )
                # el.nodes = [self.allNodes[label] for label in el.nodes]

            for elSetDef in self.elSetDefinitions.values():
                self.elSets[elSetDef.name] = ElSet(
                    elSetDef.name,
                    [self.elements[label] for label in elSetDef.elementLabels],
                )

            for nSetDef in self.nSetDefinitions.values():
                self.nSets[nSetDef.name] = NSet(
                    nSetDef.name, [self.nodes[n] for n in nSetDef.nodeLabels]
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

            self.ensightExporter.setupModel(
                self.nodes, self.nSets, self.elements, self.elSets
            )
            self.ensightExporter.exportGeometry()

        elif self.currentState == "surface definition":
            pass

        elif self.currentState == "increment parsing":
            self.nIncrements += 1
            self.ensightExporter.setCurrentTime(self.currentIncrement["tTotal"])
            self.timeHistory.append(self.currentIncrement["tTotal"])

            print("*" * 80)
            print("increment contains element results for")
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
            print("increment contains node results for")
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

            # operate on elemetal results (e.g. compute average over quadraturePoint )
            for uelSdvToQpJob in self.uelSdvToQpJobs:
                self.computeUelSdvToQp(uelSdvToQpJob)

            for qpAverageJob in self.qpAverageJobs:
                self.computeQpAverage(qpAverageJob)

            self.ensightExporter.exportPerNodeVariables(
                self.currentIncrement["nodeResults"]
            )
            self.ensightExporter.exportPerElementVariables(
                self.currentIncrement["elementResults"]
            )

            if self.nIncrements % 10 == 0:
                # intermediate saving ...

                self.ensightExporter.finalize(
                    closeFileHandles=False,
                )

            # data might consume a lot of memory, so we delete it (explicitly)
            del self.currentIncrement
            # and create a new dict for a new increment; it is filled with the next start increment entry in the fil file
            self.currentIncrement = dict()

    def finalize(self):
        self.ensightExporter.finalize()

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

    def outputDefinition(self, recordContent: np.ndarray):
        """Initialize a new output we are working on.

        Parameters
        ----------
        recordContent
            The fil record. Contains the element label, and the current quadrature point.
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

    def elementHeaderRecord(self, recordContent: np.ndarray):
        """Initialize the element we are working on.

        Parameters
        ----------
        recordContent
            The fil record. Contains the element label, and the current quadrature point.
        """

        elNum = filFlag(recordContent[0])
        self.currentElementNum = elNum
        self.currentIpt = filFlag(recordContent[1])

    def handlePerElementOutput(self, rec, result: str):
        """Data for an element.

        Parameters
        ----------
        recordContent
            The fil record. Contains the data.
        result
            The result type (e.g., S,E,SDV...).
        """

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

    def handlePerNodeOutput(self, recordContent: np.ndarray, result: str):
        """Data for a node.

        Parameters
        ----------
        recordContent
            The fil record. Contains the node label, and the data.
        result
            The result (e.g., U, NT, ...)
        """

        node = filInt(recordContent[0])[0]
        vals = filDouble(recordContent[1:])

        if result not in self.currentIncrement["nodeResults"]:
            self.currentIncrement["nodeResults"][result] = {}
        self.currentIncrement["nodeResults"][result][node] = vals

    def addNode(self, recordContent: np.ndarray):
        """Definition of a node.

        Parameters
        ----------
        recordContent
            The fil record. Contains the label, and the coordinates.
        """

        label = filInt(recordContent[0])[0]
        coords = filDouble(recordContent[1:4])

        # make coords 3D, always!
        if coords.shape[0] < 3:
            coords = np.pad(coords, (0, 3 - coords.shape[0]), mode="constant")

        if label in self.nodes:
            print("Node {:} already exists!".format(label))
            exit(0)

        self.nodes[label] = Node(label, coords)

    def addElementDefinition(self, recordContent):
        """Definition of an element.

        Parameters
        ----------
        recordContent
            The fil record. Contains the label, the type and the node labels.
        """

        elNum = filInt(recordContent[0])[0]
        elType = filStrippedString(recordContent[1])
        elabqNodes = filInt(recordContent[2:])

        nodes = [n for n in elabqNodes]

        if elType in self.ignoreLastNodesForElType:
            nodes = nodes[0 : -self.ignoreLastNodesForElType[elType]]

        self.elementDefinitions[elNum] = ElementDefinition(
            elNum, self.ensightElementTypeMappings[elType], nodes
        )

    def addElsetDefinition(self, recordContent):
        """Definition of an Abaqus element set.

        Parameters
        ----------
        recordContent
            The fil record. Contains the name and the element labels.
        """

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
        """Add more element labels to a element set definition.

        Parameters
        ----------
        recordContent
            The fil record. Contains the element labels.
        """

        abqElements = filInt(recordContent)
        self.elSetDefinitions[self.currentSetName].appendElementLabels(
            [e for e in abqElements]
        )

    def createNodeSetDefinition(self, recordContent):
        """Definition of an Abaqus node set.

        Parameters
        ----------
        recordContent
            The fil record. Contains the name and the node labels.
        """

        setName = filStrippedString(recordContent[0])

        if not setName:
            setName = "ALL"
        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]

        self.currentSetName = setName
        abqNodes = filInt(recordContent[1:])

        self.nSetDefinitions[setName] = NSetDefinition(setName, [n for n in abqNodes])

    def contNodeSetDefinition(self, recordContent):
        """Add more nodes labels to a node set definition.

        Parameters
        ----------
        recordContent
            The fil record. Contains The node labels.
        """

        abqNodes = filInt(recordContent)
        self.nSetDefinitions[self.currentSetName].appendNodeLabels(
            [n for n in abqNodes]
        )

    def addIncrement(self, recordContent):
        """A .fil increment (or the model setup) starts. We prepare everything.

        Parameters
        ----------
        recordContent
            The fil record. Contains time information.
        """

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
        """Reference to a label using an integer.

        Parameters
        ----------
        recordContent
            The fil record. Contains the label information.
        """

        r = recordContent
        intKey = filFlag(r[0])
        label = filStrippedString(r[1:])
        self.labelCrossReferences[str(intKey)] = label

    def surfaceDefHeader(self, recordContent):
        self.currentState = "surface definition"
