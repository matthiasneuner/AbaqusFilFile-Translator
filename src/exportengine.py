"""
Created on Fri Oct 28 13:43:53 2016

Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np
from collections import defaultdict
from src.ensight.ensightexporter import EnsightExporter
from src.modeldatabase import Node, NSet, Element, ElSet
from src.misc import RecursiveDefaultDict
from prettytable import PrettyTable


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


class _ElementDefinition:
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


class _ElSetDefinition:
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


class _NSetDefinition:
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


class ExportEngine:
    def __init__(self, inputFile: dict, exportName: str):
        """This is the export engine. It parses a .fil file record wise,
        and exports results based on user defined jobs.

        Parameters
        ----------
        inputFile
            The dictionary containing the input file.
        exportName
            The export name.
        """

        self.uelSdvToQpJobs = self.collectUelSDVToQpJobs(inputFile["*UELSDVToQuadraturePoints"])
        self.qpAverageJobs = self.collectQpAverageJobs(inputFile["*computeAverageOverQuadraturePoints"])

        self.ignoreLastNodesForElType = {x["element"]: x["number"] for x in inputFile["*ignoreLastNodesForElementType"]}

        self.ensightExporter = EnsightExporter(exportName, inputFile)

        self.nodes = {}
        # add a default node, to which abaqus falls back if it creates node in place (e.g, for hex27 elements in contact)
        self.nodes[0] = Node(0, np.array([0.0, 0.0, 0.0]))

        self.elementDefinitions = {}
        self.elSetDefinitions = {}
        self.nSetDefinitions = {}

        self.elements = {}
        self.nSets = {}
        self.elSets = {}
        self._substituteElSets = self._assembleSubsitutionElSets(inputFile)

        self.currentState = "model setup"
        self.currentIncrement = {}
        self.currentAbqElSet = None
        self.currentSetName = "ALL"
        self.currentIpt = 1
        self.nIncrements = 0
        self.timeHistory = []
        self.labelCrossReferences = {}

        self.knownRecords = {
            1: ("Element header record", self._elementHeaderRecord),
            5: ("SDV output", lambda x: self._handlePerElementOutput(x, "SDV")),
            11: ("S output", lambda x: self._handlePerElementOutput(x, "S")),
            21: ("E output", lambda x: self._handlePerElementOutput(x, "E")),
            22: ("PE output", lambda x: self._handlePerElementOutput(x, "PE")),
            85: ("Local coordinate system(?)", lambda x: None),
            89: ("LE output", lambda x: self._handlePerElementOutput(x, "LE")),
            101: ("U output", lambda x: self._handlePerNodeOutput(x, "U")),
            102: ("V output", lambda x: self._handlePerNodeOutput(x, "V")),
            103: ("A output", lambda x: self._handlePerNodeOutput(x, "A")),
            108: ("POR output", lambda x: self._handlePerNodeOutput(x, "POR")),
            104: ("RF output", lambda x: self._handlePerNodeOutput(x, "RF")),
            201: ("NT output", lambda x: self._handlePerNodeOutput(x, "NT")),
            1501: ("Surface definition header", self._surfaceDefHeader),
            1502: ("Surface facet", lambda x: None),
            1900: ("element definition", self._addElementDefinition),
            1901: ("node definition", self._addNode),
            1902: ("active dof", lambda x: None),
            1911: ("output request definition", self._outputDefinition),
            1921: ("heading", self._printHeading1921),
            1922: ("heading", lambda x: None),
            1931: ("node set definition", self._createNodeSetDefinition),
            1932: ("node set definition cont.", self._contNodeSetDefinition),
            1933: ("element set definition", self._addElsetDefinition),
            1934: ("element set definition cont.", self._contAddElset),
            1940: ("label cross reference", self._addLabelCrossReference),
            1999: ("total energies", self._printEnergies),
            2000: ("start increment", self._addIncrement),
            2001: ("end increment", self._finishAndParseIncrement),
        }

    def computeRecord(self, recordLength: int, recordType: int, recordContent: np.ndarray):
        """The main function of the export engine. It computes a .fil file record.

        Parameters
        ----------
        recordLength
            The length of the .fil record content.
        recordType
            The unique integer defining the type of record according to the Abaqus documentation.
        recordContent
            The data of the .fil record.

        Returns
        -------
        bool
            Success of the record computation.
        """

        if recordType in self.knownRecords:
            doc, action = self.knownRecords[recordType]
            action(recordContent)
            return True
        else:
            print("{:<20}{:>6}{:>10}{:>4}".format("unknown record:", recordType, " of length", recordLength))
            return False

    def _finishAndParseIncrement(self, recordContent: np.ndarray):
        """A .fil increment (or the model setup) is finished. Time to write results and geometry!

        Parameters
        ----------
        recordContent
            The fil record. It is empty.
        """

        if self.currentState == "model setup":
            # we always create the 'ALL' set
            ALLSet = _ElSetDefinition("ALL", list(self.elementDefinitions.keys()))
            self.elSetDefinitions["ALL"] = ALLSet

            # Time to create Elements, ElSets and NodeSets from the definitions!
            for elDef in self.elementDefinitions.values():
                self.elements[elDef.label] = Element(
                    elDef.label,
                    elDef.shape,
                    [self.nodes[label] for label in elDef.nodeLabels],
                )
            # replace all label references by the respective labels
            for key, label in self.labelCrossReferences.items():
                strKey = key  # str(key)
                if strKey in self.nSetDefinitions:
                    self.nSetDefinitions[label] = self.nSetDefinitions[strKey]
                    self.nSetDefinitions[label].name = label
                    del self.nSetDefinitions[strKey]
                if strKey in self.elSetDefinitions:
                    self.elSetDefinitions[label] = self.elSetDefinitions[strKey]
                    self.elSetDefinitions[label].name = label
                    del self.elSetDefinitions[strKey]

            self.elSetDefinitions.update(self._substituteElSets)
            for elSetDef in self.elSetDefinitions.values():

                try:
                    self.elSets[elSetDef.name] = ElSet(
                        elSetDef.name,
                        [self.elements[label] for label in elSetDef.elementLabels],
                    )
                except KeyError as e:
                    print("Element set {:} not created!".format(elSetDef.name))
                    print("Please check following element labels: {:}".format([int(i) for i in elSetDef.elementLabels]))
                    print(
                        "For Abaqus/Explicit, it is observed that the definition of some element sets is faulty if multiple CPUs are used."
                    )
                    print("Error: {:}".format(e))
                    continue

            for nSetDef in self.nSetDefinitions.values():
                self.nSets[nSetDef.name] = NSet(nSetDef.name, [self.nodes[n] for n in nSetDef.nodeLabels])

            self.ensightExporter.setupModel(self.nodes, self.nSets, self.elements, self.elSets)
            self.ensightExporter.exportGeometry()

        elif self.currentState == "surface definition":
            pass

        elif self.currentState == "increment parsing":
            self.nIncrements += 1
            self.ensightExporter.setCurrentTime(self.currentIncrement["tTotal"])
            self.timeHistory.append(self.currentIncrement["tTotal"])

            print("increment contains element results for")
            print(
                "\n".join(
                    [
                        " {:5} [{:}]".format(resName, ", ".join([s for s in resEntries]))
                        for resName, resEntries in self.currentIncrement["elementResults"].items()
                    ]
                )
            )
            print("")
            print("increment contains node results for")
            print(
                "\n".join(
                    [
                        " {:5} [{:10} nodes]".format(resName, len(resEntries))
                        for resName, resEntries in self.currentIncrement["nodeResults"].items()
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

            self.ensightExporter.exportPerNodeVariables(self.currentIncrement["nodeResults"])
            self.ensightExporter.exportPerElementVariables(self.currentIncrement["elementResults"])

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
        self.ensightExporter.finalize(closeFileHandles=True)

    def collectQpAverageJobs(self, entries):
        """Commonly, the average of a result over all quadrature points per element should be computed.
        This function gathers all jobs.

        Parameters
        ----------
        entries
            The list of job definitions."""

        jobs = []
        for entry in entries:
            jobs.append(entry)
        return jobs

    def computeQpAverage(self, job: dict):
        """Compute the average of an  elemental variable over all quadrature points.

        Parameters
        ----------
        job
            The job defintion."""

        result = job["result"]
        setName = job["set"]

        setResults = self.currentIncrement["elementResults"][result][setName]

        for elTypeResults in setResults.values():
            for elResults in elTypeResults.values():
                elResults["computed"]["average"] = np.mean([qpRes for qpRes in elResults["qps"].values()], axis=0)

    def collectUelSDVToQpJobs(self, entries: list):
        """Abaqus UEL SDVs commonly should be computed to something resonable!
        This function gathers the respective jobs from the input file.

        Parameters
        ----------
        entries
            The list of job definitions."""

        jobs = []

        for entry in entries:
            offset = entry["qpInitialOffset"]
            nQps = entry["qpCount"]
            qpDistance = entry["qpDistance"]
            entry["qpSlices"] = [slice(offset + i * qpDistance, offset + (i + 1) * qpDistance) for i in range(nQps)]

            jobs.append(entry)

        return jobs

    def computeUelSdvToQp(self, job: dict):
        """Abaqus UEL SDVs commonly should be computed to something resonable!

        Parameters
        ----------
        job
            The job definition containing the task info."""

        setName = job["set"]
        destination = job["destination"]
        qpSlices = job["qpSlices"]

        source = self.currentIncrement["elementResults"]["SDV"][setName]

        destination = self.currentIncrement["elementResults"][job["destination"]][setName]

        for ensElType, elements in source.items():
            for elLabel, uelResults in elements.items():
                uelSdv = uelResults["qps"][1]
                qpsData = [uelSdv[qpSlice] for qpSlice in qpSlices]

                destination[ensElType][elLabel]["qps"] = {(i + 1): qpData for i, qpData in enumerate(qpsData)}

    def _outputDefinition(self, recordContent: np.ndarray):
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
            self.currentElementType = elType

        elif flag == 1:
            setName = filStrippedString(recordContent[1])

        if not setName:
            setName = "ALL"

        elif setName in self.labelCrossReferences:
            setName = self.labelCrossReferences[setName]

        self.currentSetName = setName

    def _elementHeaderRecord(self, recordContent: np.ndarray):
        """Initialize the element we are working on.

        Parameters
        ----------
        recordContent
            The fil record. Contains the element label, and the current quadrature point.
        """

        elLabel = filFlag(recordContent[0])
        self.currentElementLabel = elLabel
        self.currentIpt = filFlag(recordContent[1])

    def _handlePerElementOutput(self, recordContent: np.ndarray, result: str):
        """Data for an element.

        Parameters
        ----------
        recordContent
            The fil record. Contains the data.
        result
            The result type (e.g., S,E,SDV...).
        """

        res = filDouble(recordContent)
        currentIncrement = self.currentIncrement
        currentSetName = self.currentSetName
        currentElementType = self.currentElementType
        qp = self.currentIpt
        currentElementLabel = self.currentElementLabel

        targetLocation = currentIncrement["elementResults"][result][currentSetName][currentElementType]

        if qp not in targetLocation[currentElementLabel]["qps"]:
            targetLocation[currentElementLabel]["qps"][qp] = res
        else:
            # continuation of an existing record
            targetLocation[currentElementLabel]["qps"][qp] = np.concatenate(
                (targetLocation[currentElementLabel]["qps"][qp], res)
            )

    def _handlePerNodeOutput(self, recordContent: np.ndarray, result: str):
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

    def _addNode(self, recordContent: np.ndarray):
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
            print(
                "Node {:<6} {:} alrdy dfnd at {:}; ignoring".format(
                    label,
                    np.array2string(coords, precision=2, floatmode="fixed"),
                    np.array2string(self.nodes[label].coords, precision=2, floatmode="fixed"),
                )
            )

        self.nodes[label] = Node(label, coords)

    def _addElementDefinition(self, recordContent):
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

        self.elementDefinitions[elNum] = _ElementDefinition(elNum, elType, nodes)

    def _addElsetDefinition(self, recordContent):
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

        self.elSetDefinitions[setName] = _ElSetDefinition(setName, [e for e in abqElements])

    def _contAddElset(self, recordContent):
        """Add more element labels to a element set definition.

        Parameters
        ----------
        recordContent
            The fil record. Contains the element labels.
        """

        abqElements = filInt(recordContent)
        self.elSetDefinitions[self.currentSetName].appendElementLabels([e for e in abqElements])

    def _createNodeSetDefinition(self, recordContent):
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

        self.nSetDefinitions[setName] = _NSetDefinition(setName, [n for n in abqNodes])

    def _contNodeSetDefinition(self, recordContent):
        """Add more nodes labels to a node set definition.

        Parameters
        ----------
        recordContent
            The fil record. Contains The node labels.
        """

        abqNodes = filInt(recordContent)
        self.nSetDefinitions[self.currentSetName].appendNodeLabels([n for n in abqNodes])

    def _addIncrement(self, recordContent):
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

        # a level 4 RecursiveDefaultDict
        currentIncrement["elementResults"] = RecursiveDefaultDict(4)  # result / set / shape / element number / location

        # a level 2 RecursiveDefaultDict
        currentIncrement["nodeResults"] = RecursiveDefaultDict(2)  # result / node

        print("+" + "-" * 78 + "+")
        print(
            "| processing increment {:>5} | step time:{:>11.5f} | total time:{:>12.5f} |".format(
                self.nIncrements,
                self.currentIncrement["tStep"],
                self.currentIncrement["tTotal"],
            )
        )
        print("+" + "-" * 78 + "+")

    def _addLabelCrossReference(self, recordContent):
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

    def _surfaceDefHeader(self, recordContent):
        """A definition of a side set follows.
        Currently not used!

        Parameters
        ----------
        recordContent
            The fil record. Containts the surface information.
        """
        self.currentState = "surface definition"

    def _printEnergies(self, recordContent):
        """Print the total energies.

        Parameters
        ----------
        recordContent
            The fil record. Contains the energy information.
        """
        energyTypes = [
            "Total kinetic energy (ALLKE).",
            "Total recoverable (elastic) strain energy (ALLSE).",
            "Total external work (ALLWK).",
            "Total plastic dissipation (ALLPD).",
            "Total viscoelastic dissipation (ALLCD).",
            "Total viscous dissipation (ALLVD).",
            "Total loss of kinetic energy at impacts (ALLKL) (S).",
            "Total artificial strain energy (ALLAE).",
            "Total distortion control dissipation energy (ALLDC).",
            "Total electrostatic energy (ALLEE) (S).",
            "Total strain energy (ALLIE).",
            "Total energy balance (ETOTAL).",
            "Total energy dissipated through frictional effects (ALLFD).",
            "Total electrical energy dissipated in conductors (ALLJD) (S).",
            "Percent change in mass (DMASS).",
            "Total damage dissipation (ALLDMD).",
            "Internal heat energy (ALLIHE) (E).",
            "External heat energy (ALLHF) (E).",
        ]

        # make a pretty table with a width of 40 characters
        t = PrettyTable(max_width=80, max_table_width=80)
        t.field_names = ["Energy type", "Value"]
        for value, energyType in zip(filDouble(recordContent), energyTypes):
            t.add_row([energyType, value])

        # format values with only 5 digits
        t.float_format = "5.2"

        print(t)

    def _printHeading1921(self, recordContent):
        """Print the heading of the .fil file.

        Parameters
        ----------
        recordContent
            The fil record. Contains the heading information.
        """
        abqRelease = filStrippedString(recordContent[0])
        date = filStrippedString(recordContent[1:3])
        time = filStrippedString(recordContent[3])
        nElements = filInt(recordContent[4])[0]
        nNodes = filInt(recordContent[5])[0]
        elLength = filDouble(recordContent[6])[0]

        t = PrettyTable(max_width=80, max_table_width=80)
        t = PrettyTable(min_width=80, min_table_width=80)
        t.field_names = ["Abaqus release", "Date", "Time", "elements", "nodes"]
        t.add_row([abqRelease, date, time, nElements, nNodes])
        t.float_format = "5.2"

        print(t)

    def _assembleSubsitutionElSets(self, inputFile: dict):

        elementSets = dict()
        for elSetDefinition in inputFile["*substituteElSet"]:
            name = elSetDefinition["elSet"]

            data = elSetDefinition["data"]
            elNumbers = [int(num) for line in data for num in line]

            elementSets[name] = _ElSetDefinition(name, elNumbers)

        return elementSets
