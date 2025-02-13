"""
Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np
import src.ensight.ensightgoldformat as es
from src.modeldatabase import Node, NSet, Element, ElSet
from src.misc import sliceFromString, makeExtractionFunction


class _EnsightExportJob:
    def __init__(self, name: str, dimensions: int, timeSetID: int, writeEmptyTimeSteps=True):
        """An Ensight export defines a variable, which might be working on multiple blocks/parts.

        Parameters
        ----------
        name
            The name (of the variable).
        dimensions
            The dimensions of the variable.
        timeSetID
            The ID of the associated time set.
        writeEmptyTimeSteps
            Determine wheter this variable is written even if there is no time associated.
        """

        self.exportName = name
        self.dimensions = dimensions
        self.timeSetID = timeSetID
        self.entries = {}
        self.writeEmptyTimeSteps = writeEmptyTimeSteps


class _EnsightPerSetJobEntry:
    def __init__(
        self,
        job: _EnsightExportJob,
        setType: str,
        setName: str,
        result: str,
        location: str,
        which: str,
        extractionSlice: slice = None,
        extractionFunction: callable = None,
        offset: int = None,
        fillMissingValuesTo=None,
    ):
        """An entry (for a block/part/elset) for an Ensight export job.

        Parameters
        ----------
        job
            The parent EnsightExportJob.
        setType
            'nSet' or 'elSet'.
        setName
            The dimensions of the variable.
        result
            The name of the result.
        location
            Where to find this result (computed, qp, ...).
        which
            Which (number of quadrature point, ...).
        extractionSlice
            (Optional) A slice for extracting the result from a larger array.
        extractionFunction
            (Optional) A function for extracting the result from a larger array.
        offset
            (Optional) An offset.
        fillMissingValuesTo
            (Optional) If the strored results do not match the output dimension, results are filled to this value.
        """

        # self.job = job
        self.setType = setType
        self.setName = setName
        self.extractionSlice = extractionSlice
        self.extractionFunction = extractionFunction
        self.offset = offset
        self.result = result
        self.location = location
        self.which = which
        self.fillMissingValuesTo = fillMissingValuesTo if fillMissingValuesTo is None else float(fillMissingValuesTo)


class EnsightExporter:
    def __init__(self, caseName, inputFile):
        self.ensightCase = es.EnsightChunkWiseCase(".", caseName)
        self.ensightCaseDiscardTimeMarks = False

        self.perElementJobs = self._collectExportJobs(inputFile["*ensightPerElementVariableJob"])
        self.perNodeJobs = self._collectExportJobs(inputFile["*ensightPerNodeVariableJob"])

        self.perElementJobs = self._collectPerElementJobEntries(
            inputFile["*ensightPerElementVariableJobEntry"],
            self.perElementJobs,
        )
        self.perNodeJobs = self._collectPerNodeJobEntries(
            inputFile["*ensightPerNodeVariableJobEntry"],
            self.perNodeJobs,
        )

        self.ensightElementTypeMappings = {x["element"]: x["shape"] for x in inputFile["*defineElementType"]}

        self.ensightElementTypeMappings["node"] = "point"

        self._setToPartIDMapping = {}

        self._nodes = None
        self._elements = None
        self._nSets = None
        self._elSets = None

    def setupModel(
        self,
        nodes: list[Node],
        nSets: dict[str, NSet],
        elements: dict[int, Element],
        elSets: dict[str, ElSet],
    ):
        self._nodes = nodes
        self._elements = elements
        self._nSets = nSets
        self._elSets = elSets

    def setCurrentTime(self, currentTime: float):
        self.ensightCase.setCurrentTime(currentTime)

    def exportGeometry(
        self,
    ):
        geometryTimesetNumber = None
        geometry = self._createEnsightGeometryFromModel(self._nodes, self._nSets, self._elements, self._elSets)
        self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)

    def exportPerNodeVariables(self, nodeResults):
        for exportJob in self.perNodeJobs.values():
            enSightVar = self._createEnsightPerNodeVariableFromPerNodeJob(exportJob, nodeResults)
            if enSightVar:
                self.ensightCase.writeVariableTrendChunk(enSightVar, exportJob.timeSetID)
                del enSightVar

    def exportPerElementVariables(self, elementResults):
        for exportJob in self.perElementJobs.values():
            enSightVar = self._createEnsightPerElementVariableFromPerElementJob(exportJob, elementResults)
            if enSightVar:
                self.ensightCase.writeVariableTrendChunk(enSightVar, exportJob.timeSetID)
                del enSightVar

    def finalize(self, closeFileHandles):
        self.ensightCase.finalize(self.ensightCaseDiscardTimeMarks, closeFileHandles)

    def _collectExportJobs(self, jobDefinitions):
        """Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        jobs = {}
        for jobDef in jobDefinitions:
            jobName = jobDef["name"]
            dimensions = int(jobDef["dimensions"])
            timeSet = jobDef.get("timeSet", 1)

            jobs[jobName] = _EnsightExportJob(jobName, dimensions, timeSet, writeEmptyTimeSteps=True)

        return jobs

    def _collectPerElementJobEntries(self, entryDefinitions, perElementVariableJobs):
        """Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        for entry in entryDefinitions:
            job = perElementVariableJobs[entry["job"]]

            result = entry["result"]
            setName = entry["set"]
            location = entry["location"]

            which = entry["which"]

            if location == "qps":
                which = int(w)

            perSetJob = _EnsightPerSetJobEntry(
                job,
                "elSet",
                setName,
                result=result,
                location=location,
                which=which,
                extractionSlice=sliceFromString(entry["values"]) if "values" in entry else None,
                extractionFunction=makeExtractionFunction(entry["f(x)"]) if "f(x)" in entry else None,
                offset=None,
            )

            job.entries[setName] = perSetJob

        return perElementVariableJobs

    def _collectPerNodeJobEntries(self, entryDefinitions, perNodeVariableJobs):
        """Collect all defined per node jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""

        for entry in entryDefinitions:
            job = perNodeVariableJobs[entry["job"]]
            setType = entry.get("setType", "elSet")

            setName = entry["set"]

            if setType.lower() != "elset" and setType.lower() != "nset":
                raise Exception(
                    "Ensight per node job {:} entry {:}: set type must be either 'elSet' or 'nSet'!".format(
                        job, setName
                    )
                )

            jobEntry = _EnsightPerSetJobEntry(
                job,
                setType,
                setName,
                result=entry["result"],
                location=None,
                which=None,
                extractionSlice=sliceFromString(entry["values"]) if "values" in entry else None,
                extractionFunction=makeExtractionFunction(entry["f(x)"]) if "f(x)" in entry else None,
                offset=None,  # currently not used
                fillMissingValuesTo=entry.get("fillMissingValuesTo", None),
            )

            job.entries[setName] = jobEntry

        return perNodeVariableJobs

    def _createEnsightPerNodeVariableFromPerNodeJob(self, exportJob, nodeResults):
        partsDict = {}
        for i, (setName, jobEntry) in enumerate(exportJob.entries.items()):
            print(" {:<20} ... {:<28}".format(exportJob.exportName if not i else "", setName))
            theSet = None

            if jobEntry.setType == "elSet":
                elSet = self._elSets[setName]
                theSet = elSet

                # collect all result, do not yet make a numpy array, as the results array might be ragged, or not present for all nodes
                setNodeIndices = elSet.reducedNodes.keys()

                results = [nodeResults[jobEntry.result].get(node, None) for node in setNodeIndices]

            else:  # it"s a node set !
                nSet = self._nSets[setName]
                theSet = nSet
                results = [nodeResults[jobEntry.result].get(node.label, None) for node in nSet.nodes]

            if jobEntry.extractionSlice is not None:
                results = [r[jobEntry.extractionSlice] if r is not None else r for r in results]

            if jobEntry.extractionFunction is not None:
                results = [perSetJobEntry.extractionFunction(r) if r is not None else r for r in results]

            if jobEntry.fillMissingValuesTo is not None:
                defaultResult = np.full((exportJob.dimensions,), jobEntry.fillMissingValuesTo)
                d = exportJob.dimensions

                # first fill up all the results we have
                results = [
                    (
                        np.append(
                            r,
                            (jobEntry.fillMissingValuesTo,) * (d - r.shape[0]),
                        )
                        if r is not None and r.shape[0] != d
                        else r
                    )
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

            partsDict[self._setToPartIDMapping[theSet]] = ("coordinates", results)

        if partsDict or exportJob.writeEmptyTimeSteps:
            return es.EnsightPerNodeVariable(exportJob.exportName, exportJob.dimensions, partsDict)
        else:
            return None

    def _createEnsightPerElementVariableFromPerElementJob(self, exportJob, elementResults):
        partsDict = {}
        for i, (setName, perSetJobEntry) in enumerate(exportJob.entries.items()):
            elSet = self._elSets[setName]
            result = perSetJobEntry.result
            location = perSetJobEntry.location
            which = perSetJobEntry.which

            setVariableDimensions = None

            incrementVariableResults = elementResults[result][setName]
            incrementVariableResultsArrays = {}

            print(" {:<20} ... {:<28}".format(exportJob.exportName if not i else "", setName))

            for elType, elDict in incrementVariableResults.items():
                try:
                    results = np.asarray(
                        [elDict[el.label][location][which] for el in elSet.elementsByShape[elType]],
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

                if perSetJobEntry.extractionSlice:
                    results = results[:, perSetJobEntry.extractionSlice]

                if perSetJobEntry.extractionFunction:
                    results = np.apply_along_axis(perSetJobEntry.extractionFunction, axis=1, arr=results)
                    results = np.reshape(results, (results.shape[0], -1))  # ensure that dimensions are kept

                incrementVariableResultsArrays[elType] = results
                setVariableDimensions = results.shape[1]

            if not setVariableDimensions:
                raise Exception(
                    "No results for set {:} in job '{:}'. ".format(
                        setName,
                        exportJob.exportName,
                    )
                )

            if setVariableDimensions != exportJob.dimensions:
                raise Exception(
                    "Variable dimension {:} in set {:} does not match the defined job dimension of {:} in job '{:}'. Consider using the 'fillMissingValuesTo' option for the export entry.".format(
                        setVariableDimensions,
                        setName,
                        exportJob.dimensions,
                        exportJob.exportName,
                    )
                )

            partsDict[self._setToPartIDMapping[elSet]] = incrementVariableResultsArrays

        if partsDict or exportJob.writeEmptyTimeSteps:
            return es.EnsightPerElementVariable(
                exportJob.exportName,
                exportJob.dimensions,
                partsDict,
                self.ensightElementTypeMappings,
            )
        else:
            return None

    def _createEnsightGeometryFromModel(self, nodes: list[Node], nSets, elements: dict, elSets: dict[str, ElSet]):
        partList = []
        partNumber = 1

        for elSet in elSets.values():
            elSetPart = es.EnsightUnstructuredPart(
                elSet.name,
                partNumber,
                elSet.reducedElements,
                elSet.reducedNodeCoords3D,
                list(elSet.reducedNodes.keys()),
                self.ensightElementTypeMappings,
            )
            # elSet.ensightPartID = partNumber
            self._setToPartIDMapping[elSet] = partNumber
            partList.append(elSetPart)
            partNumber += 1

        for nSet in nSets.values():
            nSetPart = es.EnsightUnstructuredPart(
                "NSET_" + nSet.name,
                partNumber,
                {
                    "node": [
                        (
                            n.label,
                            [
                                i,
                            ],
                        )
                        for i, n in enumerate(nSet.nodes)
                    ]
                },
                np.array([n.coords for n in nSet.nodes]),
                [n.label for n in nSet.nodes],
                self.ensightElementTypeMappings,
            )
            # nSet.ensightPartID = partNumber
            partList.append(nSetPart)
            self._setToPartIDMapping[nSet] = partNumber
            partNumber += 1

        geometry = es.EnsightGeometry("geometry", "-", "-", partList, "given", "given")
        return geometry
