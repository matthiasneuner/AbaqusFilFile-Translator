import numpy as np
import utils.ensightgoldformat as es


def sliceFromString(string, shift=0):
    """generate a slice from a string, which can represent a slice or an index"""
    if ":" in string:
        a, b = string.split(":")
        return slice(int(a) + shift, int(b) + shift)
    else:
        return slice(int(string) + shift, int(string) + 1 + shift)


def makeExtractionFunction(expression, symbol="x"):
    """make a simple f(x) expression from string"""
    return lambda x: eval(expression, globals(), {symbol: x})


class _EnsightExportJob:
    def __init__(self, name, dimensions, timeSetID, writeEmptyTimeSteps=True):
        self.exportName = name
        self.dimensions = dimensions
        self.timeSetID = timeSetID
        self.entries = {}
        self.writeEmptyTimeSteps = writeEmptyTimeSteps


class _EnsightPerSetJobEntry:
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


class EnsightExporter:
    def __init__(self, caseName, inputFile):
        self.ensightCase = es.EnsightChunkWiseCase(".", caseName)
        self.ensightCaseDiscardTimeMarks = False

        self.perElementJobs = self._collectExportJobs(
            inputFile["*ensightPerElementVariableJob"]
        )
        self.perNodeJobs = self._collectExportJobs(
            inputFile["*ensightPerNodeVariableJob"]
        )

        self.perElementJobs = self._collectPerElementJobEntries(
            inputFile["*ensightPerElementVariableJobEntry"],
            self.perElementJobs,
        )
        self.perNodeJobs = self._collectPerNodeJobEntries(
            inputFile["*ensightPerNodeVariableJobEntry"],
            self.perNodeJobs,
        )

        self.ensightElementTypeMappings = {
            x["element"]: x["shape"] for x in inputFile["*defineElementType"]
        }

        # TODO determine if we can do that at instantiation
        self.nodes = None
        self.elements = None
        self.nSets = None
        self.elSets = None

    def setupModel(self, nodes, nSets, elements, elSets):
        self.nodes = nodes
        self.elements = elements
        self.nSets = nSets
        self.elSets = elSets

    def setCurrentTime(self, currentTime: float):
        self.ensightCase.setCurrentTime(currentTime)

    def exportGeometry(
        self,
    ):
        geometryTimesetNumber = None
        geometry = self._createEnsightGeometryFromModel(
            self.nodes, self.nSets, self.elements, self.elSets
        )
        self.ensightCase.writeGeometryTrendChunk(geometry, geometryTimesetNumber)

    def exportPerNodeVariables(self, nodeResults):
        for exportJob in self.perNodeJobs.values():
            enSightVar = self._createEnsightPerNodeVariableFromPerNodeJob(
                exportJob, nodeResults
            )
            if enSightVar:
                self.ensightCase.writeVariableTrendChunk(
                    enSightVar, exportJob.timeSetID
                )
                del enSightVar

    def exportPerElementVariables(self, elementResults):
        for exportJob in self.perElementJobs.values():
            enSightVar = self._createEnsightPerElementVariableFromPerElementJob(
                exportJob, elementResults
            )
            if enSightVar:
                self.ensightCase.writeVariableTrendChunk(
                    enSightVar, exportJob.timeSetID
                )
                del enSightVar

    def finalize(self):
        self.ensightCase.finalize(discardTimeMarks=self.ensightCaseDiscardTimeMarks)

    def _collectExportJobs(self, jobDefinitions):
        """Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        jobs = {}
        for jobDef in jobDefinitions:
            jobName = jobDef["name"]
            dimensions = int(jobDef["dimensions"])
            timeSet = jobDef.get("timeSet", 1)

            jobs[jobName] = _EnsightExportJob(
                jobName, dimensions, timeSet, writeEmptyTimeSteps=True
            )

        return jobs

    def _collectPerElementJobEntries(self, entryDefinitions, perElementVariableJobs):
        """Collect all defined per element jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""
        for entry in entryDefinitions:
            job = perElementVariableJobs[entry["job"]]

            result = entry["result"]
            setName = entry["set"]
            loc = entry["location"]

            which = entry["which"]

            if loc == "qps":
                which = int(w)

            perSetJob = _EnsightPerSetJobEntry(
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

    def _collectPerNodeJobEntries(self, entryDefinitions, perNodeVariableJobs):
        """Collect all defined per node jobs in a dictionary
        a job can consist of multiple Ensight variable definitions
        with the same name defined on different element sets"""

        for entry in entryDefinitions:
            job = perNodeVariableJobs[entry["job"]]
            setName = entry["set"]

            jobEntry = _EnsightPerSetJobEntry(
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

    def _createEnsightPerNodeVariableFromPerNodeJob(self, exportJob, nodeResults):
        partsDict = {}
        for setName, jobEntry in exportJob.entries.items():
            elSet = self.elSets[setName]

            print(" {:<20} / {:<28}".format(exportJob.exportName, setName))

            # collect all result, do not yet make a numpy array, as the results array might be ragged, or not present for all nodes
            setNodeIndices = elSet.reducedNodes.keys()

            results = [
                nodeResults[jobEntry.result].get(node, None) for node in setNodeIndices
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

    def _createEnsightPerElementVariableFromPerElementJob(
        self, exportJob, elementResults
    ):
        partsDict = {}
        for setName, perSetJobEntry in exportJob.entries.items():
            elSet = self.elSets[setName]
            result = perSetJobEntry.result
            location = perSetJobEntry.location
            which = perSetJobEntry.which

            incrementVariableResults = elementResults[result][setName]
            incrementVariableResultsArrays = {}

            print(" {:<20} / {:<28}".format(exportJob.exportName, setName))

            for elType, elDict in incrementVariableResults.items():
                try:
                    results = np.asarray(
                        [
                            elDict[el.label][location][which]
                            for el in elSet.elementsByShape[elType]
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

                incrementVariableResultsArrays[elType] = results
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
            return es.EnsightPerElementVariable(
                exportJob.exportName,
                exportJob.dimensions,
                partsDict,
                self.ensightElementTypeMappings,
            )
        else:
            return None

    def _createEnsightGeometryFromModel(self, nodes, nodeSets, elements, elSets):
        partList = []
        partNumber = 1

        for elSet in elSets.values():
            elSetPart = es.EnsightUnstructuredPart(
                elSet.name,
                partNumber,
                elSet._getEnsightCompatibleElements(),
                elSet._getEnsightCompatibleReducedNodeCoords(),
                list(elSet._getEnsightCompatibleReducedNodes().keys()),
                self.ensightElementTypeMappings,
            )
            elSet.ensightPartID = partNumber
            partList.append(elSetPart)
            partNumber += 1

        geometry = es.EnsightGeometry("geometry", "-", "-", partList, "given", "given")
        return geometry
