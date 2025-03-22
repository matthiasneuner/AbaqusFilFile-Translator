"""
Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import argparse
import sys
import os
import numpy as np
import math
from src.exportengine import ExportEngine, filInt
from src.inputfileparser import parseInputFile, printKeywords
from src.misc import fileSizeHumanReadable, getCurrentFileSize
from src.filfileformat import FIL_BATCHSIZE, getCurrentMaxIdxEnd, getFilFileWords
import time
import textwrap

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A translator for Abaqus .fil files.")

    parser.add_argument(
        "fil",
        metavar="FILFILE.fil",
        help="The Abaqus .fil file",
        type=str,
    )
    parser.add_argument(
        "expDef",
        metavar="EXPORTDEFINITION.inp",
        help="The .inp export definition file",
        type=str,
    )
    parser.add_argument("--keywords", dest="kw", action="store_true", help="print keywords")
    parser.add_argument("--verbose", dest="verbose", action="store_true", help="print verbose output")
    args = parser.parse_args()

    if args.kw:
        printKeywords()
        exit(0)

    fn = args.fil
    jobFile = args.expDef

    exportJobs = parseInputFile(jobFile)

    exportName = "".join(fn.split("/")[-1].split(".")[-2])
    lockFile = fn.split(".")[0] + ".lck"
    print("+" + "-" * 78 + "+")
    print("| Opening file {:<64}|".format(os.path.basename(fn)))
    print("+" + "-" * 78 + "+")

    exportEngine = ExportEngine(exportJobs, exportName, verbose=args.verbose)

    currentFileSize = getCurrentFileSize(fn)
    numberOfBatchSteps = math.ceil(currentFileSize / FIL_BATCHSIZE)

    print("file has a size of {:}".format(fileSizeHumanReadable(currentFileSize)))
    print("file will be processed in {:} batch(es)".format(numberOfBatchSteps))

    currentFileIdx = 0
    wordIdx = 0

    parseFile = True
    while parseFile:
        try:
            currentFileSize = getCurrentFileSize(
                fn,
            )

            if currentFileIdx < currentFileSize:
                idxEnd = getCurrentMaxIdxEnd(fn, currentFileIdx, currentFileSize)
                words = getFilFileWords(fn, currentFileIdx, idxEnd)

                while wordIdx < len(words):
                    recordLength = filInt(words[wordIdx])[0]
                    if recordLength <= 2:
                        print("found a record with 0 length content, possible an aborted Abaqus analysis")
                        if os.path.exists(lockFile):
                            print("found .lck file, waiting for new result .fil data")
                            time.sleep(5)
                            currentFileSize = getCurrentFileSize(
                                fn,
                            )
                            idxEnd = getCurrentMaxIdxEnd(fn, currentFileIdx)
                            words = getFilFileWords(fn, currentFileIdx, idxEnd)
                            continue
                        else:
                            parseFile = False
                            break

                    # the next record exceeds our batchChunk, so we do some trick:
                    # - set the wordIdx to the end of the so far progressed frame
                    # - move the frame to the wordIDx
                    if wordIdx + recordLength > len(words):
                        bytesProgressedInCurrentBatch = int(math.floor(wordIdx / 512)) * 513 * 8

                        if bytesProgressedInCurrentBatch == 0:  # indicator for an aborted analysis
                            print("terminated file, possible an aborted Abaqus analysis")
                            if os.path.exists(lockFile):
                                print("found .lck file, waiting for new result .fil data")
                                time.sleep(5)

                                currentFileSize = getCurrentFileSize(
                                    fn,
                                )
                                idxEnd = getCurrentMaxIdxEnd(fn, currentFileIdx)
                                words = getFilFileWords(fn, currentFileIdx, idxEnd)

                                continue
                            else:
                                parseFile = False
                                break

                        currentFileIdx += bytesProgressedInCurrentBatch  # move to beginning of the current 512 word block in the batchChunk and restart with a new bathChunk
                        wordIdx = wordIdx % 512  # of course, restart at the present index
                        break

                    # Time to process the record!
                    recordType = filInt(words[wordIdx + 1])[0]
                    recordContent = words[wordIdx + 2 : wordIdx + recordLength]
                    success = exportEngine.computeRecord(recordLength, recordType, recordContent)
                    wordIdx += recordLength

                # clean finish of a batchChunk
                if wordIdx == len(words):
                    wordIdx = 0
                    currentFileIdx = idxEnd
                del words

            else:
                if os.path.exists(lockFile):
                    print("found .lck file, waiting for new result .fil data or CTRL-C to finish...")
                    time.sleep(10)
                else:
                    break

        except KeyboardInterrupt:
            print("Interrupted by user")
            break

    exportEngine.finalize()

    print("+" + "-" * 78 + "+")
    print("| Summary of {:<66}|".format(os.path.basename(fn)))
    print("+" + "-" * 78 + "+")
    print("|{:<60}{:>18}|".format("nodes:", len(exportEngine.nodes)))
    print("|{:<60}{:>18}|".format("elements:", len(exportEngine.elements)))
    print("|{:<60}{:>18}|".format("element sets:", len(exportEngine.elSets)))
    for setName, elSet in exportEngine.elSets.items():
        for i, (elType, elements) in enumerate(elSet.elementsByShape.items()):
            print("|{:<4}{:<46}{:10}{:>9} elements|".format(" ", setName if not i else "", elType, len(elements)))
    print("|{:<60}{:>18}|".format("node sets:", len(exportEngine.nSets)))
    for setName, nSet in exportEngine.nSets.items():
        print("|{:<4}{:<46}{:10}{:>9}    nodes|".format(" ", setName, "", len(nSet.nodes)))
    print("|{:<60}{:>18}|".format("increments:", exportEngine.nIncrements))
    print("+" + "-" * 78 + "+")
    print("|{:78}|".format(" Finished"))
    print("+" + "-" * 78 + "+")
