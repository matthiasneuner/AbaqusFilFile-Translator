"""
Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import sys
import os
import numpy as np
import math
import utils.exportEngine as eE
from utils.inputfileparser import parseInputFile, printKeywords
import time

FIL_WORDSIZE = 8
FIL_CHUNKSIZE = 513 * FIL_WORDSIZE
FIL_BATCHSIZE = FIL_CHUNKSIZE * 4096 * 32  # = ~ 538 MByte  ... size in BYTES


def getCurrentFileSize(
    fn,
):
    fileStat = os.stat(fn)
    fileSize = fileStat.st_size
    return fileSize


def getCurrentMaxIdxEnd(fn, fileIdx):
    fileRemainder = fileSize - fileIdx  # remaining file size in BYTES
    idxEnd = fileIdx + (FIL_BATCHSIZE if fileRemainder >= FIL_BATCHSIZE else fileRemainder)  # get end index
    # in case we are operating on an unfinished file and 'catch' an unfinished chunk
    idxEnd -= idxEnd % FIL_CHUNKSIZE
    return idxEnd


def getWords(fn, fileIdx, idxEnd):
    fnMap = np.memmap(
        fn,
        dtype="b",
        mode="r",
    )
    batchChunk = np.copy(fnMap[fileIdx:idxEnd])  # get chunk of file
    words = batchChunk.reshape(-1, FIL_CHUNKSIZE)  # get words
    words = words[:, 4:-4]
    words = words.reshape(-1, 8)
    return words


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] == "--help":
        print("Usage: filConverter.py  FILFILE.fil  EXPORTDEFINITION.inp")
        print("")
        print("Available Keywords:")
        print("")
        printKeywords()
        exit()

    fn = sys.argv[1]
    jobFile = sys.argv[2]

    exportJobs = parseInputFile(jobFile)

    exportName = "".join(fn.split("/")[-1].split(".")[-2])
    lockFile = fn.split(".")[0] + ".lck"
    print("{:<20}{:>60}".format("opening file", fn))
    print("*" * 80)

    print(exportJobs)

    exportEngine = eE.ExportEngine(exportJobs, exportName)

    fileSize = getCurrentFileSize(fn)
    numberOfBatchSteps = math.ceil(fileSize / FIL_BATCHSIZE)

    print("file has a size of {:} bytes".format(fileSize))
    print("file will be processed in {:} steps".format(numberOfBatchSteps))

    fileIdx = 0
    wordIdx = 0

    parseFile = True
    while parseFile:
        try:
            fileSize = getCurrentFileSize(
                fn,
            )

            if fileIdx < fileSize:
                idxEnd = getCurrentMaxIdxEnd(fn, fileIdx)
                words = getWords(fn, fileIdx, idxEnd)

                while wordIdx < len(words):
                    recordLength = eE.filInt(words[wordIdx])[0]
                    if recordLength <= 2:
                        print("found a record with 0 length content, possible an aborted Abaqus analysis")
                        if os.path.exists(lockFile):
                            print("found .lck file, waiting for new result .fil data")
                            time.sleep(5)
                            fileSize = getCurrentFileSize(
                                fn,
                            )
                            idxEnd = getCurrentMaxIdxEnd(fn, fileIdx)
                            words = getWords(fn, fileIdx, idxEnd)
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

                                fileSize = getCurrentFileSize(
                                    fn,
                                )
                                idxEnd = getCurrentMaxIdxEnd(fn, fileIdx)
                                words = getWords(fn, fileIdx, idxEnd)

                                continue
                            else:
                                parseFile = False
                                break

                        fileIdx += bytesProgressedInCurrentBatch  # move to beginning of the current 512 word block in the batchChunk and restart with a new bathChunk
                        wordIdx = wordIdx % 512  # of course, restart at the present index
                        break

                    recordType = eE.filInt(words[wordIdx + 1])[0]
                    recordContent = words[wordIdx + 2 : wordIdx + recordLength]
                    success = exportEngine.computeRecord(recordLength, recordType, recordContent)
                    wordIdx += recordLength

                # clean finish of a batchChunk
                if wordIdx == len(words):
                    wordIdx = 0
                    fileIdx = idxEnd
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

    print("*" * 80)
    print("Summary of .fil file:")
    print("{:<60}{:>20}".format("nodes:", len(exportEngine.allNodes)))
    print("{:<60}{:>20}".format("elements:", len(exportEngine.allElements)))
    print("{:<60}{:>20}".format("element sets:", len(exportEngine.elSets)))
    for setName, elSet in exportEngine.elSets.items():
        for elType, elements in elSet.elements.items():
            print("{:<4}{:<46}{:10}{:>11} elements".format(" ", setName, elType, len(elements)))
    print("{:<60}{:>20}".format("increments:", exportEngine.nIncrements))
