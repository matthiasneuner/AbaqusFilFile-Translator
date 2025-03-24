"""
Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np

# a word in a .fil file has a size of 8 bytes
FIL_WORDSIZE = 8
# a chunk in a .fil file consists of 513 words
FIL_CHUNKSIZE = 513 * FIL_WORDSIZE
# .fil files may become huge. We do not load them at once, but
# but rather we split them into multiple batch sizes
FIL_BATCHSIZE = FIL_CHUNKSIZE * 4096 * 32  # = ~ 538 MByte  ... size in BYTES


def getCurrentMaxIdxEnd(fn: str, fileIdx: str, fileSize: int):
    """Determine the maximum index in the .fil file depending on the current position.
    It may be a complete chunk, or less if we are already near the end of the .fil file.

    Parameters
    ----------
    fn
        The .fil file name.
    fileIdx
        The current file index.

    Returns
    -------
    type
        The maximum allowed index in the file.
    """

    fileRemainder = fileSize - fileIdx  # remaining file size in BYTES
    idxEnd = fileIdx + (FIL_BATCHSIZE if fileRemainder >= FIL_BATCHSIZE else fileRemainder)  # get end index
    # in case we are operating on an unfinished file and 'catch' an unfinished chunk
    idxEnd -= idxEnd % FIL_CHUNKSIZE
    return idxEnd


def getFilFileWords(fn, fileIdx, idxEnd):
    """Get readable words between the fileIdx and idxEnd.

    Parameters
    ----------
    fileIdx
        The current file index.
    idxEnd
        The end index.

    Returns
    -------
    type
        The words.
    """
    fnMap = np.memmap(
        fn,
        dtype="b",
        mode="r",
    )
    batchChunk = np.copy(fnMap[fileIdx:idxEnd])  # get chunk of file
    words = batchChunk.reshape(-1, FIL_CHUNKSIZE)  # get words
    # strip unused bytes. Probably they contain checksums,
    # so we may leverage that feature in a future version.
    words = words[:, 4:-4]
    words = words.reshape(-1, 8)
    return words
