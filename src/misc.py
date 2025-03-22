"""
Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

from collections import defaultdict
import os


class RecursiveDefaultDict(dict):
    # return defaultdict(RecursiveDefaultDict)

    def __init__(self, maxLevels):
        self.level = maxLevels
        return super().__init__()

    def __getitem__(self, key):
        if key in self:
            return self.get(key)

        if self.level > 0:
            return self.setdefault(key, RecursiveDefaultDict(self.level - 1))

        return self.setdefault(key, dict())


def sliceFromString(string: str, shift: int = 0):
    """Generate a slice from a string, which can represent a slice or an index.

    Parameters
    ----------
    string
        The string representing the slice/index.
    shift
        Potential shift.

    Returns
    -------
    type
        The slice!
    """

    if ":" in string:
        a, b = string.split(":")
        return slice(int(a) + shift, int(b) + shift)
    else:
        return slice(int(string) + shift, int(string) + 1 + shift)


def makeExtractionFunction(expression, symbol="x"):
    """make a simple f(x) expression from string"""
    return lambda x: eval(expression, globals(), {symbol: x})


def fileSizeHumanReadable(num: int, suffix: str = "B"):
    """Pretty format bytes to a human readable form.

    Parameters
    ----------
    int
        The number.
    suffix
        The suffix.

    Returns
    -------
    type
        The pretty formatted string.
    """

    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1000:
            return f"{num:3.1f} {unit}{suffix}"
        num /= 1000
    return f"{num:.1f} Y{suffix}"


def getCurrentFileSize(
    fn: str,
):
    """Determine the size of the .fil file.

    Parameters
    ----------
    fn
        The .fil file name.

    Returns
    -------
    type
        The file size.
    """

    fileStat = os.stat(fn)
    fileSize = fileStat.st_size
    return fileSize
