"""
Created on Tue Sep 22 18:56:03 2015

Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np
from os.path import dirname, join
import textwrap
import shlex


class InputSyntaxException(Exception):
    pass


dTypes = {
    int: "integer",
    float: "float",
    str: "string",
}

typeMappings = {
    "*defineElementType": (
        "assign an ensight Shape to an Abaqus Element",
        {
            "element": (str, "Abaqus (User) Element"),
            "shape": (
                str,
                "Ensight Shape, can be any of: point bar2 bar3 tria3 tria6 quad4 quad8 tetra4 tetra10 pyramid5 pyramid13 penta6 penta15 hexa8 hexa20 nsided nfaced ",
            ),
        },
    ),
    "*ignoreLastNodesForElementType": (
        "Ignore trailing nodes to be ignored (e.g, make a hexa27 to a hex20 with number=7)",
        {
            "element": (str, "Abaqus (User) Element"),
            "number": (
                int,
                "The number of nodes to be ignored",
            ),
        },
    ),
    "*ensightCaseOptions": (
        "modify Ensight export options",
        {
            "discardTime": (
                str,
                "discard Time values and replace by enumeration of time steps",
            ),
        },
    ),
    "*computeAverageOverQuadraturePoints": (
        "perform a computation on an elemental result",
        {
            "set": (str, "Abaqus element set"),
            "result": (str, "Abaqus variable identifier"),
        },
    ),
    "*UELSDVToQuadraturePoints": (
        "relate SDV data to quadrature points",
        {
            "set": (str, "Abaqus element set"),
            "destination": (str, "new name of the result"),
            "qpCount": (
                int,
                "define a periodical pattern for a repeated extraction for results at quadrature points",
            ),
            "qpDistance": (
                int,
                "define a periodical pattern: data distance between qps",
            ),
            "qpInitialOffset": (
                int,
                "define a periodical pattern: initial constant offset before qp data begins",
            ),
        },
    ),
    "*ensightPerNodeVariableJob": (
        "define an Ensight per node variable for export",
        {
            "name": (str, "export name of the variable"),
            "dimensions": (
                int,
                "(optional), 1/3/6/9 for scalar/vector/tensor/tensor asym; missing components will be zero filled",
            ),
            "timeSet": (int, "(optional), define a timeset, for 'different' timelines"),
        },
    ),
    "*ensightPerNodeVariableJobEntry": (
        "define an Ensight per node variable for an element set",
        {
            "job": (str, "The associated export job"),
            "setType": (str, "elSet or nSet, default=elSet"),
            "set": (str, "Abaqus setname"),
            "result": (str, "Abaqus variable identifier"),
            "values": (
                str,
                "(optional), define a index/slice to extract a subarray from the total result array (per Element)",
            ),
            "f(x)": (
                str,
                "(optional), apply a mathematical/array expression on the result array (per Element, slow!)",
            ),
            "fillMissingValues": (
                float,
                "(optional), fill missing nodal values with a constant values, requires specified dimensions (slow!)",
            ),
        },
    ),
    "*ensightPerElementVariableJob": (
        "define an Ensight per element variable for export",
        {
            "name": (str, "export name of the variable"),
            "dimensions": (
                int,
                "(optional), 1/3/6/9 for scalar/vector/tensor/tensor9; missing components will be zero filled",
            ),
            "timeSet": (int, "(optional), define a timeset, for 'different' timelines"),
        },
    ),
    "*ensightPerElementVariableJobEntry": (
        "define an Ensight per element variable entry for an element set",
        {
            "job": (str, "export name of the variable"),
            "set": (str, "Abaqus element set"),
            "result": (str, "Abaqus variable identifier"),
            "location": (str, "where is the result ? qps | computed "),
            "which": (
                str,
                'which one? e.g. quadrature point numbers or "average" for average computed results',
            ),
            "values": (
                str,
                "(optional), define a index/slice to extract a subarray from the total result array (per Element)",
            ),
            "f(x)": (
                str,
                "(optional), apply a mathematical/array expression on the result array (per Element, slow!)",
            ),
        },
    ),
    "*include": (
        "(optional) load extra .inp file (fragment), use relative path to current .inp",
        {"input": (str, "filename")},
    ),
}


def getMapType(kw, varName):
    kwDict = typeMappings.get(kw, (None, {}))[1]
    mType = kwDict.get(varName, [str])[0]
    return mType


def parseInputFile(fileName, currentKeyword=None, existingFileDict=None):
    """parse an Abaqus like input file to generate an dictionary with its content"""
    if not existingFileDict:
        fileDict = {key: [] for key in typeMappings.keys()}
    else:
        fileDict = existingFileDict
    keyword = currentKeyword
    with open(fileName) as f:
        for l in f:
            if l.startswith("**"):
                # line is a comment
                continue
            lexer = shlex.shlex(l.strip(), posix=True)
            lexer.whitespace_split = True
            lexer.whitespace = ","
            lineElements = [x.strip() for x in lexer]
            if not lineElements:
                pass
            elif lineElements[0].startswith("*"):
                # line is keywordline
                lastkeyword = keyword
                keyword = lineElements[0]

                if keyword not in fileDict:
                    raise InputSyntaxException(
                        "Invalid keyword {:} provided".format(keyword)
                    )

                optionAssignments = lineElements[1:]

                objectentry = {}
                objectentry[
                    "inputFile"
                ] = fileName  # save also the filename of the original inputfile!

                for ass in optionAssignments:
                    opts = ass.split("=")
                    optKey = opts[0].strip()
                    val = opts[1].strip()
                    try:
                        mType = getMapType(keyword, optKey)
                    except KeyError:
                        raise InputSyntaxException(
                            "'{:}' is not a valid option for keyword {:}".format(
                                optKey, keyword
                            )
                        )

                    try:
                        objectentry[optKey] = mType(val)
                    except:
                        raise InputSyntaxException(
                            "'{:}' is not of type {:} for option {:} in keyword {:}".format(
                                val, mType, optKey, keyword
                            )
                        )

                # special treatment for *include:
                if keyword == "*include":
                    includeFile = objectentry["input"]
                    parseInputFile(
                        join(dirname(fileName), includeFile),
                        currentKeyword=lastkeyword,
                        existingFileDict=fileDict,
                    )
                    keyword = lastkeyword
                else:
                    fileDict[keyword].append(objectentry)
    return fileDict


def printKeywords():
    """print the input file language set"""
    kwString = "    {:}    "
    kwDataString = "        {:30}{:14}"

    wrapper = textwrap.TextWrapper(width=120, replace_whitespace=False)
    for kw, (kwDoc, optiondict) in sorted(typeMappings.items()):
        wrapper.initial_indent = kwString.format(str(kw))
        wrapper.subsequent_indent = " " * len(wrapper.initial_indent)
        print(wrapper.fill(kwDoc))
        print("")

        for key in sorted(optiondict.keys()):
            optionName = key
            dType, description = optiondict[key]
            wrapper.initial_indent = kwDataString.format(str(optionName), dTypes[dType])
            wrapper.subsequent_indent = " " * len(wrapper.initial_indent)
            print(wrapper.fill(description))
        print("\n")
