"""
Created on Fri Oct 28 13:43:53 2016

Copyright (C) 2019 Matthias Neuner <matthias.neuner@uibk.ac.at>

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""

import numpy as np
from collections import defaultdict


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
