import numpy as np
import copy
import math

L = 0 # Lower
U = 1 # Upper

X = 0
Y = 1

class Mesh:

    def __init__(self, coordinates, elements, grid):
        self.coordinates = coordinates
        self.elements = elements
        
        if "elemId" in grid:
            box = self.elemIdToBox(grid["elemId"])
        else:
            raise ValueError("Grid data must contain \"elemId\" or \"box\".")

        Lxy = box[U] - box[L]
        dxy = grid["steps"]
        self.pos = np.array([
            np.linspace(box[L][X], box[U][X], num=Lxy[X]/dxy[X]+1), \
            np.linspace(box[L][Y], box[U][Y], num=Lxy[Y]/dxy[Y]+1)])

        self.bounds = grid["bounds"]

    def steps(self):
        return self.pos[1]-self.pos[0]  # Assumes everything is equispaced.

    def elemIdToBox(self, id):

        return ( np.array(self.coordinates[ self.elements[id][0] ]), \
                 np.array(self.coordinates[ self.elements[id][1] ]) )

    def toIds(self, coords):
        if type(coords) != tuple and type(coords) != list:
            coords = np.array(coords)
        ids = np.empty((0,2), int)
        nearest = np.zeros((1,2))
        for coord in coords:
            for x in [X,Y]:
                nearest[0,x] = (np.abs(self.pos[x] - coord[x])).argmin()
            ids = np.vstack((ids, nearest))
        return ids