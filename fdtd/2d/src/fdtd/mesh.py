from __future__ import division

import numpy as np
import copy
import math


X = 0 # Cartesian indices
Y = 1

L = 0 # Lower
U = 1 # Upper
class Mesh:

    class Bound:
        def __init__(self, ids=None):
            self.ids = ids
            
        def orientation(self):
            for xy in range(2):
                x = xy
                y = (xy + 1) % 2
                if self.ids[L][x] == self.ids[U][x] and \
                   self.ids[L][y] != self.ids[U][y]:
                    return y
            raise ValueError("Error getting orientation")

        def lu(self):
            if (self.ids[L][X], self.ids[U][X]) == (0, 0) or \
               (self.ids[L][Y], self.ids[U][Y]) == (0, 0) :
                return L
            elif (self.ids[L][X], self.ids[U][X]) == (-1, -1) or \
                 (self.ids[L][Y], self.ids[U][Y]) == (-1, -1):
                return U
            else:
                raise ValueError("Error getting lower/upper bound")

        def arrayIdx(self, lu, xy):
            if lu == L:
                return self.ids[lu][xy]
            else:
                if self.ids[U][xy] == -1:
                    return None
                else:
                    return self.ids[U][xy] + 1

        def idsAs(self, lu, xy):
            if xy is X and lu is L:
                self.ids = (np.array([ 0,  0], int), 
                            np.array([ 0, -1], int))
            if xy is X and lu is U:
                self.ids = (np.array([-1,  0], int), 
                            np.array([-1, -1], int))
            if xy is Y and lu is L:
                self.ids = (np.array([ 0,  0], int), 
                            np.array([-1,  0], int))
            if xy is Y and lu is U:
                self.ids = (np.array([ 0, -1], int), 
                            np.array([-1, -1], int))
            
            return self

    class BoundPEC(Bound):
        def __init__(self, ids=None):
            Mesh.Bound.__init__(self, ids)

    def __init__(self, coordinates, elements, grid):
        self.coordinates = coordinates
        self.elements = elements
        
        if "elemId" in grid:
            box = self.elemIdToBox(grid["elemId"])
        elif "box" in grid:
            box = grid["box"]
        else:
            raise ValueError("Grid data must contain \"elemId\" or \"box\".")

        (Lx, Ly) = abs(box[U] - box[L])
        (dx, dy) = grid["steps"]
        self.pos =  \
            (np.linspace(box[L][X], box[U][X], num=Lx/dx+1, endpoint=True),
             np.linspace(box[L][Y], box[U][Y], num=Ly/dy+1, endpoint=True) )

        self.bounds = []
        if "bounds" in grid:
            for xy in range(len(grid["bounds"])):
                for lu in range(len(grid["bounds"][xy])):
                    if grid["bounds"][xy][lu] == "pec":
                        self.bounds.append(Mesh.BoundPEC().idsAs(lu, xy))
        
                
    def steps(self):
        return (self.pos[X][1]-self.pos[X][0], self.pos[Y][1]-self.pos[Y][0])


    def origin(self):
        return (self.pos[X][0], self.pos[Y][0])

    def elemIdToBox(self, id):
        return ( np.array(self.coordinates[ self.elements[id][0] ]), \
                 np.array(self.coordinates[ self.elements[id][1] ]) )

    def toIdx(self, coords):
        if type(coords) != tuple and type(coords) != list:
            coords = [coords]
        idx = np.empty((0,2), int)
        for coord in coords:
            nearest = ((np.abs(self.pos[X] - coord[X])).argmin(),
                       (np.abs(self.pos[Y] - coord[Y])).argmin())
            idx = np.vstack((idx, np.asarray(nearest).astype(int)))
        return idx

    def snap(self, coords):
        res = []
        for coord in coords:
            id = self.toIdx(coord)[0]
            res.append( np.array([ self.pos[X][id[X]], self.pos[Y][id[Y]] ]) )
        return res