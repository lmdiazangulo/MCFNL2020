import numpy as np
import copy
import math

L = 0 # Lower
U = 1 # Upper

class Mesh:

    def __init__(self, coordinates, elements, grid):
        self.coordinates = coordinates
        self.elements = elements
        
        _,pos = self.positions(grid[0])
        for g in grid[1:]:
            box, paux = self.positions(g)
            bottom = pos[pos < box[L]]
            top = pos[pos > box[U]]
            pos = np.concatenate([bottom,paux,top])

        self.pos = pos
        self.bounds = grid[0]["bounds"]

    def positions(self,grid):
        if "elemId" in grid:
            box = self.elemIdToBox(grid["elemId"])
        else:
            raise ValueError("Grid data must contain \"elemId\" or \"box\".")
        Lx = abs(box[U] - box[L])
        dx = grid["steps"]
        return  box, np.linspace(box[L], box[U], num=int(Lx/dx+1), endpoint=True)

    def steps(self):
        return np.array(self.pos[1:]-self.pos[0:-1])  

    def hsteps(self):
        return np.array(0.5*(self.steps()[0:-1]+self.steps()[1:]))

    def elemIdToBox(self, id):

        return ( self.coordinates[ self.elements[id][0] ], \
                 self.coordinates[ self.elements[id][-1] ] )

    def toIds(self, coords): # Will send source term to nearest coord
        if type(coords) != tuple and type(coords) != list:
            coords = np.array(coords)
        ids = np.empty((0,1), int)
        for coord in np.unique(coords):
            nearest = (np.abs(self.pos - coord)).argmin()
            ids = np.append(ids, nearest)
        return ids

    def snap(self, coords):
        res = []
        for coord in coords:
            id = np.array(self.toIds(coord)[0], dtype=int)
            res.append( np.array([ self.pos[id] ]) )
        return res