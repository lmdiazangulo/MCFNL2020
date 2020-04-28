import numpy as np
import copy
import math

L = 0 # Lower
U = 1 # Upper

class Mesh:

    def __init__(self, coordinates, elements, grid):
        self.coordinates = coordinates
        self.elements = elements
        
        if "elemId" in grid:
            box = self.elemIdToBox(grid["elemId"])
        else:
            raise ValueError("Grid data must contain \"elemId\" or \"box\".")

        Lx = abs(box[U] - box[L])
<<<<<<< HEAD
        dx = grid["steps"]
        self.pos =  np.linspace(box[L], box[U], num=int(Lx/dx+1), endpoint=True)
=======
        if "type" not in grid or grid["type"] == "uniform":
            dx = grid["steps"]
            self.pos =  np.linspace(box[L], box[U], num=int(Lx/dx+1), endpoint=True)
        elif grid["type"] == "custom":
            self.pos = np.array([self.coordinates[el] for el in self.elements[grid["elemId"]]])
        else:
            raise ValueError("unknown grid type")

>>>>>>> 1073cfed35496dea875dc9aad8caa4d441f5aa22
        self.bounds = grid["bounds"]

    def steps(self):
        return np.array([self.pos[i+1]-self.pos[i] for i in range(len(self.pos)-1)])  

    def hsteps(self):
        steps = self.steps()
        return np.array([0.5*(steps[i]+steps[i+1]) for i in range(len(steps)-1)])

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