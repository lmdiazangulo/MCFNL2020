import math
import numpy as np
import scipy.constants as sp
import copy
import time

X = 0 # Cartesian indices
Y = 1

L = 0 # Lower
U = 1 # Upper

def gaussian(x, delay, spread):
    return 4*np.exp( - ((x-delay)**2 / (2*spread**2)) )

def sinusoidal(t,amplitude,omega):
    return amplitude*np.sin(omega*t)

def subsId(id):
    if id is None:
        return -1
    else:
        return id-1

def isinring(x, y, rin, rout, cent):
    r=np.sqrt((float(x)-float(cent[0]))**2.0+(float(y)-float(cent[1]))**2)
    return (r>=rin and r<rout)

def isindiamond(x, y, d, cent):
    return (abs(float(x)-float(cent[0])) + abs(float(y)- float(cent[1])))<d

def isinrectangle(x, y, h, b, cent):
    return (abs(x-cent[0])<=b/2 and abs(y-cent[1])<=h/2)



class Solver:
    
    class Fields: 
        def __init__(self, ex, ey, hz):
            self.ex = ex
            self.ey = ey
            self.hz = hz
        
        def get(self):
            return (self.ex, self.ey, self.hz)

    __timeStepPrint = 5000

    def __init__(self, mesh, options, probes, sources, material):
        self.options = options
        self._mesh = copy.deepcopy(mesh)
        self._material = copy.deepcopy(material)

        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIdx(box)
            Nxy = abs(ids[Y] - ids[X])
            p["mesh"] = {"origin": box[L], "steps": abs(box[U]-box[L]) / Nxy}
            p["indices"] = ids
            p["time"]   = [0.0]
            p["values"] = [np.zeros((Nxy[X], Nxy[Y]))]

        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIdx(box)
            source["index"] = ids
            

        self.old = self.Fields(
            ex = np.zeros( (mesh.pos[X].size-1, mesh.pos[Y].size  ) ),
            ey = np.zeros( (mesh.pos[X].size,   mesh.pos[Y].size-1) ),
            hz = np.zeros( (mesh.pos[X].size-1, mesh.pos[Y].size-1) ) )
        
        # definimos la matriz de propiedades
        self.mate = [(np.ones((mesh.pos[X].size, mesh.pos[Y].size))),
            (np.ones((mesh.pos[X].size, mesh.pos[Y].size)))]
        # determinamos las coordenadas del centro
        cent=self._mesh.elemIdToBox(material["elemId"])
        idc=mesh.toIdx(cent)
        for i in range(mesh.pos[X].size):
            for j in range(mesh.pos[Y].size):
                if "shape" in material:
                    if material["shape"] == "ring":
                        # determinamos los valores Rin y Rout
                        Ri=self._mesh.toIdl(material["Rin"])
                        Ro=self._mesh.toIdl(material["Rout"])
                        if isinring(float(i), float(j), Ri, Ro, idc[0]):
                            self.mate[0][i,j]=material["mu"]
                            self.mate[1][i,j]=material["epsilon"]
                    elif material["shape"] == "rectangle":
                        Base = self._mesh.toIdl(material["base"])
                        Height = self._mesh.toIdl(material["height"])
                        if isinrectangle(float(i), float(j), Height, Base, idc[0]):
                            self.mate[0][i,j]=material["mu"]
                            self.mate[1][i,j]=material["epsilon"]
                    elif material["shape"] == "diamond":
                        diag = self._mesh.toIdl(material["diagonal"])
                        if isindiamond(float(i), float(j), diag, idc[0]):
                            self.mate[0][i,j]=material["mu"]
                            self.mate[1][i,j]=material["epsilon"]



    def getproperties(self):
        return self.mate

    def getmaterial(self):
       return self._material
    
    def _dt(self):
        return self.options["cfl"] * min(self._mesh.steps()) / math.sqrt(2.0)  

    def timeStep(self):
        return self._dt() / sp.speed_of_light

    def getProbes(self):
        res = self._probes
        return res

    # ======================= UPDATE E =============================
    def _updateE(self, t, dt, overFields = None):
        eNew = (np.zeros( self.old.ex.shape ),
                np.zeros( self.old.ey.shape ) )
        (ex, ey, h) = self.old.get()
        e = (ex, ey)

        (dX, dY) = self._mesh.steps()
        A = dX * dY
        prop=self.getproperties()
        material = self.getmaterial()

        eNew[X][:,1:-1] = e[X][:,1:-1] + dt/A*dX * (h[:,1:] - h[:,:-1])
        eNew[Y][1:-1,:] = e[Y][1:-1,:] - dt/A*dY * (h[1:,:] - h[:-1,:])

        # Boundary conditions
        for bound in self._mesh.bounds:
            xy = bound.orientation()
            (lx, ux) = (bound.arrayIdx(L,X), \
                        bound.arrayIdx(U,X))
            (ly, uy) = (bound.arrayIdx(L,Y), \
                        bound.arrayIdx(U,Y))
            if isinstance(bound, self._mesh.BoundPEC):
                eNew[xy][lx:ux,ly:uy] = 0.0
            else:
                raise ValueError("Unrecognized boundary type")
             # ----- PEC (?) ---------
            if material["boundary"] == "pec":
                for i in range(prop[0][1:,1].size):
                    for j in range(prop[0][1,1:].size):
                        if prop[0][i,j] != 1:
                            eNew[X][i,j] = 0
                            eNew[Y][i,j] = 0

        
        # Subgridding and updating
        e[X][:] = eNew[X][:]
        e[Y][:] = eNew[Y][:]  

    # ======================= UPDATE H =============================
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.hz.shape )
        (ex, ey, h) = self.old.get()
        
        (dX, dY) = self._mesh.steps()
        A = dX * dY
        prop=self.getproperties()
        #material = self.getmaterial()
              
        hNew[:,:] = h[:,:] \
                     - (dt/A * dY / prop[0][1:,1:]) * ey[1:,  :] \
                     + (dt/A * dX / prop[0][1:,1:]) * ex[ :, 1:] \
                     + (dt/A * dY / prop[0][1:,1:]) * ey[:-1,   :] \
                     - (dt/A * dX / prop[0][1:,1:]) * ex[  :, :-1]
        
        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    c0 = sp.speed_of_light
                    delay  = c0 * magnitude["gaussianDelay"]
                    spread = c0 * magnitude["gaussianSpread"]
                    id = source["index"]
                    hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                        gaussian(t, delay, spread)*dt
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])

            elif source["type"] == "planewave":
                magnitude = source["magnitude"]
                if magnitude["type"] == "sinusoidal":
                    amplitude = magnitude["sinAmplitude"]
                    omega = magnitude["sinFreq"]
                    id = source["index"]
                    hNew[id[L][X]:id[U][X], id[L][Y]:id[U][Y]] += \
                    amplitude*np.sin(omega*t)
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])
                    
            else:
                raise ValueError("Invalid source type: " + source["type"])

        h[:] = hNew[:]
            
    def _updateProbes(self, t):
        for p in self._probes:
            dimensionalTime = t/sp.speed_of_light
            writeStep = "samplingPeriod" in p \
                and (dimensionalTime/p["samplingPeriod"] >= len(p["time"]))
            writeStep = writeStep or "samplingPeriod" not in p
            if writeStep:
                p["time"].append(dimensionalTime)
                idx = p["indices"]
                values = np.zeros(tuple(idx[U]-idx[L]))
                values[:,:] = \
                    self.old.hz[ idx[L][X]:idx[U][X], idx[L][Y]:idx[U][Y] ]
                p["values"].append(values)

    def solve(self, dimensionalFinalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        numberOfTimeSteps = \
            int(dimensionalFinalTime * sp.speed_of_light / dt)
        for n in range(numberOfTimeSteps):
        
            self._updateE(t, dt, self.old)
            t += dt/2.0

            self._updateH(t, dt)
            t += dt/2.0

            self._updateProbes(t)
    
            if n % self.__timeStepPrint == 0 or n+1 == numberOfTimeSteps:
                remaining = (time.time() - tic) * \
                    (numberOfTimeSteps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("    Step: %6d of %6d. Remaining: %2.0f:%02.0f"% (n, \
                    numberOfTimeSteps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))

    @staticmethod
    def movingGaussian(x,y,t,c,center,A,spread):
        return A*np.exp(-(((x-center)-c*t)**2 /(2*spread**2)))