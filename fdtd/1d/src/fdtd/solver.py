
import math
import numpy as np
import scipy.constants as sp
import copy
import time
#import matplotlib.pyplot as plt

L = 0 # Lower
U = 1 # Upper

class Fields: 
    def __init__(self, e, h):
        self.e = e
        self.h = h

    def get(self):
        return (self.e, self.h)

class ComplexField: 

    '''
        For dispersive media
    '''
    def __init__(self, Jp_old,eDispersive,hDispersive):
        self.Jp_old= Jp_old
        self.eDispersive = eDispersive      #To save wave in dispersive and free space to compare
        self.hDispersive = hDispersive      #To save wave in dispersive and free space to compare
    
    def get(self):
        return (self.Jp_old,self.eDispersive,self.hDispersive)        #Define for the whole mesh but only use it
                                    # where the dispersive layers is placed


class Solver:
    
    _timeStepPrint = 100

    def __init__(self, mesh, options, probes,sources, initialCond=[{"type": "none"}],dispLayer=None):
        self.options = options
        
        self._mesh = copy.deepcopy(mesh)
        self._initialCond = copy.deepcopy(initialCond)
        self._dispLayer = None          #need to check along the code if there is a layer
 
        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIds(box)
            Nx = abs(ids)

            p["mesh"] = {"origin": box[L], "steps": abs(box[U]-box[L]) / Nx}
            p["indices"] = ids
            p["time"]   = [0.0]
            
            for initial in self._initialCond:
                if initial["type"]=="none":
                    values=np.zeros( mesh.pos.size )
                    p["values"] = [np.zeros((1,Nx[1]))]          
                elif ( initial["type"] == "gaussian"):
                    position=self._mesh.pos
                    #print(source["index"])  #Lugar del pico
                    values=Solver.movingGaussian(position, 0, \
                       sp.speed_of_light,initial["peakPosition"],\
                       initial["gaussianAmplitude"], \
                       initial["gaussianSpread"] )  
                    p["values"]= [values[ids[0]:ids[1]]]
                        #plt.plot(position,eNew)
                else:
                    raise ValueError(\
                    "Invalid initial condition type: " + initial["type"] )


        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIds(box)
            source["index"] = ids

        self.old = Fields(e = values.copy(),
                          h = np.zeros( mesh.pos.size-1 ) )

       #Get dispersivee media properties in case it is defined
        if dispLayer is not None:

            self._dispLayer=copy.deepcopy(dispLayer)
            self._epsilon=self._dispLayer.epsilon
            self._mu=self._dispLayer.mu

            #Copy outside since it is needed in every loop. Speed up the code
            self._ap=self._dispLayer.ap
            self._cp=self._dispLayer.cp 
             #Changed?
            self.oldDispersive = ComplexField(Jp_old = np.zeros(( self._dispLayer.coords.size, len(self._dispLayer.ap))), eDispersive= values.copy(), hDispersive = np.zeros( mesh.pos.size-1 ) )      #Takes the size of the layer, not the full grid
            self._layerIndices = self._dispLayer.indices
            #Save values as if there is no layer to layer compare the results 
            p["valuesFree"] = p["values"].copy()

    def solve(self, finalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        print ('time step:',dt)
        numberOfTimeSteps = int(finalTime / dt)

        if self._dispLayer is not None:
            self._kp,self._bp=self._dispLayer ._calcDispersionVar(dt,self._dispLayer.ap,self._dispLayer.cp)

        for n in range(numberOfTimeSteps):
            self._updateE(t, dt)
            t += dt/2.0
            self._updateH(t, dt)
            t += dt/2.0
            self._updateProbes(t)
    
            if n % self._timeStepPrint == 0 or n+1 == numberOfTimeSteps:
                remaining = (time.time() - tic) * \
                    (numberOfTimeSteps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("    Step: %6d of %6d. Remaining: %2.0f:%02.0f"% (n, \
                    numberOfTimeSteps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))


    def _dt(self):
        print ('dt for cfl condition: ', self.options["cfl"] * self._mesh.steps() / sp.speed_of_light)
        input("Press Enter to continue...")
        return self.options["cfl"] * self._mesh.steps() / sp.speed_of_light  

    def timeStep(self):
        return self._dt()

    def getProbes(self):
        res = self._probes
        return res

    def _updateE(self, t, dt):
        (e, h) = self.old.get()
        eNew = np.zeros( self.old.e.shape )


        cE = dt / sp.epsilon_0 / self._mesh.steps()
        eNew[1:-1] = e[1:-1] +    cE * (h[1:] - h[:-1])

        if self._dispLayer is not None:
            (Jp_old, eDisp, hDisp ) = self.oldDispersive.get()
            JpNew = np.zeros( Jp_old.shape )
            eNewDisp = np.zeros( eDisp.shape )
            eNewDisp[1:-1]= eDisp[1:-1] +    cE * (hDisp[1:] - hDisp[:-1])

            cE2 = (2*dt/((2*self._epsilon*sp.epsilon_0+np.sum(2*np.real(self._bp)))*self._mesh.steps()))
            #Term multiplying e[1:-1] is 1 since conductivity=0
            #Need to add an extra term to indices for dh/dx
            indH= np.concatenate(([self._layerIndices[0]-1],self._layerIndices))
            eNewDisp[self._layerIndices] = eDisp[self._layerIndices] + cE2 * ((hDisp[indH[1:]] \
               - hDisp[indH[:-1]])-np.real(np.sum((1+self._kp)*Jp_old[:,:],1)))


        #add mur conditions to layer
        # Boundary conditions
        #WARNING: assumes dispersive layer is not at the boundary
        for lu in range(2):
            if lu == 0:
                pos = 0
            else:
                pos = -1
            if self._mesh.bounds[lu] == "pec":
                eNew[pos] = 0.0
                if self._dispLayer is not None:
                    eNewDisp[pos] = 0.0
            elif self._mesh.bounds[lu] == 'pmc':
                eNew[pos] = e[pos] + 2*cE*(h[pos] if pos == 0 else -h[pos])
                if self._dispLayer is not None:
                    eNewDisp[pos] = eDisp[pos] + 2*cE*(hDisp[pos] if pos == 0 else -h[pos])
            elif self._mesh.bounds[lu] == 'mur':
                if pos == 0:
                    eNew[0] =  e[ 1]+(sp.speed_of_light*dt-self._mesh.steps())* (eNew[ 1]-e[ 0]) / (sp.speed_of_light*dt+self._mesh.steps())
                    if self._dispLayer is not None:
                        eNewDisp[0] = eDisp[ 1]+(sp.speed_of_light*dt-self._mesh.steps())* (eNewDisp[ 1]-eDisp[ 0]) / (sp.speed_of_light*dt+self._mesh.steps())
                else:
                    eNew[-1] = e[-2]+(sp.speed_of_light*dt-self._mesh.steps())* (eNew[-2]-e[ 1]) / (sp.speed_of_light*dt+self._mesh.steps())
                    if self._dispLayer is not None:
                        eNewDisp[-1] = eDisp[-2]+(sp.speed_of_light*dt-self._mesh.steps())* (eNewDisp[-2]-e[ 1]) / (sp.speed_of_light*dt+self._mesh.steps())
 
            else:
                raise ValueError("Unrecognized boundary type")

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    eNew[source["index"]] += Solver._gaussian(t, \
                        magnitude["gaussianDelay"], \
                        magnitude["gaussianSpread"] )  

                if self._dispLayer is not None:
                    eNewDisp[source["index"]] += Solver._gaussian(t, \
                            magnitude["gaussianDelay"], \
                            magnitude["gaussianSpread"] )                         
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])

            elif source["type"] == "none":
                continue
            else:
                raise ValueError("Invalid source type: " + source["type"])
        
        if self._dispLayer is not None:

            #Get new JpNew
            for i in range(0,np.shape(Jp_old)[1]):
                
                #Changed?
              # JpNew[1:-1,i] = self._kp[i]*Jp_old[1:-1,i]+ self._bp[i] * (eNew[self._layerIndices[1:-1]] - e[self._layerIndices[1:-1]]) / dt
               JpNew[:,i] = self._kp[i]*Jp_old[:,i]+ self._bp[i] * (eNewDisp[self._layerIndices[:]] - eDisp[self._layerIndices[:]]) / dt
            Jp_old[:,:] = JpNew[:,:]
            eDisp[:] = eNewDisp[:]
        e[:] = eNew[:]
        
        
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.h.shape )
        (e, h) = self.old.get()
        cH = dt / sp.mu_0 / self._mesh.steps()
          #Get also field for free particle      
        hNew[:] = h[:] + cH * (e[1:] - e[:-1])
        h[:] = hNew[:]
                
        if self._dispLayer is not None:
            hNewDisp = np.zeros(self.oldDispersive.hDispersive.shape )
            (J_old,eDisp, hDisp) = self.oldDispersive.get()
            cH2 = dt / (self._mu*sp.mu_0) / self._mesh.steps()
            hNewDisp[:] = hDisp[:] + cH2 * (eDisp[1:] - eDisp[:-1])
            hDisp[:] = hNewDisp[:]
            
    def _updateProbes(self, t):
        for p in self._probes:
            if "samplingPeriod" not in p or \
               "samplingPeriod" in p and \
               (t/p["samplingPeriod"] >= len(p["time"])):
                p["time"].append(t)
                ids = p["indices"]
                values = np.zeros(ids[U]-ids[L])
                values[:] = self.old.e[ ids[0]:ids[1] ]
                 #Save values as if there is no layer and layer compare the results
               
                if self._dispLayer is not None:
                    valuesDisp = np.zeros(ids[U]-ids[L])
                    valuesDisp[:] = self.oldDispersive.eDispersive[ ids[0]:ids[1] ]
                    p["values"].append(valuesDisp)
                    p["valuesFree"].append(values)
                else:
                    p["values"].append(values)

    def _calcDispersionVar(self,dt,ap,cp):

        dem=(1-ap*dt/2.0)
        self._kp = (1+ap*dt/2.0)/dem
        self._bp =  sp.epsilon_0*cp*dt/dem
        return self._kp,self._bp

    @staticmethod
    def _gaussian(x, delay, spread):
        return np.exp( - ((x-delay)**2 / (2*spread**2)) )
    
    def movingGaussian(x,t,c,center,A,spread):
        return A*np.exp(-(((x-center)-c*t)**2 /(2*spread**2)))