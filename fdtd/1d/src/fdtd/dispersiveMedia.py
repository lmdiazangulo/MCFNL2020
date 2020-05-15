import numpy as np
import copy
import math
import scipy.constants as sp

class DispersiveMedia:
    def __init__(self,mesh,layer):
        #Copy outside since it is needed in every loop. Speed up the code
        self._media=copy.deepcopy(layer)
        self.ap=np.empty(len(self._media["ap"]),dtype=complex)
        self.cp=np.empty(len(self._media["cp"]),dtype=complex)
        self.epsilon=self._media['permittivity'] 
        self.mu=self._media['permeability'] 
        self.pos = self._media['startPosition']
        self.width = self._media['width']
        for i in range(0,len(self._media["ap"])):
            self.ap[i]=complex(self._media["ap"][i])
            self.cp[i]=complex(self._media["cp"][i])
        

        #Change units
        if (self._media['unitsFreq']!="Hz"):
            if (self._media['unitsFreq']=="eV"):
                self.ap = self.ap * sp.e/sp.h
                self.cp = self.cp * sp.e/sp.h

            elif (self._media['unitsFreq']=="kHz"):
                self.ap = self.ap * 10e3
                self.cp = self.cp * 10e3
            elif (self._media['unitsFreq']=="MHz"):
                self.ap = self.ap * 10e6
                self.cp = self.cp * 10e6  
            elif (self._media['unitsFreq']=="GHz"):
                self.ap = self.ap * 10e9
                self.cp = self.cp * 10e9  
            else:
                raise ValueError("Invalid frequency units. Frequency must be in multiples of Hz or eV")

        #Get layer coord
        self.indices=self.layerIndices(mesh)
        self.coords=self.layerCoords(mesh)


    def layerIndices(self, mesh):
        '''
            Gets indices of the layer
            Could take only first and last index, but it's
            better if we have all of them
        '''
        dx=mesh.steps()

        #get index of start and end position of layer in the mesh. Check a small interval of [pos-dx/2,pos+dx/2]
        indexStart = np.where((mesh.pos>=(self.pos-dx/2.0)) & (mesh.pos<(self.pos+dx/2.0)) )[0]        
        indexEnd = np.where((mesh.pos>=((self.pos+self.width)-dx/2.0)) & (mesh.pos<((self.pos+self.width)+dx/2.0)) )[0]        
 
        
        return np.arange(indexStart,indexEnd+1,1)
    
    def layerCoords(self,mesh):
        '''
            Gets coordinates of the layer. Copy from mesh to 
            avoid rounding problems when comparing
        '''
        dx=mesh.steps()
        return mesh.pos[self.indices]
    

    def _calcDispersionVar(self,dt,ap,cp):
        '''
            Computes kp and bp needed to calculate E and Jp in dispersive media
        '''
        dem=(1-ap*dt/2.0)
        kp = (1+ap*dt/2.0)/dem
        bp =  sp.epsilon_0*cp*dt/dem
        return kp,bp

