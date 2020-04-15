# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 00:53:52 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""
import scipy.constants as sp
import numpy as np

class AnalyticComp:

    def __init__(self, mesh, probe):
        

        ids = probe["indices"]
        self.gridE = mesh.pos[ids[0]:ids[1]]
        self.probeTime = probe["time"][:]
        self.values    = probe["values"][:]

    #Get data directly from json
    '''
    self._mesh = copy.deepcopy(mesh)
       
    self._probes = copy.deepcopy(probes)
    def getProbes(self):
        res = self._probes
        return res   
    
     def _getTime(self, t):
                 t = 0.0
        dt = self._dt()
        numberOfTimeSteps = int(finalTime / dt)
        for n in range(numberOfTimeSteps):
            t += dt/2.0
            t += dt/2.0
            self._updateProbes(t)
        for p in self._probes:
            if "samplingPeriod" not in p or \
               "samplingPeriod" in p and \
               (t/p["samplingPeriod"] >= len(p["time"])):
                p["time"].append(t)   
                '''
    #Get analytical solution for the same times and positions
    def AnalyticalSol(self,gridE,probeTime):
        sol=[]
        H0=1
        E0y=1
        E0z=1
        k=1
        w=k*sp.speed_of_light 
        sig=1
        
        for time in probeTime:
            #return np.exp( - ((x-delay)**2 / (2*spread**2)) )
            E=np.real([E0y*np.exp(k*self.gridE-w*time*1j),E0z*np.exp(k*self.gridE-w*time*1j)])
            sol.append(np.sqrt(E[0,:]**2+E[1,:]**2))
        return sol

    def L2Error(self):
        err=[]
        solReal=self.AnalyticalSol(self.gridE, self.probeTime)

        for i in range(0,len(self.probeTime),1):
            err.append(np.sum((self.values[i][:]-solReal[i][:])**2)/len(self.values[i][:]))
        return err
    
    def PrintTogether(self, mesh,probe):
        sol=self.AnalyticalSol(self.gridE, self.probeTime)
        Animator(mesh, solver.getProbes()[0], analytical=sol)
        return        