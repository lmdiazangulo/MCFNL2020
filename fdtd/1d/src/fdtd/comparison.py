# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 00:53:52 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""

class AnalyticComp:

    def __init__(self, mesh, probe):
        

        ids = probe["indices"]
        gridE = mesh.pos[ids[0]:ids[1]]
        probeTime = probe["time"][:]
        values    = probe["values"][:]

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
    def AnalyticalSol(self, gridE, probeTime):
        sol=[]
        A=1
        H0=1
        sig=1
        
        for time in probeTime:
            #return np.exp( - ((x-delay)**2 / (2*spread**2)) )
            #Ey=np.real([E0*np.exp()])
            sol.append(A*np.exp(-((gridE-c*time)/(sig*np.sqrt(2)))**2))
        return sol

    def L2Error(self):
        err=[]
        solReal=self.AnalyticalSol(gridE, probeTime)

        for i in range(0,len(probeTime),1):
            err.append(np.sum((values[i][:]-solReal[i][:])**2)/len(solNum[i][:]))
        return err
    
    def PrintTogether(self, mesh,probe):
        sol=self.AnalyticalSol(gridE, probeTime)
        Animator(mesh, solver.getProbes()[0], analytical=sol)
        return        