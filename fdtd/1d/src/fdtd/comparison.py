# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 00:53:52 2020

@author: María Pedrosa Bustos (mpedrosab@gmail.com)
"""
import scipy.constants as sp
import numpy as np
from fdtd.solver import Solver
from fdtd.viewer import Animator
import copy
import matplotlib.pyplot as plt

class AnalyticComp:

    def __init__(self, mesh, probe, initialCond):
        

        ids = probe["indices"]
        self.gridE = mesh.pos[ids[0]:ids[1]]
        self.probeTime = probe["time"][:]
        self.values    = probe["values"][:]
        self._initialCond = copy.deepcopy(initialCond)
        self.mesh=mesh
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
                # Source terms
        sol=[]
        for initial in self._initialCond:
            if (initial["type"] == "gaussian"):   
                for time in probeTime:
                    #return np.exp( - ((x-delay)**2 / (2*spread**2)) )
                    Ereal_right=Solver.movingGaussian(gridE,time, \
                        sp.speed_of_light,initial["peakPosition"],\
                        initial["gaussianAmplitude"]/2, \
                        initial["gaussianSpread"])
                    Ereal_left=Solver.movingGaussian(gridE,time, \
                        -sp.speed_of_light,initial["peakPosition"],\
                        initial["gaussianAmplitude"]/2, \
                        initial["gaussianSpread"])
                    sol.append(Ereal_right+Ereal_left)
            else:
                raise ValueError("Invalid source type: " + initial["type"])
           # plt.plot(gridE,sol[9])
        sol=[sol[0]]+sol[:-1]            #Repeat the initial because it is repeated in Solver
        return sol

    def L2Error(self, solReal=None):

        err=[]
        if solReal is None:
            solReal=self.AnalyticalSol(self.gridE, self.probeTime)
        for i in range(0,len(self.probeTime),1):
            err.append(np.sum((self.values[i][:]-solReal[i][:])**2)/len(self.values[i][:]))
        return err
    
    def AnimatorTogether(self,probe,SolReal=None):
        if SolReal is None:
            SolReal=self.AnalyticalSol(self.gridE, self.probeTime)
        Animator( self.mesh, probe, analytical=SolReal)
        plt.figure()
        #plt.plot(self.gridE,probe["values"][50],label='Numerical')
        #plt.plot(self.gridE,SolReal[50],label='Analytic')
        #plt.legend()
        return 

    def PrintErr(self,probe,err=None):
        if err is None:
            err=self.L2Error()
        plt.figure()
        plt.plot(self.probeTime, err)
        plt.xlabel('t')
        plt.ylabel('$L^2 error$')
        return             