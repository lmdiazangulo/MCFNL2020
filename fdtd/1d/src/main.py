
import json
import argparse
import os.path
import sys

from fdtd.mesh import Mesh
from fdtd.solver import Solver
from fdtd.viewer import Animator
from fdtd.comparison import AnalyticComp
from fdtd.dispersiveMedia import DispersiveMedia
from measure.Transmittance import MeasureTransmittance, AnalyticTransmittance, PlotTransmittance
print("=== Python FDTD 1D")

'''
parser = argparse.ArgumentParser(description='Python FDTD 1D')
parser.add_argument('-i', '--input', nargs=1, type=str)
args = parser.parse_args()0
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

inputFilename = ''.join(args.input).strip()
'''
inputFilename='..\\tests\\cavity_dispersive_test_NoComplex.json'
print("--- Reading file: %s"%(inputFilename))
data = json.load(open(inputFilename))



print('--- Initializing mesh')
mesh = Mesh(data["coordinates"], data["elements"], data["grid"])
#layer = None
layer = DispersiveMedia(mesh,data["dispersiveLayers"])
#indeces = layer.layerIndices(mesh)
print('--- Initializing solver')
#layer=None
solver = Solver(mesh, data["options"], data["probes"], data["sources"],dispLayer = layer)

print('--- Solving')
solver.solve(data["options"]["finalTime"])

print('--- Visualizing')
solNum=solver.getProbes()[0]

#Measure transmittance
transmission = MeasureTransmittance(layer,solNum['time'],solNum['valuesFree'],solNum['values'])
freq, transCoef= transmission.T()
freq2, reflecCoef= transmission.R()
PlotTransmittance(freq, [transCoef,reflecCoef],labels=['T','R'])
freq2, transmittance,reflectance= transmission.TransmittanceReflect()
PlotTransmittance(freq2, [transmittance,reflectance],labels=['Transmittance','Reflectance'])

#Analytical transmittance
#layer.width = 100e-9
transmittanceReal = AnalyticTransmittance(layer)
import numpy as np
#freq  = np.linspace(1e2/(2 * np.pi), 1e10/(2 * np.pi), int(1e2+1)) * 2 * np.pi
#freq = np.linspace(1e14,12e14,1000)

transReal = transmittanceReal.T(freq)
reflecReal = transmittanceReal.R(freq)
PlotTransmittance(freq, [np.abs(transReal),np.abs(reflecReal)],labels=['T','R'])

Animator(mesh, solNum,layer=layer,dispAndFree=True, fps=100)
input("Press [enter] to finish.")
#%%
'''
print('--- Comparison with analytical solution')
comparison=AnalyticComp(mesh, solNum, data["initialCond"])
solReal=comparison.AnalyticalSol(comparison.gridE,comparison.probeTime)
err=comparison.L2Error(solReal)
comparison.AnimatorTogether(solNum,solReal)
comparison.PrintErr(solNum,err)

print('=== Program finished')
'''