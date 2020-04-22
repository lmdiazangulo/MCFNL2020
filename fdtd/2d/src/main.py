
import json
import argparse
import os.path
import sys

from fdtd.mesh import Mesh
from fdtd.solver import Solver
from fdtd.viewer import Animator
from fdtd.comparison import AnalyticComp

print("=== Python FDTD 2D")


parser = argparse.ArgumentParser(description='Python FDTD 2D')
parser.add_argument('-i', '--input', nargs=1, type=str)
args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

inputFilename = ''.join(args.input).strip()
print("--- Reading file: %s"%(inputFilename))
data = json.load(open(inputFilename))

print('--- Initializing mesh')
mesh = Mesh(data["coordinates"], data["elements"], data["grid"])

print('--- Initializing solver')
solver = Solver(mesh, data["options"], data["probes"], data["sources"],data["initialCond"])

print('--- Solving')
solver.solve(data["options"]["finalTime"])

print('--- Visualizing')
solNum=solver.getProbes()[0]
Animator(mesh, solNum)

print('=== Program finished')