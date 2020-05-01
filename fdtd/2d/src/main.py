
import json
import argparse
import os.path
import sys

from fdtd.xdmf import Xdmf
from fdtd.mesh import Mesh
from fdtd.solver import Solver

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
solver = Solver(mesh, data["options"], data["probes"], data["sources"])

print('--- Solving')
solver.solve(data["options"]["finalTime"])

print('--- Writing output files')
(folder, file) = os.path.split(inputFilename)
caseName = os.path.splitext(file)[0]
xdmf = Xdmf(basename = caseName, format = "XML")
for p in solver.getProbes():
    xdmf.add(p)

open(xdmf.basename + ".xmf", "w+").write(xdmf.tostring().decode('utf-8'))

print('=== Program finished')