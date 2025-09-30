#!/usr/bin/env python3
from getSimInfo import *
from gridToVTK import *
from getSimInfo import __macroNames__
import sys
import os

def usage_and_exit():
    print("Usage: python3 processSteps.py <velocity_set> <id>")
    print("Example: python3 processSteps.py D3Q19 000")
    sys.exit(1)

if len(sys.argv) != 3:
    usage_and_exit()

VELOCITY_SET = sys.argv[1]
SIM_ID       = sys.argv[2]

path = os.path.join("..", "bin", VELOCITY_SET, SIM_ID, "")

if not os.path.exists(path):
    print(f"Error: Path {path} does not exist.")
    sys.exit(1)

macroSteps = getMacroSteps(path)
info = retrieveSimInfo(path)

for step in macroSteps:
    macr = getMacrosFromStep(step, path)
    print("Processing step", step)
    saveVTK3D(macr, path, info['ID'] + "macr" + str(step).zfill(6), points=True)

    # optional cleanup
    # for macr_name in __macroNames__:
    #     filenames = get_filenames_macr(macr_name, path)
    #     step_filename = [f for f in filenames if f"{macr_name}{step:06d}.bin" in f]
    #     for file_to_delete in step_filename:
    #         os.remove(file_to_delete)
    #         print(f"Deleted {file_to_delete}")
