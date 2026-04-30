import numpy as np
import matplotlib.pyplot as plt

import os, sys
#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import helperFunctions

def extractInfo(filename):
    ch_nums = []
    rho_values = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            # skip lines that don't have at least 4 columns (Chan, Freq, Rmaj, Rho)
            if len(parts) == 4 and parts[0].isdigit():
                try:
                    ch_num = int(parts[0])
                    rho = float(parts[-1])
                    ch_nums.append(ch_num)
                    rho_values.append(rho)
                except ValueError:
                    pass  # skip if not a number
    return ch_nums, rho_values

# Example usage:
shot = '203917'
time = '.03000'
filename = f'/home/grantr/codes/GENRAY_CQL3D_scripts/exp_scripts/ECElocs/ECE_{shot}{time}.txt'  # Replace with your actual file path
ch_nums, rho_ts = extractInfo(filename)
rho_p = helperFunctions.convertRhotorToRhopol(np.abs(rho_ts))
for i in range(len(ch_nums)):
    print(f'ch_num: {ch_nums[i]}, rho_ts: {rho_ts[i]}, rho_p: {rho_p[i]}')
