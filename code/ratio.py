import argparse
import numpy as np

#########################################################################
#             create ratio's from gc-corrected counts                   #
#########################################################################

parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("output", help="output filename", type=str)
parser.add_argument("windows", help="Input npy file counts per windows", type=str)

args = parser.parse_args()

windows = np.load(args.windows)

ratios = []

for i in range(len(windows)):
    ratio = []
    #Ratio without GC-correction
    if float(windows[i, 9]) > 0:
        ratio.append(float(windows[i,8])/float(windows[i, 9]))
    else:
        ratio.append(float('nan'))
    #Ratio whole correction
    if float(windows[i, 2]) > 0:
        ratio.append(float(windows[i,0])/float(windows[i, 2]))
    else:
        ratio.append(float('nan'))
    #Ratio separete correction
    if float(windows[i, 4]) > 0:
        ratio.append(float(windows[i,3])/float(windows[i, 4]))
    else:
        ratio.append(float('nan'))
    ratios.append(ratio)

np.save(args.output + "_ratio", ratios)