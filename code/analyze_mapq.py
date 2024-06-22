import argparse
import numpy as np
from scipy import stats

# Get file path from arguments
parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("median", help="Input npy file ratio of median", type=str)
parser.add_argument("measure", help="Input npy file counts per windows", type=str)
args = parser.parse_args()

# Load files
measurements = np.load(args.measure)
median = np.load(args.median)

# Split up median and sample into unfiltered, mapq 5, 20 and 30
normal_ratio_meas = measurements[:,0]
normal_ratio_median = median[0, :]

meas_5 = measurements[:,1]
median_5 = median[1,:]

meas_20 = measurements[:,2]
median_20 = median[2,:]

meas_30 = measurements[:,3]
median_30 = median[3,:]

#Print statistics
print(stats.ks_2samp(normal_ratio_meas, meas_5, nan_policy="omit"))
print(stats.ks_2samp(normal_ratio_meas, meas_20, nan_policy="omit"))
print(stats.ks_2samp(normal_ratio_meas, meas_30, nan_policy="omit"))
print("corr:" + str(stats.spearmanr(normal_ratio_median, normal_ratio_meas, nan_policy="omit")))
print("corr:" + str(stats.spearmanr(median_5, meas_5, nan_policy="omit")))
print("corr:" + str(stats.spearmanr(median_20, meas_20, nan_policy="omit")))
print("corr:" + str(stats.spearmanr(median_30, meas_30, nan_policy="omit")))