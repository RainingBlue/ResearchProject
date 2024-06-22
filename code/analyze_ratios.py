import argparse
import numpy as np
from scipy import stats

# Get median healthy profile, the counts after correction and the deeptools file from args
parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("file1", help="Input npy file ratio of median", type=str)
parser.add_argument("file2", help="Input npy file counts per windows", type=str)
parser.add_argument("deeptools", help="Input npy file with ratio for deeptools", type=str)
args = parser.parse_args()

# Load files
measurements = np.load(args.file2)
median = np.load(args.file1)
deeptools = np.load(args.deeptools)

#Split nedian and counts into different methods, z-score normalize the ratio's
normal_ratio_meas = (measurements[:,0]-np.nanmean(measurements[:,0]))/np.nanstd(measurements[:,0])
normal_ratio_median = median[0, :]

whole_ratio_meas = (measurements[:,1]-np.nanmean(measurements[:,1]))/np.nanstd(measurements[:,1])
whole_ratio_median = median[1,:]

seperate_ratio_meas = (measurements[:,2]-np.nanmean(measurements[:,2]))/np.nanstd(measurements[:,2])
seperate_ratio_median = median[2,:]

deeptools = (deeptools - np.nanmean(deeptools))/np.nanstd(deeptools)
deeptools_median = median[3,:]

#Print the statistics for further processing
print(stats.ks_2samp(normal_ratio_meas, whole_ratio_meas, nan_policy="omit"))
print(stats.ks_2samp(normal_ratio_meas, seperate_ratio_meas, nan_policy="omit"))
print(stats.ks_2samp(deeptools, normal_ratio_meas, nan_policy="omit"))
print("corr:" + str(stats.spearmanr(normal_ratio_median, normal_ratio_meas, nan_policy="omit")))
print("corr:" + str(stats.spearmanr(whole_ratio_median, whole_ratio_meas, nan_policy="omit")))
print("corr:" + str(stats.spearmanr(seperate_ratio_median, seperate_ratio_meas, nan_policy="omit")))
print("corr:" + str(stats.spearmanr(deeptools_median, deeptools, nan_policy="omit")))