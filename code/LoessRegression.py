import numpy as np
names = "length, gc, count, start, end, chr"
import statsmodels.api as sm
# import matplotlib.pyplot as plt
from scipy import interpolate, stats
import argparse

# Get file from collect_bin_stats.py
parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("file", help="Input counts file name", type=str)

args = parser.parse_args()
my_bam_file = args.file

print(my_bam_file)

#Load into array
arr = np.load(my_bam_file + "_counts.npy")

# Extract data to seperate arrays
gc = arr[:,1]
short = arr[:,0]
long = arr[:,2]
chr = arr[:,5]
start = arr[:,3]
end = arr[:, 4]
short = np.asarray(short).astype(float)
long = np.asarray(long).astype(float)
count = np.add(short, long)

short_f = np.asarray(short).astype(float)
long_f = np.asarray(long).astype(float)
start_f = np.asarray(start).astype(int)
end_f = np.asarray(end).astype(int)
gc_f = np.asarray(gc).astype(float)
count_f = np.add(short_f, long_f)

# Get 1,99th percentile
perecntile = (np.percentile(count_f, [1, 99]))
perecntile_short = (np.percentile(short_f, [1, 99]))
perecntile_long = (np.percentile(long_f, [1, 99]))

gc_filtered = []
count_filtered = []
short_filtered = []
gc_short = []
long_filtered = []
gc_long = []

# Keep only bins wich fall between 1,99th percentile for creating curve
for i in range(len(count_f)):
    if short_f[i] > perecntile_short[0] and short_f[i] < perecntile[1] and long_f[i] > perecntile_long[0] and long_f[i] < perecntile_long[1] and gc_f[i] > 0.35:
        count_filtered.append(count_f[i])
        short_filtered.append(short_f[i])
        long_filtered.append(long_f[i])
        gc_filtered.append(gc_f[i])

# Get medians
median = np.median(count_f)
median_short = np.median(short_f)
median_long = np.median(long_f)

# Create loess regression curves for all fragments, short fragments and long fragments
z = sm.nonparametric.lowess(count_filtered,gc_filtered, frac=0.75)
z_short = sm.nonparametric.lowess(short_filtered,gc_filtered, frac=0.75)
z_long = sm.nonparametric.lowess(long_filtered,gc_filtered, frac=0.75)

# Create interpolation models for all fragments, short fragments and long fragments
f = interpolate.interp1d(z[:,0],z[:,1], fill_value='extrapolate')
f_short = interpolate.interp1d(z_short[:,0],z_short[:,1], fill_value='extrapolate')
f_long = interpolate.interp1d(z_long[:,0],z_long[:,1], fill_value='extrapolate')

short_whole = []
long_whole = []
short_corrected = []
long_corrected = []
scale_w = []
scale_short_w = []
scale_long_w = []

# For all reads generate a scale LOESS whole value, scale LOESS short value and scale LOESS long value
# And apply scales to counts for GC-bias calculation purposes
for i in range(len(count_f)):
    # Predict what count should be, given the GC-content of a bin
    predict = f(gc_f[i])
    predict_short = f_short(gc_f[i])
    predict_long = f_long(gc_f[i])

    # Get the actual count of the bin
    coverage = count_f[i]
    coverage_short = short_f[i]
    coverage_long = long_f[i]
    # If the count is larger than 0, subtract predicted from the coverage, leaving only the coverage which is not caused by 
    # GC-content. Add the median back to transform back to original scale, divide it by the found coverage to give the scaling factor
    # this produces. If 0 counts, add nan to scale.
    if count_f[i] > 0:
        scale = ((coverage - predict + median)/count_f[i])
    else:
        scale = float('nan')
    if (short[i] > 0):
        scale_short = ((coverage_short - predict_short + median_short)/short_f[i])
    else:
        scale_short = float('nan')
    if (long[i] > 0):
        scale_long = ((coverage_long - predict_long + median_long)/long_f[i])
    else:
        scale_long = float('nan')
    scale_w.append(scale)
    scale_short_w.append(scale_short)
    scale_long_w.append(scale_long)
    short_whole.append(short_f[i]*scale)
    long_whole.append(long_f[i]*scale)
    short_corrected.append(short_f[i]*scale_short)
    long_corrected.append(long_f[i]*scale_long)


count_proc = [short_whole[i] + long_whole[i] for i in range(len(short_f))]
perecntile = (np.percentile(count_proc, [1, 99]))
gc_filtered = []
count_filtered = []
short_filtered = []

# GC-bias before processing
print(stats.spearmanr(gc_f, short_f, nan_policy='omit'))
print(stats.spearmanr(gc_f, long_f, nan_policy='omit'))
# GC-bias LOESS whole
print(stats.spearmanr(gc_f, short_whole, nan_policy='omit'))
print(stats.spearmanr(gc_f, long_whole, nan_policy='omit'))
# GC-bias LOESS seperate
print(stats.spearmanr(gc_f, short_corrected, nan_policy='omit'))
print(stats.spearmanr(gc_f, long_corrected, nan_policy='omit'))

scale_w = np.asarray(scale_w)
scale_short_w = np.asarray(scale_short_w)
scale_long_w = np.asarray(scale_long_w)

# For all bins, save the different scaling factors
results = []
for i in range(len(count)):
    cur_bin = (scale_w[i], scale_short_w[i], scale_long_w[i], start[i], end[i], chr[i])
    results.append(cur_bin)

results = np.array(results)
np.save(my_bam_file + "_counts_gc_corrected", results)