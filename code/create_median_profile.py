import os
import numpy as np
import random

# Directory where the profiles of the healthy samples are stored
rootdir = './healthy'

# Add all healthy profiles to a list
samples = []
for file in os.scandir(rootdir):
    if str(file).endswith('proc.bam_ratio.npy\'>'):
        filename = str(file).removeprefix('<DirEntry \'').removesuffix('_proc.bam_ratio.npy\'>')
        deeptools_file = rootdir + "/" + filename + ".bam_deeptools.bam_stats.npy"
        sample = np.load(file)
        deeptools = np.load(deeptools_file)
        sample = np.c_[sample, deeptools]
        sample = (sample - np.nanmean(sample, axis=0))/np.nanstd(sample, axis=0)
        samples.append(sample)

sample_indices = []
ratio_normal = []
ratio_whole = []
ratio_separete = []
ratio_deeptools = []
samples = np.asarray(samples)

# Pick 30 random indices of samples in the list, no repeats
for i in range(30):
    r_sample = random.randrange(0,samples.shape[0])
    while r_sample in sample_indices:
        r_sample = random.randrange(0,samples.shape[0])
    sample_indices.append(r_sample)
print(sample_indices)

# For every bin, append the median of the samples with the chosen indices, this is done for every method
for i in range(samples.shape[1]):
    ratio_normal.append(np.nanmedian(samples[sample_indices,i,0]))
    ratio_whole.append(np.nanmedian(samples[sample_indices,i,1]))
    ratio_separete.append(np.nanmedian(samples[sample_indices,i,2]))
    ratio_deeptools.append(np.nanmedian(samples[sample_indices,i,3]))

# Save the profile
np.save("median_profile", [ratio_normal, ratio_whole, ratio_separete, ratio_deeptools])