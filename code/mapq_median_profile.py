import os
import numpy as np
import random
rootdir = './out'

# Same logic as creating normal median profile, different file name

samples = []
for file in os.scandir(rootdir):
    if str(file).endswith('_mapq.npy\'>'):
        sample = np.load(file)
        samples.append(sample)

sample_indices = []
ratio_normal = []
ratio_5 = []
ratio_20 = []
ratio_30 = []
samples = np.asarray(samples)
print(samples.shape)

for i in range(30):
    r_sample = random.randrange(0,samples.shape[0])
    while r_sample in sample_indices:
        r_sample = random.randrange(0,samples.shape[0])
    sample_indices.append(r_sample)
print(sample_indices)

for i in range(samples.shape[1]):
    ratio_normal.append(np.nanmedian(samples[sample_indices,i,0]))
    ratio_5.append(np.nanmedian(samples[sample_indices,i,1]))
    ratio_20.append(np.nanmedian(samples[sample_indices,i,2]))
    ratio_30.append(np.nanmedian(samples[sample_indices,i,3]))

np.save("mapq_profile", [ratio_normal, ratio_5, ratio_20, ratio_30])