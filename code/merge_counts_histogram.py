import os
import numpy as np
import matplotlib.pyplot as plt

# Max MAPQ in sample was 60
bins = []
hist = np.zeros(61)

# Get all samples binned to MAPQ
for filename in os.listdir("out"):
    f = os.path.join("out", filename)
    # checking if it is a file
    if f.endswith("_filtered.npy") and os.path.isfile(f):
        bins.append(np.load(f))

# Add all reads from bin to corresponfing bin
bins = np.asarray(bins)
for bin in bins:
    for j in range(61):
        hist[j] += bin[j]

# Make seperate array for MAPQ 0-4
mapq_5 = np.copy(hist[0:5])

# Normalize to get fractions
mapq_5 = mapq_5/np.sum(mapq_5)
hist = hist/np.sum(hist)

# Plot histograms
plt.bar(range(61), hist, width=1, align='center')
plt.xticks(range(61))
plt.ylabel("Percentage of all reads")
plt.ylim([0,1])
plt.xlim([-1, 61])
plt.show()

plt.bar(range(5), mapq_5, width=1, align='center')
plt.xticks(range(5))
plt.ylabel("Percentage of reads with MAPQ < 5")
plt.xlim([-1, 5])
plt.ylim([0,1])
plt.show()