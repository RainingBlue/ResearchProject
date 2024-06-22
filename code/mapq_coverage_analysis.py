from scipy import stats
import os
import numpy as np
import random


directory = ["luad", "crc", "brca"]

healthy = []
cancer = []
# Output of MAPQ ratio reads retained
for filename in os.listdir("healthy/mapq/all_reads"):
    f = os.path.join("healthy/mapq/all_reads", filename)
    # checking if it is a file
    if f.endswith(".out") and os.path.isfile(f):
        sample = np.genfromtxt(f, delimiter=',')
        if len(sample) > 0:
            healthy.append(sample)
for dir in directory:
    for filename in os.listdir(dir + "/mapq/all_reads"):
        f = os.path.join(dir + "/mapq/all_reads", filename)
        # checking if it is a file
        if f.endswith(".out") and os.path.isfile(f):
            sample = np.genfromtxt(f, delimiter=',')
            if len(sample) > 0:
                cancer.append(sample)

# Print mean amount of fragments for healthy and cancer
healthy_mean = np.mean(healthy, axis=0)
cancer_mean = np.mean(cancer, axis=0)
print(healthy_mean)
print(cancer_mean)
h_perc = []
c_perc = []
# Print median percentage of fragments retained for all MAPQ values
for i in healthy:
    MAPQ_5 = 100*i[1]/i[0]
    MAPQ_20 = 100*i[2]/i[0]
    MAPQ_30 = 100*i[3]/i[0]
    h_perc.append([MAPQ_5,MAPQ_20,MAPQ_30])
print(np.median(h_perc, axis=0))
for i in cancer:
    MAPQ_5 = 100*i[1]/i[0]
    MAPQ_20 = 100*i[2]/i[0]
    MAPQ_30 = 100*i[3]/i[0]
    c_perc.append([MAPQ_5,MAPQ_20,MAPQ_30])
print(np.median(c_perc, axis=0))
    

