from scipy import stats
import os
import numpy as np
import random

##############################################################################
#      Works same as nearest_neigbour.py, other file names for profiles      #
##############################################################################
directory = ["luad", "crc", "brca"]

correlation = []
healthy = []
cancer = []
for filename in os.listdir("Healthy"):
    f = os.path.join("Healthy", filename)
    # checking if it is a file
    if f.endswith("_mapq.npy") and os.path.isfile(f):
        healthy.append(np.load(f))
for dir in directory:
    for filename in os.listdir(dir):
        f = os.path.join(dir, filename)
        # checking if it is a file
        if f.endswith("_mapq.npy") and os.path.isfile(f):
            cancer.append(np.load(f))

indices_cancer_test = []
indices_healthy_test = []
for i in range(int(len(cancer)*0.3)):
    index = random.randrange(0,len(cancer))
    while (index in indices_cancer_test):
        index = random.randrange(0,len(cancer))
    indices_cancer_test.append(index)

for i in range(int(len(healthy)*0.3)):
    index = random.randrange(0,len(healthy))
    while (index in indices_healthy_test):
        index = random.randrange(0,len(healthy))
    indices_healthy_test.append(index)

cancer_test = []
cancer_train = []
healthy_test = []
healthy_train = []

for i in range(len(cancer)):
    if i in indices_cancer_test:
        cancer_test.append(cancer[i])
    else:
        cancer_train.append(cancer[i])
    
for i in range(len(healthy)):
    if i in indices_healthy_test:
        healthy_test.append(healthy[i])
    else:
        healthy_train.append(healthy[i])

cancer_test = np.asarray(cancer_test)
cancer_train = np.asarray(cancer_train)
healthy_test = np.asarray(healthy_test)
healthy_train = np.asarray(healthy_train)

print(cancer_test.shape)
print(cancer_train.shape)
print(healthy_test.shape)
print(healthy_train.shape)

cancer_classified_cancer_un = 0
cancer_classified_health_un = 0
healthy_classified_healthy_un = 0
healthy_classified_cancer_un = 0

cancer_classified_cancer_5 = 0
cancer_classified_health_5 = 0
healthy_classified_healthy_5 = 0
healthy_classified_cancer_5 = 0

cancer_classified_cancer_20 = 0
cancer_classified_health_20 = 0
healthy_classified_healthy_20 = 0
healthy_classified_cancer_20 = 0

cancer_classified_cancer_30 = 0
cancer_classified_health_30 = 0
healthy_classified_healthy_30 = 0
healthy_classified_cancer_30 = 0

for i in cancer_test:
    max_cancer_corr_un = 0
    max_healthy_corr_un = 0
    max_cancer_corr_5 = 0
    max_healthy_corr_5 = 0
    max_cancer_corr_20 = 0
    max_healthy_corr_20 = 0
    max_cancer_corr_30 = 0
    max_healthy_corr_30 = 0
    for j in cancer_train:
        corr_un = stats.spearmanr(i[:,0], j[:,0], nan_policy="omit").statistic
        corr_5 = stats.spearmanr(i[:,1], j[:,1], nan_policy="omit").statistic
        corr_20 = stats.spearmanr(i[:,2], j[:,2], nan_policy="omit").statistic
        corr_30 = stats.spearmanr(i[:,3],j[:,3], nan_policy="omit").statistic
        max_cancer_corr_un = max(max_cancer_corr_un, corr_un)
        max_cancer_corr_5 = max(max_cancer_corr_5, corr_5)
        max_cancer_corr_20 = max(max_cancer_corr_20, corr_20)
        max_cancer_corr_30 = max(max_cancer_corr_30, corr_30)
    for j in healthy_train:
        corr_un = stats.spearmanr(i[:,0], j[:,0], nan_policy="omit").statistic
        corr_5 = stats.spearmanr(i[:,1], j[:,1], nan_policy="omit").statistic
        corr_20 = stats.spearmanr(i[:,2], j[:,2], nan_policy="omit").statistic
        corr_30 = stats.spearmanr(i[:,3],j[:,3], nan_policy="omit").statistic
        max_healthy_corr_un = max(max_healthy_corr_un, corr_un)
        max_healthy_corr_5 = max(max_healthy_corr_5, corr_5)
        max_healthy_corr_20 = max(max_healthy_corr_20, corr_20)
        max_healthy_corr_30 = max(max_healthy_corr_30, corr_30)
    if max_cancer_corr_un > max_healthy_corr_un:
        cancer_classified_cancer_un += 1
    else:
        cancer_classified_health_un += 1
    if max_cancer_corr_5 > max_healthy_corr_5:
        cancer_classified_cancer_5 += 1
    else:
        cancer_classified_health_5 += 1
    if max_cancer_corr_20 > max_healthy_corr_20:
        cancer_classified_cancer_20 += 1
    else:
        cancer_classified_health_20 += 1
    if max_cancer_corr_30 > max_healthy_corr_30:
        cancer_classified_cancer_30 += 1
    else:
        cancer_classified_health_30 += 1

for i in healthy_test:
    max_cancer_corr_un = 0
    max_healthy_corr_un = 0
    max_cancer_corr_5 = 0
    max_healthy_corr_5 = 0
    max_cancer_corr_20 = 0
    max_healthy_corr_20 = 0
    max_cancer_corr_30 = 0
    max_healthy_corr_30 = 0
    for j in cancer_train:
        corr_un = stats.spearmanr(i[:,0], j[:,0], nan_policy="omit").statistic
        corr_5 = stats.spearmanr(i[:,1], j[:,1], nan_policy="omit").statistic
        corr_20 = stats.spearmanr(i[:,2], j[:,2], nan_policy="omit").statistic
        corr_30 = stats.spearmanr(i[:,3], j[:,3], nan_policy="omit").statistic
        max_cancer_corr_un = max(max_cancer_corr_un, corr_un)
        max_cancer_corr_5 = max(max_cancer_corr_5, corr_5)
        max_cancer_corr_20 = max(max_cancer_corr_20, corr_20)
        max_cancer_corr_30 = max(max_cancer_corr_30, corr_30)
    for j in healthy_train:
        corr_un = stats.spearmanr(i[:,0], j[:,0], nan_policy="omit").statistic
        corr_5 = stats.spearmanr(i[:,1], j[:,1], nan_policy="omit").statistic
        corr_20 = stats.spearmanr(i[:,2], j[:,2], nan_policy="omit").statistic
        corr_30 = stats.spearmanr(i[:,3], j[:,3], nan_policy="omit").statistic
        max_healthy_corr_un = max(max_healthy_corr_un, corr_un)
        max_healthy_corr_5 = max(max_healthy_corr_5, corr_5)
        max_healthy_corr_20 = max(max_healthy_corr_20, corr_20)
        max_healthy_corr_30 = max(max_healthy_corr_30, corr_30)
    if max_cancer_corr_un > max_healthy_corr_un:
        healthy_classified_cancer_un += 1
    else:
        healthy_classified_healthy_un += 1
    if max_cancer_corr_5 > max_healthy_corr_5:
        healthy_classified_cancer_5 += 1
    else:
        healthy_classified_healthy_5 += 1
    if max_cancer_corr_20 > max_healthy_corr_20:
        healthy_classified_cancer_20 += 1
    else:
        healthy_classified_healthy_20 += 1
    if max_cancer_corr_30 > max_healthy_corr_30:
        healthy_classified_cancer_30 += 1
    else:
        healthy_classified_healthy_30 += 1

print("Uncorrelated: Healthy classified healthy {}; Healthy classified cancer {}; Cancer classified healthy {}; Cancer classified cancer {}"
      .format(healthy_classified_healthy_un, healthy_classified_cancer_un, cancer_classified_health_un, cancer_classified_cancer_un))

print("MAPQ 5: Healthy classified healthy {}; Healthy classified cancer {}; Cancer classified healthy {}; Cancer classified cancer {}"
      .format(healthy_classified_healthy_5, healthy_classified_cancer_5, cancer_classified_health_5, cancer_classified_cancer_5))

print("MAPQ 20: Healthy classified healthy {}; Healthy classified cancer {}; Cancer classified healthy {}; Cancer classified cancer {}"
      .format(healthy_classified_healthy_20, healthy_classified_cancer_20, cancer_classified_health_20, cancer_classified_cancer_20))

print("MAPQ 30: Healthy classified healthy {}; Healthy classified cancer {}; Cancer classified healthy {}; Cancer classified cancer {}"
      .format(healthy_classified_healthy_30, healthy_classified_cancer_30, cancer_classified_health_30, cancer_classified_cancer_30))


        