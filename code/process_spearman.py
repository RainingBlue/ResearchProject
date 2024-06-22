import argparse
import re
import csv

parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("file", help="Input bam file name", type=str)

args = parser.parse_args()
loess = args.file

# Wrote results of file to excel file

with open("spearman.csv", "a") as csv_file:
    writer = csv.writer(csv_file)
    with open(loess) as file:
        arr = []
        for line in file:
            if line.startswith("processing"):
                arr.append(re.findall('^processing (.*)', line)[0])
            if line.startswith("SignificanceResult(statistic"):
                arr.append(re.findall('^SignificanceResult\(statistic=(-?\d.\d*).*', line)[0])
    if len(arr) > 0:
        writer.writerow(arr)

with open("correlation.csv", "a") as csv_file:
    writer = csv.writer(csv_file)
    with open(loess) as file:
        arr = []
        for line in file:
            if line.startswith("processing"):
                arr.append(re.findall('^processing (.*)', line)[0])
            if line.startswith("corr"):
                arr.append(abs(float(re.findall('^.*SignificanceResult\(statistic=(-?\d.\d*).*', line)[0])))
    if len(arr) > 0:
        writer.writerow(arr)

with open("KS-test.csv", "a") as csv_file:
    writer = csv.writer(csv_file)
    with open(loess) as file:
        arr = []
        for line in file:
            if line.startswith("processing"):
                arr.append(re.findall('^processing (.*)', line)[0])
            if line.startswith("KstestResult(statistic="):
                print(line)
                arr.append(re.findall('^KstestResult\(statistic=(-?\d+.\d*).*', line)[0])
    if len(arr) > 0:
        writer.writerow(arr)

with open("nn.csv", "a") as csv_file:
    writer = csv.writer(csv_file)
    with open(loess) as file:
        arr = []
        for line in file:
            if line.startswith("processing"):
                arr.append(re.findall('^processing (.*)', line)[0])
            if line.startswith("['"):
                sample = line[1:-2].split(",")
                writer.writerow([sample[0][1:-1],sample[1][2:-1],sample[2][2:-1],sample[3][2:-1], sample[4][2:-1], sample[5][2:-1],sample[6][2:-1], sample[7][2:-1]])
    if len(arr) > 0:
        writer.writerow(arr)