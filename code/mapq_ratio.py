import numpy as np
import pysam
import argparse

############################################################################
#               Functions similarly to collect_bin_stats.py                #
############################################################################  

def calcGC(seq):
    return ((seq.count('G') + seq.count('C'))/ len(seq))
def addRead(read, cur_start, cur_end):
    scale_window = window_size/(int(cur_end)-int(cur_start))
    #short frag count
    if 100 <= abs(read.template_length) <= 150:
        if read.mapping_quality >= 30:
            cur_bin[0] += 1*scale_window
            cur_bin[5] += 1*scale_window
            cur_bin[7] += 1*scale_window
            cur_bin[9] += 1*scale_window
        elif read.mapping_quality >= 20:
            cur_bin[0] += 1*scale_window
            cur_bin[5] += 1*scale_window
            cur_bin[7] += 1*scale_window
        elif read.mapping_quality >= 5:
            cur_bin[0] += 1*scale_window
            cur_bin[5] += 1*scale_window
        else:
            cur_bin[0] += 1*scale_window
    #long frag count
    elif 150 < abs(read.template_length) <= 220:
        if read.mapping_quality >= 30:
            cur_bin[1] += 1*scale_window
            cur_bin[6] += 1*scale_window
            cur_bin[8] += 1*scale_window
            cur_bin[10] += 1*scale_window
        elif read.mapping_quality >= 20:
            cur_bin[1] += 1*scale_window
            cur_bin[6] += 1*scale_window
            cur_bin[8] += 1*scale_window
        elif read.mapping_quality >= 5:
            cur_bin[1] += 1*scale_window
            cur_bin[6] += 1*scale_window
        else:
            cur_bin[1] += 1*scale_window

valid_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
              'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

# Only calc ratio if more then 10 short and long fragments otherwise nan
def calc_ratio(short, long):
    return short/long if long > 10 and short > 10 else np.nan


parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("file", help="Input bam file to get distribution of", type=str)
parser.add_argument("windows", help="Input bed file containing windows", type=str)

args = parser.parse_args()
my_bam_file = args.file
windows = open(args.windows, 'r')
cur_chr,cur_start,cur_end = windows.readline().split("\t")
next_chr,next_start,next_end = windows.readline().split("\t")
imported = pysam.AlignmentFile(my_bam_file, mode = 'rb')
bam_it = imported.fetch()
ratios = []
window_size = 5000000
mapq_30 = 0
mapq_20 = 0
mapq_5 = 0
all_reads = 0

# Short count, Long count, start window, end window, chromosome, short_5, long_5, short_20, long_20, short_30, long_30
cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
bin_counter = 0
for read in bam_it:
    #Filter unpaired, secondary and other chromosomes
    if read.is_paired and not read.is_secondary and read.reference_name in valid_chr:
        if read.mapping_quality >= 30:
            mapq_30 += 1
        if read.mapping_quality >= 20:
            mapq_20 += 1
        if read.mapping_quality >= 5:
            mapq_5 += 1
        all_reads += 1
        if read.reference_name == cur_chr: 
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                addRead(read, cur_start, cur_end)
            elif (next_chr != read.reference_name):
                print("error")
            else:
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0)):
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        ratios.append([calc_ratio(cur_bin[0], cur_bin[1]), calc_ratio(cur_bin[5], cur_bin[6]), calc_ratio(cur_bin[7],cur_bin[8]), calc_ratio(cur_bin[9], cur_bin[10])])
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
                addRead(read, cur_start, cur_end)
        else:
            while(read.reference_name != cur_chr):
                    ratios.append([calc_ratio(cur_bin[0], cur_bin[1]), calc_ratio(cur_bin[5], cur_bin[6]), calc_ratio(cur_bin[7],cur_bin[8]), calc_ratio(cur_bin[9], cur_bin[10])])
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                addRead(read, cur_start, cur_end)
            else:
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0)):
                    ratios.append([calc_ratio(cur_bin[0], cur_bin[1]), calc_ratio(cur_bin[5], cur_bin[6]), calc_ratio(cur_bin[7],cur_bin[8]), calc_ratio(cur_bin[9], cur_bin[10])])
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
                addRead(read, cur_start, cur_end)

ratios.append([calc_ratio(cur_bin[0], cur_bin[1]), calc_ratio(cur_bin[5], cur_bin[6]), calc_ratio(cur_bin[7],cur_bin[8]), calc_ratio(cur_bin[9], cur_bin[10])])
cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
while cur_chr in valid_chr:
    cur_chr = next_chr
    cur_start = next_start
    cur_end = next_end
    next_window = windows.readline()
    ratios.append([calc_ratio(cur_bin[0], cur_bin[1]), calc_ratio(cur_bin[5], cur_bin[6]), calc_ratio(cur_bin[7],cur_bin[8]), calc_ratio(cur_bin[9], cur_bin[10])])
    if next_window == "":
        cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]
        break
    else:
        next_chr,next_start,next_end = next_window.split("\t")
        cur_bin = [0,0, cur_start, cur_end, cur_chr, 0,0,0,0,0,0]

# Z-score normalize the ratio's before saving them
ratios = (ratios - np.nanmean(ratios,axis=0))/np.nanstd(ratios,axis=0)
print("{},{},{},{}".format(all_reads, mapq_5, mapq_20, mapq_30))
np.save("out/" + my_bam_file.split('/')[-1] + "_mapq", ratios)


