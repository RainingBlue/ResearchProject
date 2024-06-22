import numpy as np
import pysam
import argparse

parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("file", help="Input bam file to get distribution of", type=str)
parser.add_argument("windows", help="Input bed file containing windows", type=str)
parser.add_argument("scale", help="Input npy file containing scaling", type=str)

args = parser.parse_args()
my_bam_file = args.file
windows = open(args.windows, 'r')
cur_chr,cur_start,cur_end = windows.readline().split("\t")
next_chr,next_start,next_end = windows.readline().split("\t")
imported = pysam.AlignmentFile(my_bam_file, mode = 'rb')
bam_it = imported.fetch()
scale_array = np.load(args.scale)
window_size = 5000000
chromos = []

bins = []
# Short count whole, total GC, Long count whole, short count seperate, long count seperate, start window, end window, chromosome, short, long
cur_bin = [0,0,0, 0 ,0 , cur_start, cur_end, cur_chr, 0, 0]
bin_counter = 0

valid_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
              'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

def calcGC(seq):
    return ((seq.count('G') + seq.count('C'))/ len(seq))

# Bin counter is counter of current bin in scale file
def addRead(read, bin_counter):
    scale_start = int(scale_array[bin_counter][3])
    scale_end = int(scale_array[bin_counter][4])
    scale_chr = scale_array[bin_counter][5]
    # Get next scale bin until correct
    while scale_chr != read.reference_name:
        bin_counter += 1
        scale_start = int(scale_array[bin_counter][3])
        scale_end = int(scale_array[bin_counter][4])
        scale_chr = scale_array[bin_counter][5]
    while read.get_overlap(scale_start, scale_end+1) == 0:
        bin_counter += 1
        scale_start = int(scale_array[bin_counter][3])
        scale_end = int(scale_array[bin_counter][4])
        scale_chr = scale_array[bin_counter][5]
    scale_window = window_size/(scale_end-scale_start)
    # Add reads, scale using found scales through loess correction
    if 100 <= abs(read.template_length) <= 150:
        if scale_array[bin_counter][0] == 'nan':
            cur_bin[0] += scale_window
        else:
            cur_bin[0] += 1*float(scale_array[bin_counter][0])*scale_window
        if scale_array[bin_counter][1] == 'nan':
            cur_bin[3] += scale_window
        else:
            cur_bin[3] += 1*float(scale_array[bin_counter][1])*scale_window
        cur_bin[8] += 1*scale_window
    elif 150 < abs(read.template_length) <= 220:
        if scale_array[bin_counter][0] == 'nan':
            cur_bin[2] += scale_window
        else:
            cur_bin[2] += 1*float(scale_array[bin_counter][0])*scale_window
        if scale_array[bin_counter][2] == 'nan':
            cur_bin[4] += scale_window
        else:
            cur_bin[4] += 1*float(scale_array[bin_counter][2])*scale_window
        cur_bin[9] += 1*scale_window
    cur_bin[1] += calcGC(read.query_sequence.upper())
    return bin_counter

# Logic works the same as collect_bin_stats.py
for read in bam_it:
    #Filter bad reads
    if read.mapping_quality >= 30 and read.is_paired and not read.is_secondary and read.reference_name in valid_chr:
        if read.reference_name == cur_chr:
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                bin_counter = addRead(read, bin_counter)
            elif (next_chr != read.reference_name):
                print("error")
            else:
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0) and next_chr == cur_chr):
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        if (cur_bin[8] > 0 or cur_bin[9] > 0):
                            bins.append((cur_bin[0], cur_bin[1]/(cur_bin[8]+cur_bin[9]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
                        else:
                            bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
                        cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0, 0]
                bin_counter = addRead(read, bin_counter)
        else:
            chromos.append(read.reference_name)
            while(read.reference_name != cur_chr):
                    if (cur_bin[8] > 0 or cur_bin[9] > 0):
                        bins.append((cur_bin[0], cur_bin[1]/(cur_bin[8]+cur_bin[9]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
                    else:
                        bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0 , 0]
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                bin_counter = addRead(read, bin_counter)
            else:
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0) and cur_chr == read.reference_name):
                    if (cur_bin[8] > 0 or cur_bin[9] > 0):
                        bins.append((cur_bin[0], cur_bin[1]/(cur_bin[8]+cur_bin[9]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
                    else:
                        bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0, 0]
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0, 0]
                bin_counter = addRead(read, bin_counter)
if (cur_bin[8] > 10 or cur_bin[9] > 10):
    bins.append((cur_bin[0], cur_bin[1]/(cur_bin[8]+cur_bin[9]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
else:
    bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0, 0]

while cur_chr in valid_chr:
    cur_chr = next_chr
    cur_start = next_start
    cur_end = next_end
    next_window = windows.readline()
    bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5], cur_bin[6], cur_bin[7], cur_bin[8], cur_bin[9]))
    if next_window == "":
        cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0, 0]
        break
    else:
        next_chr,next_start,next_end = next_window.split("\t")
        cur_bin = [0,0,0,0,0, cur_start, cur_end, cur_chr, 0, 0]

arr = np.array(bins)
np.save(args.scale.removesuffix("_counts_gc_corrected.npy") + "_stats", arr)