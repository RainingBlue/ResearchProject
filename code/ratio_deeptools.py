import numpy as np
import pysam
import argparse


###################################################################################
#             Add all reads to corresponding bin, then save ratio's               #
###################################################################################    

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
window_size = 5000000
chromos = []

bins = []
# Short count whole, total GC, Long count whole, short count seperate, long count seperate, start window, end window, chromosome, short, long
cur_bin = [0,0,0, cur_start, cur_end, cur_chr]

valid_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
              'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

def calcGC(seq):
    return ((seq.count('G') + seq.count('C'))/ len(seq))

def addRead(read, end, start):
    scale_window = window_size/(int(end)-int(start))
    if 100 <= abs(read.template_length) <= 150:
        cur_bin[0] += 1*scale_window
    #long frag count
    elif 150 < abs(read.template_length) <= 220:
        cur_bin[2] += 1*scale_window
    cur_bin[1] += calcGC(read.query_sequence.upper())


for read in bam_it:
    #Filter bad reads or too long reads
    if read.mapping_quality >= 30 and read.template_length > 0 and read.reference_name in valid_chr:
        if read.reference_name == cur_chr:
            #If fragment overlaps window, add the read to the window, if next_chr is not the 
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                addRead(read, cur_end, cur_start)
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
                        if (cur_bin[0] > 0 or cur_bin[2] > 0):
                           bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                        else:
                            bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
                addRead(read, cur_end, cur_start)
        else:
            chromos.append(read.reference_name)
            while(read.reference_name != cur_chr):
                    if (cur_bin[0] > 0 or cur_bin[2] > 0):
                       bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                    else:
                        bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                addRead(read, cur_end, cur_start)
            else:
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0) and cur_chr == read.reference_name):
                    if (cur_bin[0] > 0 or cur_bin[2] > 0):
                       bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                    else:
                        bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
                addRead(read, cur_end, cur_start)
if (cur_bin[0] > 0 or cur_bin[2] > 0):
   bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
else:
    bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
cur_bin = [0,0,0, cur_start, cur_end, cur_chr]

while cur_chr in valid_chr:
    cur_chr = next_chr
    cur_start = next_start
    cur_end = next_end
    next_window = windows.readline()
    bins.append((cur_bin[0], cur_bin[1], cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
    if next_window == "":
        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
        break
    else:
        next_chr,next_start,next_end = next_window.split("\t")
        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
arr = np.array(bins)
print(arr.shape)

ratio = []
for row in arr:
    if float(row[2]) > 0:
        ratio.append(float(row[0])/float(row[2]))
    else:
        ratio.append(np.nan)
np.save(my_bam_file + "_stats", ratio)
