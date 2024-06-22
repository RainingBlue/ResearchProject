import numpy as np
import pysam
import argparse

#############################################################
# Method to determine the GC-content of a sequence of bases #
#############################################################
def calcGC(seq):
    return ((seq.count('G') + seq.count('C'))/ len(seq))

#############################################################
# Method to add read to the current bin                     #
#############################################################
def addRead(read, cur_start, cur_end):
    #If the size of the window is smaller that the full window size, counts are scaled to compensate
    scale_window = window_size/(int(cur_end)-int(cur_start))
    #short frag count
    if 100 <= abs(read.template_length) <= 150:
        cur_bin[0] += 1*scale_window
    #long frag count
    elif 150 < abs(read.template_length) <= 220:
        cur_bin[2] += 1*scale_window
    cur_bin[1] += calcGC(read.query_sequence.upper())

# List of valid chromosomes, filter out any reads belonging to other regions (i.e mitochondrial / sex chromosomes)
valid_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
              'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

# Get path to output to, the path to the input BAM file and the sorted windows file to be used through arguments
parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("output", help="output npy file", type=str)
parser.add_argument("file", help="Input bam file to get distribution of", type=str)
parser.add_argument("windows", help="Input bed file containing windows", type=str)

args = parser.parse_args()
my_bam_file = args.file
windows = open(args.windows, 'r')
cur_chr,cur_start,cur_end = windows.readline().split("\t")
next_chr,next_start,next_end = windows.readline().split("\t")
imported = pysam.AlignmentFile(my_bam_file, mode = 'rb')
bam_it = imported.fetch()
bins = []
window_size = 100000

# Short count, total GC, Long count, start window, end window, chromosome
# Initialize current bin
cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
for read in bam_it:
    #Filter low mapq reads, unpaired, secondary and reads from other chromosomes
    if read.is_paired and not read.is_secondary and read.mapping_quality >= 30 and read.reference_name in valid_chr:
        # If read is not on the current chromosome, need to get bin of the read first before adding it
        if read.reference_name == cur_chr:
            #If fragment overlaps window, add the read to the window
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                addRead(read, cur_start, cur_end)
            #If read is not in the current bin and the next chromosome is different, read is not on chromosome
            elif (next_chr != read.reference_name):
                print("error")
            else:
                #Aslong as there is no overlap with the bin, get next bin and append the last bin
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0)):
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        # Bins with less than 10 reads for either short or long fragments are filtered out
                        if (cur_bin[2] > 10 and cur_bin[0] > 10):
                            bins.append((cur_bin[0], cur_bin[1]/(cur_bin[2]+cur_bin[0]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                        else:
                            bins.append((0, 0, 0, cur_bin[3], cur_bin[4], cur_bin[5]))
                        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
                addRead(read, cur_start, cur_end)
        else:
            # Get next bin until chromosome is correct, appending bins as we go
            while(read.reference_name != cur_chr):
                    if (cur_bin[2] > 10 and cur_bin[0] > 10):
                        bins.append((cur_bin[0], cur_bin[1]/(cur_bin[2]+cur_bin[0]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                    else:
                        bins.append((0, 0, 0, cur_bin[3], cur_bin[4], cur_bin[5]))
                    cur_chr = next_chr
                    cur_start = next_start
                    cur_end = next_end
                    next_window = windows.readline()
                    if next_window == "":
                        break
                    else:
                        next_chr,next_start,next_end = next_window.split("\t")
                        cur_bin = [0,0,0, cur_start, cur_end, cur_chr]
            # If there is overlap add the read to the bin, otherwise, go through bins again until overlap
            if read.get_overlap(int(cur_start), int(cur_end)+1):
                addRead(read, cur_start, cur_end)
            else:
                while((read.get_overlap(int(cur_start), int(cur_end)+1)==0)):
                    if (cur_bin[2] > 10 and cur_bin[0] > 10):
                        bins.append((cur_bin[0], cur_bin[1]/(cur_bin[2]+cur_bin[0]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
                        #print(bins[-1])
                    else:
                        bins.append((0, 0, 0, cur_bin[3], cur_bin[4], cur_bin[5]))
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
                addRead(read, cur_start, cur_end)

# Add the last bin to the list
if (cur_bin[2] > 10 and cur_bin[0] > 10):
    bins.append((cur_bin[0], cur_bin[1]/(cur_bin[2]+cur_bin[0]), cur_bin[2], cur_bin[3], cur_bin[4], cur_bin[5]))
else:
    bins.append((0, 0, 0, cur_bin[3], cur_bin[4], cur_bin[5]))

# Save counts per bin to file
arr = np.array(bins)
np.save(args.output + "_counts", arr)

