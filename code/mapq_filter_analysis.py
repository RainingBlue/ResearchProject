import numpy as np
import pysam
import argparse
# import matplotlib.pyplot as plt

valid_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
              'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']



parser = argparse.ArgumentParser("fragment_distribution")
parser.add_argument("file", help="Input bam file to get distribution of", type=str)
parser.add_argument("window", help="Input bam file to get distribution of", type=str)

args = parser.parse_args()
my_bam_file = args.file
imported = pysam.AlignmentFile(my_bam_file, mode = 'rb')

bam_it = imported.fetch()
bins = []
bins_filtered = []
windows = open(args.window, 'r')
window_size = 100000

filename = my_bam_file.split("/")[-1]

filtered_lengths = []
unfiltered_lengths = []
short_unfiltered = 0
long_unfiltered = 0
short_filtered = 0
long_filtered = 0
cur_bin_index = 0
cur_bin = []
cur_bin_filtered = []
mapq = []
bin_names = []
windows = windows.readlines()
window_index = 0
for line in windows:
    l = line.split("\t")
    bin_names.append(l[0]+":"+l[1]+"-"+l[2])

windows = [i.split("\t") for i in windows]
bin_names = np.asarray(bin_names)

for read in bam_it:
    #Filter bad reads
    if read.is_mapped and read.is_paired and read.reference_name in valid_chr and not read.is_secondary:
        # Add mapq to array for Histogram of mapq
        mapq.append(read.mapping_quality)
        # If read MAPQ < 5 add to unfiltered
        if read.mapping_quality < 5:
            if read.reference_name == windows[window_index][0] and read.get_overlap(int(windows[window_index][1]), int(windows[window_index][2])+1):
                cur_bin.append(abs(read.template_length))
                unfiltered_lengths.append(abs(read.template_length))
            else:
                while (read.reference_name != windows[window_index][0] or not read.get_overlap(int(windows[window_index][1]), int(windows[window_index][2])+1)):
                    bins.append(len(cur_bin))
                    cur_bin = []
                    bins_filtered.append(len(cur_bin_filtered))
                    cur_bin_filtered = []
                    window_index += 1
                cur_bin.append(abs(read.template_length))
            if 100 <= abs(read.template_length) <= 150:
                short_unfiltered += 1
            elif 150 < abs(read.template_length) <= 220:
                long_unfiltered += 1
        # Else add to filtered
        else:
            if read.reference_name == windows[window_index][0] and read.get_overlap(int(windows[window_index][1]), int(windows[window_index][2])+1):
                cur_bin_filtered.append(abs(read.template_length))
                filtered_lengths.append(abs(read.template_length))
            else:
                while (read.reference_name != windows[window_index][0] or not read.get_overlap(int(windows[window_index][1]), int(windows[window_index][2])+1)):
                    bins.append(len(cur_bin))
                    cur_bin = []
                    bins_filtered.append(len(cur_bin_filtered))
                    cur_bin_filtered = []
                    window_index += 1
                cur_bin_filtered.append(abs(read.template_length))
            if 100 <= abs(read.template_length) <= 150:
                short_filtered += 1
            elif 150 < abs(read.template_length) <= 220:
                long_filtered += 1
while window_index < len(windows):
    bins.append(len(cur_bin))
    bins_filtered.append(len(cur_bin_filtered))
    cur_bin_filtered = []
    cur_bin = []
    window_index += 1


# Print counts short, long for both filtered and unfiltered
print(short_filtered)
print(long_filtered)
print(short_unfiltered)
print(long_unfiltered)
# Print total amount of unfiltered, print max of a bin
print(np.sum(bins))
print(np.max(bins))

# Print stats
print("Mean filtered: {}".format(np.mean(filtered_lengths)))
if unfiltered_lengths:
    print("Mean unfiltered: {}".format(np.mean(unfiltered_lengths)))
else:
    print("Mean unfiltered: nan")

print("SD filtered: {}".format(np.std(filtered_lengths)))
if unfiltered_lengths:
    print("SD unfiltered: {}".format(np.std(unfiltered_lengths)))
else:
    print("SD unfiltered: nan")

# Save all MAPQ values
counts = np.bincount(mapq, minlength=70)
counts = counts/len(mapq)
filename = args.file.split("/")[-1]
np.save("out/" +  filename + "_filtered", counts)

###################################################################################################
###     Code below is for plotting the MAPQ histogram as well as the MAPQ < 5 reads per bin     ###
###################################################################################################

# plt.bar(range(61), counts, width=1, align='center')
# plt.xticks(range(61))
# plt.xlim([-1, 61])
# plt.show()

# fig = plt.figure()
# x = np.arange(len(bin_names))
# colors = []

# chr_1_filtered = []
# chr_1_unfiltered = []
# chr_1_names = []
# for i in range(len(bin_names)):
#     if bin_names[i].startswith("chr2"):
#         break
#     chr_1_filtered.append(bins_filtered[i])
#     chr_1_unfiltered.append(bins[i])
#     chr_1_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_1_filtered = chr_1_filtered/np.sum(chr_1_filtered)
# chr_1_unfiltered = chr_1_unfiltered/np.sum(chr_1_unfiltered)

# chr_2_filtered = []
# chr_2_unfiltered = []
# chr_2_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr2:"):
#         continue
#     chr_2_filtered.append(bins_filtered[i])
#     chr_2_unfiltered.append(bins[i])
#     chr_2_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_2_filtered = chr_2_filtered/np.sum(chr_2_filtered)
# chr_2_unfiltered = chr_2_unfiltered/np.sum(chr_2_unfiltered)

# chr_3_filtered = []
# chr_3_unfiltered = []
# chr_3_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr3"):
#         continue
#     chr_3_filtered.append(bins_filtered[i])
#     chr_3_unfiltered.append(bins[i])
#     chr_3_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_3_filtered = chr_3_filtered/np.sum(chr_3_filtered)
# chr_3_unfiltered = chr_3_unfiltered/np.sum(chr_3_unfiltered)

# chr_4_filtered = []
# chr_4_unfiltered = []
# chr_4_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr4"):
#         continue
#     chr_4_filtered.append(bins_filtered[i])
#     chr_4_unfiltered.append(bins[i])
#     chr_4_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_4_filtered = chr_4_filtered/np.sum(chr_4_filtered)
# chr_4_unfiltered = chr_4_unfiltered/np.sum(chr_4_unfiltered)

# chr_5_filtered = []
# chr_5_unfiltered = []
# chr_5_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr5"):
#         continue
#     chr_5_filtered.append(bins_filtered[i])
#     chr_5_unfiltered.append(bins[i])
#     chr_5_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_5_filtered = chr_5_filtered/np.sum(chr_5_filtered)
# chr_5_unfiltered = chr_5_unfiltered/np.sum(chr_5_unfiltered)

# chr_6_filtered = []
# chr_6_unfiltered = []
# chr_6_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr6"):
#         continue
#     chr_6_filtered.append(bins_filtered[i])
#     chr_6_unfiltered.append(bins[i])
#     chr_6_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_6_filtered = chr_6_filtered/np.sum(chr_6_filtered)
# chr_6_unfiltered = chr_6_unfiltered/np.sum(chr_6_unfiltered)

# chr_7_filtered = []
# chr_7_unfiltered = []
# chr_7_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr7"):
#         continue
#     chr_7_filtered.append(bins_filtered[i])
#     chr_7_unfiltered.append(bins[i])
#     chr_7_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_7_filtered = chr_7_filtered/np.sum(chr_7_filtered)
# chr_7_unfiltered = chr_7_unfiltered/np.sum(chr_7_unfiltered)

# chr_8_filtered = []
# chr_8_unfiltered = []
# chr_8_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr8"):
#         continue
#     chr_8_filtered.append(bins_filtered[i])
#     chr_8_unfiltered.append(bins[i])
#     chr_8_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_8_filtered = chr_8_filtered/np.sum(chr_8_filtered)
# chr_8_unfiltered = chr_8_unfiltered/np.sum(chr_8_unfiltered)

# chr_9_filtered = []
# chr_9_unfiltered = []
# chr_9_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr9"):
#         continue
#     chr_9_filtered.append(bins_filtered[i])
#     chr_9_unfiltered.append(bins[i])
#     chr_9_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_9_filtered = chr_9_filtered/np.sum(chr_9_filtered)
# chr_9_unfiltered = chr_9_unfiltered/np.sum(chr_9_unfiltered)

# chr_10_filtered = []
# chr_10_unfiltered = []
# chr_10_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr10"):
#         continue
#     chr_10_filtered.append(bins_filtered[i])
#     chr_10_unfiltered.append(bins[i])
#     chr_10_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_10_filtered = chr_10_filtered/np.sum(chr_10_filtered)
# chr_10_unfiltered = chr_10_unfiltered/np.sum(chr_10_unfiltered)

# chr_11_filtered = []
# chr_11_unfiltered = []
# chr_11_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr11"):
#         continue
#     chr_11_filtered.append(bins_filtered[i])
#     chr_11_unfiltered.append(bins[i])
#     chr_11_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_11_filtered = chr_11_filtered/np.sum(chr_11_filtered)
# chr_11_unfiltered = chr_11_unfiltered/np.sum(chr_11_unfiltered)

# chr_12_filtered = []
# chr_12_unfiltered = []
# chr_12_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr12"):
#         continue
#     chr_12_filtered.append(bins_filtered[i])
#     chr_12_unfiltered.append(bins[i])
#     chr_12_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_12_filtered = chr_12_filtered/np.sum(chr_12_filtered)
# chr_12_unfiltered = chr_12_unfiltered/np.sum(chr_12_unfiltered)


# chr_13_filtered = []
# chr_13_unfiltered = []
# chr_13_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr13"):
#         continue
#     chr_13_filtered.append(bins_filtered[i])
#     chr_13_unfiltered.append(bins[i])
#     chr_13_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_13_filtered = chr_13_filtered/np.sum(chr_13_filtered)
# chr_13_unfiltered = chr_13_unfiltered/np.sum(chr_13_unfiltered)

# chr_14_filtered = []
# chr_14_unfiltered = []
# chr_14_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr14"):
#         continue
#     chr_14_filtered.append(bins_filtered[i])
#     chr_14_unfiltered.append(bins[i])
#     chr_14_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_14_filtered = chr_14_filtered/np.sum(chr_14_filtered)
# chr_14_unfiltered = chr_14_unfiltered/np.sum(chr_14_unfiltered)

# chr_15_filtered = []
# chr_15_unfiltered = []
# chr_15_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr15"):
#         continue
#     chr_15_filtered.append(bins_filtered[i])
#     chr_15_unfiltered.append(bins[i])
#     chr_15_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_15_filtered = chr_15_filtered/np.sum(chr_15_filtered)
# chr_15_unfiltered = chr_15_unfiltered/np.sum(chr_15_unfiltered)

# chr_16_filtered = []
# chr_16_unfiltered = []
# chr_16_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr16"):
#         continue
#     chr_16_filtered.append(bins_filtered[i])
#     chr_16_unfiltered.append(bins[i])
#     chr_16_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_16_filtered = chr_16_filtered/np.sum(chr_16_filtered)
# chr_16_unfiltered = chr_16_unfiltered/np.sum(chr_16_unfiltered)

# chr_17_filtered = []
# chr_17_unfiltered = []
# chr_17_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr17"):
#         continue
#     chr_17_filtered.append(bins_filtered[i])
#     chr_17_unfiltered.append(bins[i])
#     chr_17_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_17_filtered = chr_17_filtered/np.sum(chr_17_filtered)
# chr_17_unfiltered = chr_17_unfiltered/np.sum(chr_17_unfiltered)

# chr_18_filtered = []
# chr_18_unfiltered = []
# chr_18_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr18"):
#         continue
#     chr_18_filtered.append(bins_filtered[i])
#     chr_18_unfiltered.append(bins[i])
#     chr_18_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_18_filtered = chr_18_filtered/np.sum(chr_18_filtered)
# chr_18_unfiltered = chr_18_unfiltered/np.sum(chr_18_unfiltered)

# chr_19_filtered = []
# chr_19_unfiltered = []
# chr_19_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr19"):
#         continue
#     chr_19_filtered.append(bins_filtered[i])
#     chr_19_unfiltered.append(bins[i])
#     chr_19_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_19_filtered = chr_19_filtered/np.sum(chr_19_filtered)
# chr_19_unfiltered = chr_19_unfiltered/np.sum(chr_19_unfiltered)

# chr_20_filtered = []
# chr_20_unfiltered = []
# chr_20_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr20"):
#         continue
#     chr_20_filtered.append(bins_filtered[i])
#     chr_20_unfiltered.append(bins[i])
#     chr_20_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_20_filtered = chr_20_filtered/np.sum(chr_20_filtered)
# chr_20_unfiltered = chr_20_unfiltered/np.sum(chr_20_unfiltered)

# chr_21_filtered = []
# chr_21_unfiltered = []
# chr_21_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr21"):
#         continue
#     chr_21_filtered.append(bins_filtered[i])
#     chr_21_unfiltered.append(bins[i])
#     chr_21_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_21_filtered = chr_21_filtered/np.sum(chr_21_filtered)
# chr_21_unfiltered = chr_21_unfiltered/np.sum(chr_21_unfiltered)

# chr_22_filtered = []
# chr_22_unfiltered = []
# chr_22_names = []
# for i in range(len(bin_names)):
#     if not bin_names[i].startswith("chr22"):
#         continue
#     chr_22_filtered.append(bins_filtered[i])
#     chr_22_unfiltered.append(bins[i])
#     chr_22_names.append(bin_names[i].split(":")[1].split("-")[0])

# chr_22_filtered = chr_22_filtered/np.sum(chr_22_filtered)
# chr_22_unfiltered = chr_22_unfiltered/np.sum(chr_22_unfiltered)
    

# bins = bins/np.sum(bins)
# bins_filtered = bins_filtered/np.sum(bins_filtered)

# for i in range(len(bins)):
#     if bins[i] > bins_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# labels = ["chr1:0-100000","chr2:0-100000","chr3:0-100000","chr4:0-100000","chr5:0-100000","chr6:0-100000","chr7:0-100000","chr8:0-100000","chr9:0-100000","chr11:0-100000","chr10:0-100000","chr12:0-100000","chr13:0-100000","chr14:0-100000","chr15:0-100000","chr16:0-100000","chr17:0-100000","chr18:0-100000","chr20:0-100000","chr19:0-100000","chr22:0-100000","chr21:0-100000"]
# plt.xticks([x[np.where(bin_names == "chr1:0-100000\n")[0][0]], x[np.where(bin_names == "chr2:0-100000\n")[0][0]], x[np.where(bin_names == "chr3:0-100000\n")[0][0]], x[np.where(bin_names == "chr4:0-100000\n")[0][0]], x[np.where(bin_names == "chr5:0-100000\n")[0][0]], x[np.where(bin_names == "chr6:0-100000\n")[0][0]], x[np.where(bin_names == "chr7:0-100000\n")[0][0]], x[np.where(bin_names == "chr8:0-100000\n")[0][0]], x[np.where(bin_names == "chr9:0-100000\n")[0][0]], x[np.where(bin_names == "chr11:0-100000\n")[0][0]], x[np.where(bin_names == "chr10:0-100000\n")[0][0]], x[np.where(bin_names == "chr12:0-100000\n")[0][0]], x[np.where(bin_names == "chr13:0-100000\n")[0][0]], x[np.where(bin_names == "chr14:0-100000\n")[0][0]], x[np.where(bin_names == "chr15:0-100000\n")[0][0]], x[np.where(bin_names == "chr16:0-100000\n")[0][0]], x[np.where(bin_names == "chr17:0-100000\n")[0][0]], x[np.where(bin_names == "chr18:0-100000\n")[0][0]], x[np.where(bin_names == "chr20:0-100000\n")[0][0]], x[np.where(bin_names == "chr19:0-100000\n")[0][0]], x[np.where(bin_names == "chr22:0-100000\n")[0][0]], x[np.where(bin_names == "chr21:0-100000\n")[0][0]]], labels, rotation=90)
# plt.scatter(x, bins, color=colors, s=0.9)
# fig.subplots_adjust(bottom=0.6)

# plt.show()

# fig = plt.figure(figsize=(12,6))

# colors = []
# for i in range(len(chr_1_filtered)):
#     if chr_1_unfiltered[i] > chr_1_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_1_names))
# plt.xticks(x[::50], chr_1_names[::50], rotation="vertical")
# plt.scatter(x, chr_1_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 1234, color = 'b', label = 'centromere')
# plt.title("CHR 1")

# plt.show()

# colors = []
# for i in range(len(chr_2_filtered)):
#     if chr_2_unfiltered[i] > chr_2_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_2_names))
# plt.xticks(x[::50], chr_2_names[::50], rotation="vertical")
# plt.scatter(x, chr_2_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 939, color = 'b', label = 'centromere')
# plt.title("CHR 2")

# plt.show()

# colors = []
# for i in range(len(chr_3_filtered)):
#     if chr_3_unfiltered[i] > chr_3_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_3_names))
# plt.xticks(x[::50], chr_3_names[::50], rotation="vertical")
# plt.scatter(x, chr_3_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 909, color = 'b', label = 'centromere')
# plt.title("CHR 3")

# plt.show()

# colors = []
# for i in range(len(chr_4_filtered)):
#     if chr_4_unfiltered[i] > chr_4_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_4_names))
# plt.xticks(x[::50], chr_4_names[::50], rotation="vertical")
# plt.scatter(x, chr_4_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 500, color = 'b', label = 'centromere')
# plt.title("CHR 4")

# plt.show()

# colors = []
# for i in range(len(chr_5_filtered)):
#     if chr_5_unfiltered[i] > chr_5_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_5_names))
# plt.xticks(x[::50], chr_5_names[::50], rotation="vertical")
# plt.scatter(x, chr_5_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 488, color = 'b', label = 'centromere')
# plt.title("CHR 5")

# plt.show()

# colors = []
# for i in range(len(chr_6_filtered)):
#     if chr_6_unfiltered[i] > chr_6_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_6_names))
# plt.xticks(x[::50], chr_6_names[::50], rotation="vertical")
# plt.scatter(x, chr_6_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 598, color = 'b', label = 'centromere')
# plt.title("CHR 6")

# plt.show()

# colors = []
# for i in range(len(chr_7_filtered)):
#     if chr_7_unfiltered[i] > chr_7_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_7_names))
# plt.xticks(x[::50], chr_7_names[::50], rotation="vertical")
# plt.scatter(x, chr_7_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 601, color = 'b', label = 'centromere')
# plt.title("CHR 7")

# plt.show()

# colors = []
# for i in range(len(chr_8_filtered)):
#     if chr_8_unfiltered[i] > chr_8_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_8_names))
# plt.xticks(x[::50], chr_8_names[::50], rotation="vertical")
# plt.scatter(x, chr_8_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 452, color = 'b', label = 'centromere')
# plt.title("CHR 8")

# plt.show()

# colors = []
# for i in range(len(chr_9_filtered)):
#     if chr_9_unfiltered[i] > chr_9_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_9_names))
# plt.xticks(x[::50], chr_9_names[::50], rotation="vertical")
# plt.scatter(x, chr_9_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 430, color = 'b', label = 'centromere')
# plt.title("CHR 9")

# plt.show()

# colors = []
# for i in range(len(chr_10_filtered)):
#     if chr_10_unfiltered[i] > chr_10_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_10_names))
# plt.xticks(x[::50], chr_10_names[::50], rotation="vertical")
# plt.scatter(x, chr_10_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 398, color = 'b', label = 'centromere')
# plt.title("CHR 10")

# plt.show()

# colors = []
# for i in range(len(chr_11_filtered)):
#     if chr_11_unfiltered[i] > chr_11_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_11_names))
# plt.xticks(x[::50], chr_11_names[::50], rotation="vertical")
# plt.scatter(x, chr_11_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 534, color = 'b', label = 'centromere')
# plt.title("CHR 11")

# plt.show()

# colors = []
# for i in range(len(chr_12_filtered)):
#     if chr_12_unfiltered[i] > chr_12_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_12_names))
# plt.xticks(x[::50], chr_12_names[::50], rotation="vertical")
# plt.scatter(x, chr_12_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 355, color = 'b', label = 'centromere')
# plt.title("CHR 12")

# plt.show()

# colors = []
# for i in range(len(chr_13_filtered)):
#     if chr_13_unfiltered[i] > chr_13_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_13_names))
# plt.xticks(x[::50], chr_13_names[::50], rotation="vertical")
# plt.scatter(x, chr_13_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 177, color = 'b', label = 'centromere')
# plt.title("CHR 13")
# plt.show()

# colors = []
# for i in range(len(chr_14_filtered)):
#     if chr_14_unfiltered[i] > chr_14_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_14_names))
# plt.xticks(x[::50], chr_14_names[::50], rotation="vertical")
# plt.scatter(x, chr_14_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 172, color = 'b', label = 'centromere')
# plt.title("CHR 14")

# plt.show()

# colors = []
# for i in range(len(chr_15_filtered)):
#     if chr_15_unfiltered[i] > chr_15_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_15_names))
# plt.xticks(x[::50], chr_15_names[::50], rotation="vertical")
# plt.scatter(x, chr_15_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 190, color = 'b', label = 'centromere')
# plt.title("CHR 15")

# plt.show()

# colors = []
# for i in range(len(chr_16_filtered)):
#     if chr_16_unfiltered[i] > chr_16_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_16_names))
# plt.xticks(x[::50], chr_16_names[::50], rotation="vertical")
# plt.scatter(x, chr_16_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 368, color = 'b', label = 'centromere')
# plt.title("CHR 16")

# plt.show()

# colors = []
# for i in range(len(chr_17_filtered)):
#     if chr_17_unfiltered[i] > chr_17_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_17_names))
# plt.xticks(x[::50], chr_17_names[::50], rotation="vertical")
# plt.scatter(x, chr_17_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 251, color = 'b', label = 'centromere')
# plt.title("CHR 17")

# plt.show()

# colors = []
# for i in range(len(chr_18_filtered)):
#     if chr_18_unfiltered[i] > chr_18_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_18_names))
# plt.xticks(x[::50], chr_18_names[::50], rotation="vertical")
# plt.scatter(x, chr_18_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 185, color = 'b', label = 'centromere')
# plt.title("CHR 18")

# plt.show()

# colors = []
# for i in range(len(chr_19_filtered)):
#     if chr_19_unfiltered[i] > chr_19_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_19_names))
# plt.xticks(x[::50], chr_19_names[::50], rotation="vertical")
# plt.scatter(x, chr_19_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 262, color = 'b', label = 'centromere')
# plt.title("CHR 19")

# plt.show()

# colors = []
# for i in range(len(chr_20_filtered)):
#     if chr_20_unfiltered[i] > chr_20_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_20_names))
# plt.xticks(x[::50], chr_20_names[::50], rotation="vertical")
# plt.scatter(x, chr_20_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 281, color = 'b', label = 'centromere')
# plt.title("CHR 20")

# plt.show()

# colors = []
# for i in range(len(chr_21_filtered)):
#     if chr_21_unfiltered[i] > chr_21_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_21_names))
# plt.xticks(x[::50], chr_21_names[::50], rotation="vertical")
# plt.scatter(x, chr_21_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 120, color = 'b', label = 'centromere')
# plt.title("CHR 21")

# plt.show()

# colors = []
# for i in range(len(chr_22_filtered)):
#     if chr_22_unfiltered[i] > chr_22_filtered[i]:
#         colors.append("g")
#     else: colors.append("r")

# x = np.arange(len(chr_22_names))
# plt.xticks(x[::50], chr_22_names[::50], rotation="vertical")
# plt.scatter(x, chr_22_unfiltered, color=colors)
# plt.tight_layout()
# plt.axvline(x = 150, color = 'b', label = 'centromere')
# plt.title("CHR 22")

# plt.show()
