# Takes as input:
# [1] The peak file that defines the regions of interest - last value should be p-value for highest point in peak
# [2] The number of SNPs checked
# [3] The desired false discovery rate
# Outputs all peaks that pass the false discovery threshold

import sys
fdr = float(sys.argv[3])
num = int(sys.argv[2])

# Read in all peaks, store score and full line data
peaks = []
with open(sys.argv[1]) as f0:
    for line in f0:
        ld = line.strip().split()
        key = ld[0] + ':' + ld[1] + '-' + ld[2]
        score = float(ld[-1])
        peaks.append((score, line.strip())) # Store the score, a placeholder value for if we accept this peak (none accepted yet)

peaks.sort()

# Now go through the peaks and find the last one that counts
stopkey = ''
count = 1
for peak in peaks:
	score = peak[0]
	threshold = fdr / num * count
	if score <= threshold:
		stopkey = peak[1]
	count += 1

# Now output all the peaks that passed
for peak in peaks:
	print peak[1]
	if peak[1] == stopkey:
		break


'''
key = ''
num = int(sys.argv[2])
fdr = float(sys.argv[3])
count = 0
with open(sys.argv[1]) as f1:
    for line in f1:
        # get basic info
        ld = line.strip().split()
        chrom, pos = ld[0], int(ld[1])
        
        # Are we in an unclaimed peak?
        j = 0
        while j < len(peaks):
            i = j
            j += 1
            peak = peaks[i]
            peak_key = peak[1]
            pchr, pstart, pstop = peak_key.split(':')[0], int(peak_key.split(':')[1].split('-')[0]), int(peak_key.split('-')[1])
            if pchr != chrom:
                continue
            # We are in a peak
            if pstart <= pos and pos <= pstop:
                if len(peaks[i]) == 2:
                    peaks[i].append(line)
                    peaks[i][0] = float(ld[-2])

stopkey = ''
count = 1
peaks.sort()              
for peak in peaks:
	score = peak[0]
	key = peak[1]
	ld = peak[2].strip().split()
	observed = int(ld[10])
	threshold = fdr / num * count
        n = int(ld[9])
        sd = (n * threshold * (1 - threshold)) ** 0.5
        required = n *threshold - 2 * sd
        
        count += 1

	# Each time we observe something that matches, store that key                                                   
        if observed < required:
        	stopkey = key
	
	
	
for peak in peaks:
	peak_key = peak[1]
	ld = peak[2].strip().split()	
	output = peak_key.split(':')[0] + '\t' + str(peak_key.split(':')[1].split('-')[0]) + '\t' + str(peak_key.split('-')[1])
	output = output + '\t' + str(int(output.split()[2]) - int(output.split()[1])) + '\t'
	output = output + '\t'.join(ld[1:])
	print output
	if peak[1] == stopkey:
		break
'''