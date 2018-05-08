import sys

# Takes as input:
# [1] - the bed file of regions
# [2] - the calculation output file
# [3] - bed file one
# [4] - bed file two
# [5] - the genome-wide significance threshold
# [6] - how far to extend the plots on each side - recommend 10k


size = int(sys.argv[6])

with open(sys.argv[1]) as f1:
	for line in f1:
		ld = line.strip().split()
		reg1 = str(max(0, int(ld[1]) - size))
		region = ld[0] + ':' + reg1 + '-' + str(int(ld[2]) + size)
		print "python ~/genome_scan/calculate_all/scripts/plot_chromosomes_and_peaks.py " + sys.argv[2] + ' ' + region + ' ' + sys.argv[3] + " " + sys.argv[4] + " " + sys.argv[5] 