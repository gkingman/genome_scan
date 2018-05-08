import sys

# Takes as input:
# [1] - the classification output file
# [2] - a bed file to filter against
# [3] - how close can a bed entry be before we filter it out?
# [4] - the config file to be analyzed
# Outputs only the classifications that do not intersect with the bed file, with the number of samples from each group listed

def intersects(region1, region2):
	'''
	Returns True if the two regions intersect, otherwise false.
	Each region should be in [chrom, start, stop] format. start and stop should be ints
	'''
	
	chrom1, start1, stop1 = region1
	chrom2, start2, stop2 = region2
	
	if chrom1 != chrom2:
		return False # No intersection possible if separate chromosomes
		
	# Now consider the four possible ways to have intersection
	if start2 <= start1 and start1 <= stop2:
		return True # 1 start in middle
	if start2 <= stop1 and stop1 <= stop2:
		return True # 1 stop in middle		
	if start1 <= start2 and start2 <= stop1:
		return True # 2 start in middle
	if start1 <= stop2 and stop2 <= stop1:
		return True # 2 stop in middle
		
	return False # If none of these are true, there cannot be an intersect

# Import population groups
group_indices = {} # Name to index
with open(sys.argv[4]) as f1:
    for line in f1:
        line_data = line.strip().split()
        if len(line_data) != 2:
            continue     
        group_indices[line_data[0]] = int(line_data[1])

# Read in bed file
bed_regions = []
dist = int(sys.argv[3])
with open(sys.argv[2]) as f2:
	for line in f2:
		ld = line.strip().split()
		bed_regions.append([ld[0], int(ld[1]) - dist, int(ld[2]) + dist])

# Now read in classifications
with open(sys.argv[1]) as f1:
	for line in f1:
		ld = line.strip().split()
		region = [ld[0].split(':')[0], int(ld[0].split(':')[1]), int(ld[0].split(':')[1]) + 1]
		overlaps = False
		for bed_region in bed_regions:
			if intersects(bed_region, region):
				overlaps = True
		if not overlaps:
			counts = [0, 0]
			for pop in ld[2:]:
				counts[group_indices[pop]] += 1
			if sum(counts) > 0:
				prop = 100.*counts[0] / sum(counts)
			else:
				prop = 0
			print ld[0] + '\t' + ld[1] + '\t' + '\t'.join([str(x) for x in counts]) + '\t' + str(sum(counts)) + '\t' + str(prop) + '\t' + '\t'.join(ld[2:])