import sys

# Takes as input:
# [1] - a peak calls file
# [2] - the snp calculation file
# [3] - the minimum threshold required to output a call, in bed format
# Outputs the location of all SNPs in these peaks where the significance score is over the threshold

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
	
# Read in peak calls
peak_regs = []
with open(sys.argv[1]) as f1:
	for line in f1:
		ld = line.strip().split()
		key = ld[0] + ':' + ld[4] + ':' + ld[5] + ':' + ld[6] # This is the location of the peak SNP, with the allele for absolute unambiguity
		region = [ld[0], int(ld[1]), int(ld[2]), key]
		peak_regs.append(region)
		
# Now go through calculation output file
threshold = 10**(-float(sys.argv[3]))
with open(sys.argv[2]) as f2:
	for line in f2:
		ld = line.strip().split()
		if float(ld[-1]) >= threshold:
			continue # Most SNPs will not be under the threshold
		
		region = [ld[0], int(ld[1]), int(ld[1]) + 1]
		 
		f_score = float(ld[4]) / (float(ld[4]) + float(ld[5])) # Proportion of reference calls in group 0
		m_score = float(ld[6]) / (float(ld[6]) + float(ld[7])) # Proportion of reference calls in group 1
		 
		for peak in peak_regs:
			if intersects(region, peak[:3]):
				print '\t'.join([str(x) for x in region]) + '\t' + str(f_score) + '\t' + str(m_score) + '\t' + ld[0] + ':' + ld[1] + ':' + ld[2] + ':' + ld[3] + '\t' + ld[-1] + '\t' + peak[-1] 