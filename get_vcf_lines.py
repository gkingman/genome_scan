import sys

# Takes as input: 
# [1] a peaks output file
# [2] a VCF file
# Gives as output all lines from the VCF file where the variant described is the center of a peak, along with all headers
# Normal bedtools intersect doesn't work because there are a handful of bases with multiple variant calls, and we only want a specific one

# Store bases that we care about
bases = set()
with open(sys.argv[1]) as f1:
	for line in f1:
		ld = line.strip().split()
		pos_id = ld[0] + ':' + ld[4] + ':' + ld[5] + ':' + ld[6]
		bases.add(pos_id)

with open(sys.argv[2]) as f2:
	for line in f2:
		if line[0] == '#':
			print line.strip()
			continue
		ld = line.strip().split()
		pos_id = ld[0] + ':' + ld[1] + ':' + ld[3] + ':' + ld[4]		
		if pos_id in bases:
			print line.strip()