# Takes as input a list of filtered peaks, outputs mean, median, total count
# A bed file works fine for this

import sys

# Read in the size of all the enti
sizes = []
with open(sys.argv[1]) as f1:
	for line in f1:
		ld = line.strip().split()
		sizes.append(int(ld[2]) - int(ld[1]))
		
# Now calculate and output basic statistics
sizes.sort()
print "Total number of regions: " + str(len(sizes))
print "Total bases included: " + str(sum(sizes))
print "Mean region size: " + str(sum(sizes)/len(sizes))
print "Median region size: " + str(sizes[len(sizes)/2])
print "First quartile region size: " + str(sizes[len(sizes)/4])
print "Fourth quartile region size: " + str(sizes[3*len(sizes)/4])
print "First decile region size: " + str(sizes[len(sizes)/10])
print "Tenth decile region size: " + str(sizes[9*len(sizes)/10])