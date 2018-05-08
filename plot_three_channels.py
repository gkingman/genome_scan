import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MultipleLocator
from math import log
import matplotlib.patches as mpatches

# Inputs:
# [1] the first genome analysis file with simulation p-values - Pacific
# [2] the second genome analysis file with simulation p-values - Atlantic
# [3] the third genome analysis file with simulation p-values - this one should be the "combined" file
# [4] the region of interest to plot (ie, chrIV:12000000-12875000)
# [5] Pacific peak calls
# [6] Atlantic peak calls
# [7] the genome-wide significance threshold (set to <0 to not display)

# Note: for visual clarity, this will widen every bed entry to be at least 1% the width of the plot
region_size = sys.argv[4].split(':')[1].split('-')
region_size = int(region_size[1]) - int(region_size[0])


# Get CSS/HMM/other data for comparison
peak_file = sys.argv[5]
css_bed = sys.argv[6]

# Read in BED1 regions
css_regions = {}
with open(css_bed, 'r') as f2:
	for line in f2:
		ld = line.strip().split()
		chrom, start, stop = ld[0], ld[1], ld[2]
		center = (float(start) + float(stop)) / 2
		width = (float(stop) - float(start))
		# If the width is less than .25% of the window, make it .25% of the window
		if width < region_size / 400:
			width = region_size / 400
		if chrom not in css_regions:
			css_regions[chrom] = []
		css_regions[chrom].append([(center-width/2)/1000000, (center+width/2)/1000000]) # Store coordinates in Mb

# Read in BED2 regions
somhmm = {}
with open(peak_file) as f3:
	for line in f3:
		ld = line.strip().split()
		chrom, start, stop = ld[0], ld[1], ld[2]
		center = (float(start) + float(stop)) / 2
		width = (float(stop) - float(start))
		if width < region_size / 400:
			width = region_size / 400
		if chrom not in somhmm:
			somhmm[chrom] = []
		somhmm[chrom].append([(center - width/2)/1000000, (center+width/2)/1000000])


# Begin analysis of specific input

# Parse input to determine region of interest
region = sys.argv[4]
chrom = region.split(':')[0]
start = int(region.split(':')[1].split('-')[0])
stop = int(region.split(':')[1].split('-')[1])

# Construct name of output file
output_name = sys.argv[1].split("/")[-1].split('.')[0] + '.' + region + '.png'

# Now go through genome analysis file and find just those entries that fall within this region
x1s, y1s = [], []
with open(sys.argv[1]) as f1:
	for line in f1:
		line_data = line.strip('\n').split('\t')

		# Filter out data from outside the region
		if line_data[0] != chrom:
			continue
		if int(line_data[1]) > stop or int(line_data[1]) < start:
			continue

		# For data in the region, store the coordinate and the -log10 p-value
		x1s.append(float(line_data[1])/1000000) # Store as Mb  
		y1s.append(-log(float(line_data[-1]), 10))

# Do the same with the second file
x2s, y2s = [], []
with open(sys.argv[2]) as f1:
	for line in f1:
		line_data = line.strip('\n').split('\t')

		# Filter out data from outside the region
		if line_data[0] != chrom:
			continue
		if int(line_data[1]) > stop or int(line_data[1]) < start:
			continue

		# For data in the region, store the coordinate and the -log10 p-value
		x2s.append(float(line_data[1])/1000000) # Store as Mb  
		y2s.append(-log(float(line_data[-1]), 10))

# And the third file
x3s, y3s = [], []
with open(sys.argv[3]) as f1:
	for line in f1:
		line_data = line.strip('\n').split('\t')

		# Filter out data from outside the region
		if line_data[0] != chrom:
			continue
		if int(line_data[1]) > stop or int(line_data[1]) < start:
			continue

		# For data in the region, store the coordinate and the -log10 p-value
		x3s.append(float(line_data[1])/1000000) # Store as Mb  
		y3s.append(-log(float(line_data[-1]), 10) / 1.5) # Divide by 1.5 such that the same scale is being used as in the analysis

# Create the plot
fig = plt.figure(figsize=(12,9))
fig.suptitle('Probability of observed allele distribution', fontsize=27, fontweight='bold')
ax = fig.add_subplot(111)
xlab = 'Coordinates on ' + chrom + ' (Mb)'
ax.set_xlabel(xlab, fontsize=20)
ax.set_ylabel('log10 p-value', fontsize=20)
ax.plot(x1s, y1s, '*', color="xkcd:cobalt blue", ms = 5, zorder=2)
ax.plot(x2s, y2s, '*', color="xkcd:darkgreen", ms = 5, zorder=2)
ax.plot(x3s, y3s, '*', color="xkcd:rusty orange", ms = 5, zorder=2)

ax.grid(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.grid(which = 'minor', linestyle="dashed")
plt.plot([min(x1s), max(x1s)], [float(sys.argv[7]), float(sys.argv[7])], 'r-', linewidth=1.5, zorder=1)
#plt.rc('grid', linestyle="dashed", color='gray')

# Add CSS/HMM data
for entry in somhmm:
	if entry == chrom:
		for pair in somhmm[entry]:
			plt.plot(pair, [max(y1s+y2s+y3s) * 1.05, max(y1s+y2s+y3s) * 1.05], '-', color="xkcd:cobalt blue", linewidth = 6)
for entry in css_regions:
	if entry == chrom:
		for pair in css_regions[entry]:
			plt.plot(pair, [max(y1s+y2s+y3s) * 1.025, max(y1s+y2s+y3s) * 1.025], '-', color="xkcd:darkgreen", linewidth = 6)

# Set bounds on plot
plt.ylim([0, max(y1s+y2s+y3s) * 1.08])
plt.xlim([min(x1s), max(x1s)])

# Add legend

legend1 = mpatches.Patch(color="xkcd:cobalt blue", label='Pacific calls')
legend2 = mpatches.Patch(color="xkcd:darkgreen", label='Atlantic calls')
legend3 = mpatches.Patch(color="xkcd:rusty orange", label='Combined calls')
legend4 = mpatches.Patch(color='r', label='Genome-wide significance')
plt.legend(handles=[legend1, legend2, legend3, legend4], fontsize = 12)

# Save figure
plt.savefig(output_name)
plt.show()
