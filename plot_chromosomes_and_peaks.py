import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MultipleLocator
from math import log
import matplotlib.patches as mpatches

# Inputs:
# [1] the genome analysis file with simulation p-values
# [2] the region of interest to plot (ie, chrIV:12000000-12875000)
# [3] a bed file (ie, CSS bed)
# [4] a second bed file of interest (ie, called peaks)
# [5] the genome-wide significance threshold (set to <0 to not display)

# Note: for visual clarity, this will widen every bed entry to be at least 1% the width of the plot
region_size = sys.argv[2].split(':')[1].split('-')
region_size = int(region_size[1]) - int(region_size[0])


# Get CSS/HMM/other data for comparison
peak_file = sys.argv[3]
css_bed = sys.argv[4]

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
region = sys.argv[2]
chrom = region.split(':')[0]
start = int(region.split(':')[1].split('-')[0])
stop = int(region.split(':')[1].split('-')[1])

# Construct name of output file
output_name = sys.argv[1].split('.')[0] + '.' + sys.argv[1].split('.')[1] + '.' + region + '.png'

# Now go through genome analysis file and find just those entries that fall within this region
xs, ys = [], []
with open(sys.argv[1]) as f1:
	for line in f1:
		line_data = line.strip('\n').split('\t')

		# Filter out data from outside the region
		if line_data[0] != chrom:
			continue
		if int(line_data[1]) > stop or int(line_data[1]) < start:
			continue

		# For data in the region, store the coordinate and the -log10 p-value
		xs.append(float(line_data[1])/1000000) # Store as Mb  
		ys.append(-log(float(line_data[-1]), 10))


# Create the plot
fig = plt.figure(figsize=(12,9))
fig.suptitle('Probability of observed allele distribution', fontsize=27, fontweight='bold')
ax = fig.add_subplot(111)
xlab = 'Coordinates on ' + chrom + ' (Mb)'
ax.set_xlabel(xlab, fontsize=20)
ax.set_ylabel('log10 p-value', fontsize=20)
ax.plot(xs, ys, 'k*', ms = 5, zorder=2)

ax.grid(True)
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.grid(which = 'minor', linestyle="dashed")
plt.plot([min(xs), max(xs)], [float(sys.argv[5]), float(sys.argv[5])], 'r-', linewidth=1.5, zorder=1)
#plt.rc('grid', linestyle="dashed", color='gray')

# Add CSS/HMM data
for entry in somhmm:
	if entry == chrom:
		for pair in somhmm[entry]:
			plt.plot(pair, [max(ys) * 1.05, max(ys) * 1.05], 'm-', linewidth = 6)
for entry in css_regions:
	if entry == chrom:
		for pair in css_regions[entry]:
			plt.plot(pair, [max(ys) * 1.025, max(ys) * 1.025], 'g-', linewidth = 6)

# Set bounds on plot
plt.ylim([0, max(ys) * 1.08])
plt.xlim([min(xs), max(xs)])

# Add legend

legend1 = mpatches.Patch(color='m', label='CSS regions')
legend2 = mpatches.Patch(color='g', label='Best Pacific peaks')
legend3 = mpatches.Patch(color='r', label='Genome-wide significance')
plt.legend(handles=[legend1, legend2, legend3], fontsize = 12)

# Save figure
plt.savefig(output_name)
plt.show()
