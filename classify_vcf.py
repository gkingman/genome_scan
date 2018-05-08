# This program takes as input:
# [1] a VCF file with variants called on all populations of interest and 
# [2] a config file detailing which populations belong to which groups - this does NOT need to be the same as used in generating file 3
# [3] the stereotypical region file output

import sys
from scipy.stats import fisher_exact
from scipy import stats
from math import log
from math import factorial
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import json


# Import stereotypical calls
vectors = {}
base_to_key = {} # This is the quick lookup for going from any base to the appropriate peak key
with open(sys.argv[3]) as f0:
	for line in f0:
		line_data = line.strip().split()
		peak_key = line_data[7]
		
		if peak_key not in vectors:
			vectors[peak_key] = [] # All samples will go here
			
		base_key = line_data[5]
		base_to_key[base_key] = peak_key
		g0r, g1r = float(line_data[3]), float(line_data[4])
		
		# Don't allow 100% certainty
		bump = 0.01
		if g0r < 0.001:
			g0r += bump
		if g0r > 0.999:
			g0r += -bump
		if g1r < 0.001:
			g1r += bump
		if g1r > 0.999:
			g1r += -bump		
		
		# Calculate scores now
		RR = ((g0r**2 / (g0r**2 + g1r**2)) - 0.5) * 200
		Rr = ((g0r*(1-g0r) / (g0r*(1-g0r) + g1r*(1-g1r))) - 0.5) * 200
		rr = (((1-g0r)**2 / ((1-g0r)**2 + (1-g1r)**2)) - 0.5) * 200
		
		# Add to the vector
		vectors[peak_key].append([base_key, [RR, Rr, rr], []]) # Base key, genotype values, sample data


# Import population groups
group0 = []
group1 = []
combined_group_names = [] # Index to name
group_indices = {} # Name to index
with open(sys.argv[2]) as f1:
    for line in f1:
        line_data = line.strip().split()
        if len(line_data) != 2:
            continue
        group_indices[line_data[0]] = len(combined_group_names)   
        combined_group_names.append(line_data[0])        
        if line_data[1] == '0':
            group0.append(line_data[0])
        else:
            group1.append(line_data[0])

# Prepare to store entire dataset in one structure
total_groups = len(group0) + len(group1)
all_data = [[] for x in xrange(total_groups)]

# Now go through the VCF file and analyze line by line
indices = []
pop_names = []
snp_count = 0
with open(sys.argv[1]) as f2:
    for line in f2:

        if line[:2] == "##":
            continue # Skip initial headers

        if line[0] == '#': # Analyze main header for populations
            for index, entry in enumerate(line.strip().split()):
                if entry in group0 or entry in group1:
                    indices.append(index)
                    pop_names.append(entry)
            continue      

        line_data = line.strip().split()        
        base = line_data[0] + ':' + line_data[1] + ':' + line_data[3] + ':' + line_data[4]
                
        if base not in base_to_key:
        	continue
        
        peak_key = base_to_key[base]                
        
        pop_index = 0 # This is how many populations we have counted at this locus        

		# REFERENCE is not allowed
		
		# Find index
        i = 0
        for j in xrange(len(vectors[peak_key])):						
        	if vectors[peak_key][i][0] == base:
        		break
        	i += 1

        # Go through all entries
        for entry in indices:
            call = line_data[entry]
            if call == "./.":
            	vectors[peak_key][i][2].append(0)# If there's no data, use score 0 (no information to bias in either direction) 
            	pop_index += 1
                continue # Skip over if there's no data

            # Get and store the genotype call
            GT = call.split(":")[0]            
            if GT.count("0") == 2: # If call is 2 reference alleles, it gets the first precomputed value
            	vectors[peak_key][i][2].append(vectors[peak_key][i][1][0])
            	pop_index += 1
            elif GT.count("0") == 1: # If call is heterozygous, it has the second precomputed value
            	vectors[peak_key][i][2].append(vectors[peak_key][i][1][1])
            	pop_index += 1
            else: # If call is homozygous non-reference, it has the third precomputed value
                vectors[peak_key][i][2].append(vectors[peak_key][i][1][2])
            	pop_index += 1    

# Now go through and compute the data, first on the entire block
for peak in vectors:
	peak_base_count = len(vectors[peak])
	pop_count = len(vectors[peak][0][2])
	scores = [0] * pop_count # Initialize score matrix
		
	# Now go through and sum the scores
	for base in xrange(peak_base_count):
		for pop in xrange(pop_count):
			scores[pop] += vectors[peak][base][2][pop]
	
	# FIltering for testing speed only
	#if peak != "chrII:6129647:A:G": #chrVII:21439017:C:T chrIV:12797378:A:C
	#	continue
	
	# For each population, find the subset of contiguous SNPs that maximizes score
	min_size = int(min(10 + peak_base_count / 20., peak_base_count))
	matching_pops = []
	for pop in xrange(pop_count):
		pop_scores = []
		best = [-999999999999999, 0, 0, 0] # Score, start, stop
		for base in xrange(peak_base_count):
			pop_scores.append(vectors[peak][base][2][pop])
				
		for start in xrange(peak_base_count - min_size + 1):
			for stop in xrange(start + min_size, peak_base_count + 1):				
				score = sum(pop_scores[start:stop])
				if score > best[0]:
					best = [score, start, stop, score/(stop - start)]
		#print pop_names[pop], best
		if best[-1] > 75: # Use a very conservative threshold to flag a population as matching the g0 population
			matching_pops.append(pop_names[pop])
			
	# Output format	
	# One line per peak, consisting of:
	# [0] peak key
	# [1] Number of SNPs in peak
	# [2:] all populations passing the filter (with no further information output in this summary)
		
	print peak + '\t' + str(peak_base_count) + '\t' + '\t'.join(matching_pops)
	
	'''
	if peak == "chrXIX:14090964:C:T":
		for i in xrange(len(vectors[peak])):
			print vectors[peak][i] 
			print "Total bases examined: ", len(vectors[peak])
			print "Total populations examined: ", len(vectors[peak][i][2])
	'''		




'''

g0_counts = {}
g1_counts = {}
region_calls = [[] for x in xrange(len(pop_names))] # Indexed as follows: [Population index][SNP-index]
snp_index = []
for peak_key in vectors:	
	#print peak_key, len(vectors[peak_key][0]) # Print the peak key and how many SNPs we have in this peak to work with
	if peak_key not in snp_index:
		snp_index.append(peak_key)
	# Now, go through each population vector and see which reference vector it is closest to
	unclear = []
	g0 = []
	g1 = []
	for pop_index, population_vector in enumerate(vectors[peak_key][2]):
		#print pop_index, pop_names[pop_index], population_vector[0]
		
		if pop_names[pop_index] not in g0_counts:
			g0_counts[pop_names[pop_index]] = 0
			g1_counts[pop_names[pop_index]] = 0
		
		if len(population_vector) - population_vector.count('NaN') < 1:
			unclear.append(pop_index)	
			region_calls[pop_index].append('NaN')
			continue
			
		# If we have at least one good data point, let's try to classify
		g0_vector = vectors[peak_key][0]
		g1_vector = vectors[peak_key][1]
		g0_score = vector_distance(population_vector, g0_vector)
		g1_score = vector_distance(population_vector, g1_vector)
	
		# Require being a healthy margin closer to one score than to the other
		
		if g0_score * 1.5 < g1_score:
			g0.append(pop_index)
			g0_counts[pop_names[pop_index]] += 1
			region_calls[pop_index].append(0)
		elif g1_score * 1.5 < g0_score:
			g1.append(pop_index)
			g1_counts[pop_names[pop_index]] += 1
			region_calls[pop_index].append(1)
		else:
			unclear.append(pop_index)
			region_calls[pop_index].append('NaN')
			
	#print "Unclear populations:\t" + '\t'.join([pop_names[x] for x in unclear])
	#print "G0 populations:\t" + '\t'.join([pop_names[x] for x in g0])
	#print "G1 populations:\t" + '\t'.join([pop_names[x] for x in g1])

# Output region calls in organized manner here
# Header line first
# Then data, one peak per line, with populations as the columns
trans_list = map(list, zip(*region_calls))
out_name = sys.argv[2]
if '/' in out_name:
	out_name = out_name.split('/')[-1]
if '.' in out_name:
	out_name = out_name.split('.')[0]
out_name = out_name + "_pop_region_calls.txt"
f = open(out_name, 'w')
f.write('\t'.join(["Chromosome", "Peak_lower_bound-ish", "Peak_upper_bound-ish", "Peak_center", "Peak_unambiguous_ID"]) + '\t' + '\t'.join(pop_names) + '\n')
for peak_index, peak in enumerate(region_calls[0]):
	f.write(snp_index[peak_index].split(':')[0] + '\t' + str(int(snp_index[peak_index].split(':')[1]) - 1250) + '\t' + str(int(snp_index[peak_index].split(':')[1]) + 1250) + '\t' + snp_index[peak_index].split(':')[1] + '\t' + snp_index[peak_index] + '\t' + '\t'.join([str(x) for x in trans_list[peak_index]]) + '\n')
f.close()

# Output neat summary of results per population here
g0_pop_scores = []
g1_pop_scores = []
out_name = sys.argv[2]
if '/' in out_name:
	out_name = out_name.split('/')[-1]
if '.' in out_name:
	out_name = out_name.split('.')[0]
out_name = out_name + "_pop_summary.txt"
f = open(out_name, 'w')
for pindex, pop in enumerate(region_calls):
	group_name = "Group_0"
	if pop_names[pindex] in group1:
		group_name = "Group_1"
	f.write('\t'.join([str(x) for x in [pop_names[pindex], group_name, pop.count(0), pop.count(1), pop.count(0) * 100. / (pop.count(0) + pop.count(1))]]) + '\n') # Population name proportion of called regions matching group 0 
	
	# Incorporate group data 
	if pop_names[pindex] in group0:
		g0_pop_scores.append(pop.count(0) * 100. / (pop.count(0) + pop.count(1)))
	elif pop_names[pindex] in group1:
		g1_pop_scores.append(pop.count(0) * 100. / (pop.count(0) + pop.count(1)))
f.close()
print stats.ttest_ind(g0_pop_scores, g1_pop_scores)

# Now prepare to do PCA analysis
# First we need to impute the missing data
my_impute = custom_impute(region_calls) # Find the most similar vectors and use those to infer the missing calls

# Now we have all the data stored in vectors. Let's do PCA
pca = PCA()
pca.fit(my_impute)
#print(pca.components_)
print(pca.explained_variance_)
a = pca.explained_variance_
print "My imputing"
print "Percent variance explained by PCA 1: ", a[0]/sum(a) * 100
print "Percent variance explained by PCA 2: ", a[1]/sum(a) * 100

print pop_names
print my_impute

# Transform dimensionality
projected = pca.fit_transform(my_impute)
plt.figure(figsize=(100,50))
plt.plot(projected[:, 0], projected[:, 1], 'b.')

for label, x, y in zip(pop_names, projected[:, 0], projected[:, 1]):
    plt.annotate(
        label,
        xy=(x, y), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

out_name = sys.argv[2]
if '/' in out_name:
	out_name = out_name.split('/')[-1]
if '.' in out_name:
	out_name = out_name.split('.')[0]
out_name = out_name + "_pca_plot.png"

plt.savefig(out_name)
plt.close()
'''