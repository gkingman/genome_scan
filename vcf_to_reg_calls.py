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

def vector_distance(vector1, vector2):
	'''
	Takes as input two vectors, both of dimension N. Returns the distance between the two vectors, with smaller scores obviously meaning that they are more similar.
	If either vector is missing a data point, skip over that dimension.
	'''
	dist = 0.
	for index in xrange(len(vector1)):
		if vector1[index] == 'NaN' or vector2[index] == 'NaN':
			continue
		dist += (vector1[index] - vector2[index]) ** 2
	return dist ** 0.5

snp_number = 2
comp_number = 3 # This is how mnay other vectors we will average to impute missing values

def custom_impute(data_in):
	'''
		Finds the closest three samples at nearby SNPs and imputes their value here 
	'''
	
	# Start constructing a new data set with the original and imputed values
	data_out = [[] for x in xrange(total_groups)]
	for snp_index in xrange(len(data_in[0])):
		for pop_index in xrange(len(data_in)):
			# If we have a good base call, store it, no need to do anything else
			if data_in[pop_index][snp_index] != 'NaN':
				data_out[pop_index].append(data_in[pop_index][snp_index])
			else:
				lower_bound = 0 #max(0, snp_index -  snp_number)
				upper_bound = len(data_in[0]) #min(len(data_in[0]), snp_index +  snp_number)
				
				# This is the vector of datapoints we will use to find the best match
				sample_vector = data_in[pop_index][lower_bound:upper_bound]
				
				# Now to collect other vectors to compare against
				comparison_vectors = []
				for comp_index in xrange(len(data_in)):
					if comp_index == pop_index:						
						continue # We don't want to match with ourself
					if data_in[comp_index][snp_index] == 'NaN':												
						continue # We can't use another unknown sample to inform us					
					comparison_vectors.append(data_in[comp_index][lower_bound:upper_bound])
				
				# Now to add the map distance to the end of each comparison vector				
				for comp_vector in comparison_vectors:
										
					dif_sum = 0
					bases = 0 # How many bases have we successfully examined?
					for base_index, base in enumerate(sample_vector):
						#print base_index, len(comp_vector), len(sample_vector)
						if comp_vector[base_index] == 'NaN' or sample_vector[base_index] == 'NaN':
							continue # Skip over each base with a non-call in either position
						bases += 1
						difference = abs(comp_vector[base_index] - sample_vector[base_index])
						dif_sum += difference
						
					# If we have too few bases, give this comparison a very high distance score
					# This should be true for a very small number of cases
					if bases < 10:						
						comp_vector.append(1)						
					else:
						comp_vector.append(float(dif_sum)/bases) # Now store the average distance
						
				# Now find the closest vectors and use these to impute		
				comparison_vectors.sort(key=lambda x: x[-1])
																	
				# Make sure we know the proper SNP to retrieve
				entry = snp_index # snp_number
				#if lower_bound == 0:
				#	entry = len(sample_vector) - snp_number # If we have fewer entries than usual to the left, count from the right
				imputed_value = 0.				
				for sample in xrange(min(comp_number, len(comparison_vectors))):
					imputed_value += comparison_vectors[sample][entry]
				imputed_value = imputed_value / comp_number
				data_out[pop_index].append(imputed_value)
	return data_out

# Import vectors
vectors = {}
base_to_key = {} # This is the quick lookup for going from any base to the appropriate peak key
with open(sys.argv[3]) as f0:
	for line in f0:
		line_data = line.strip().split()
		key = line_data[7]
		if key not in vectors:
			vectors[key] = [[], [], []] # Group 0 ref %, Group 1 ref %, all samples will go here
		vectors[key][0].append(float(line_data[3]))
		vectors[key][1].append(float(line_data[4]))
		base = line_data[5]
		base_to_key[base] = key

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
                    for peak in vectors:
                    	vectors[peak][2].append([]) # Initialize vectors
            continue      

        line_data = line.strip().split()
        
        base = line_data[0] + ':' + line_data[1] + ':' + line_data[3] + ':' + line_data[4]
                
        if base not in base_to_key:
        	continue
        
        peak_key = base_to_key[base]                
        
        pop_index = 0 # This is how many populations we have counted at this locus        

        # Deal with special case, 'REFERENCE', if the genome itself is assigned to one group or the other
        if "REFERENCE" in group0 or "REFERENCE" in group1:
        	vectors[peak_key][2][pop_index].append(1) # Reference is always defined as 1 (reference value)           
        	pop_index += 1

        # Go through all entries
        for entry in indices:
            call = line_data[entry]
            if call == "./.":
            	vectors[peak_key][2][pop_index].append('NaN')# If there's no data use non-numerical placeholder 
            	pop_index += 1
                continue # Skip over if there's no data

            # Get and store the genotype call
            GT = call.split(":")[0]            
            if GT.count("0") == 2: # If call is 2 reference alleles, it has value 1
            	vectors[peak_key][2][pop_index].append(1)
            	pop_index += 1
            elif GT.count("0") == 1: # If call is heterozygous, it has value 0.5
            	vectors[peak_key][2][pop_index].append(0.5)
            	pop_index += 1
            else: # If call is homozygous non-reference, it has value 0
                vectors[peak_key][2][pop_index].append(0)
            	pop_index += 1    
            	
# Now go through and compute the data
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
