# This program takes as input:
# [1] a VCF file with variants called on all populations of interest and 
# [2] a config file detailing which populations belong to which groups - ideally, should be a single population in each group
# [3] the bed file of regions
# [4] the molecular clock calibration - how much divergence is expected per 1,000,000 years? Recommend 0.007, as rough estimate derived from Colosimo et al 2010

import sys
from scipy import stats
from math import log
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

age_calibration = float(sys.argv[4])

# Import regions of interest
regions = {}
with open(sys.argv[3]) as f0:
	for line in f0:
		line_data = line.strip().split()
		key = line_data[0] + ':' + line_data[1] + '-' + line_data[2] # This defines the region
		chrom, start, stop = line_data[0], int(line_data[1]), int(line_data[2])
		regions[key] = [chrom, start, stop] # Store all regions of interest

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

# Now go through the VCF file and analyze line by line
divergent_counts = {} # This is where we will store all the data
indices = []
group0_indices = []
group1_indices = []
pop_names = []
snp_count = 0
with open(sys.argv[1]) as f2:
    for line in f2:

        if line[:2] == "##":
            continue # Skip initial headers

        if line[0] == '#': # Analyze main header for populations
            for index, entry in enumerate(line.strip().split()):
                if entry in group0:
                    group0_indices.append(index)
                if entry in group1:
                    group1_indices.append(index)
            continue

        line_data = line.strip().split()
        chrom, bp = line_data[0], int(line_data[1])
        
        # See if this falls into a region of interest
        true_key = ''
        for key in regions:
        	rchrom, rstart, rstop = regions[key]
        	if rchrom != chrom:
        		continue
        	if rstart <= bp and bp <= rstop:
        		true_key = key
        		break
        if true_key == '':
        	continue
        
        # If we're here, it's a SNP we should potentially care about
        allele_counts = [0, 0, 0, 0] # G0 Ref, G0 alt, G1 ref, G1 alt
        

        
        pop_index = 0 # This is how many populations we have counted at this locus        

        # Deal with special case, 'REFERENCE', if the genome itself is assigned to one group or the other
        if "REFERENCE" in group0 or "REFERENCE" in group1:
        	vectors[peak_key][2][pop_index].append(1) # Reference is always defined as 1 (reference value)           
        	pop_index += 1

        # Now go through the calls for each population
        
        # Start with calls in group 1
        for entry in group0_indices:
            call = line_data[entry]
            if call == "./.":
                continue # Skip over if there's no data
            # Get the genotype call
            GT = call.split(":")[0]            
            if GT.count("0") == 2: # If call is 2 reference alleles
                allele_counts[0] += 2                
            elif GT.count("0") == 1: # If call is heterozygous
                allele_counts[0] += 1
                allele_counts[1] += 1                
            else: # If call is homozygous non-reference
                allele_counts[1] += 2                

        # Now look at entries in group 1
        for entry in group1_indices:
            call = line_data[entry]
            if call == "./.":
                continue # Skip over if there's no data
            # Get the genotype call
            GT = call.split(":")[0]            
            if GT.count("0") == 2: # If call is 2 reference alleles
                allele_counts[2] += 2                
            elif GT.count("0") == 1: # If call is heterozygous
                allele_counts[2] += 1
                allele_counts[3] += 1                
            else: # If call is homozygous non-reference
                allele_counts[3] += 2                
                
        # Now score the results
        true_dif_snp = False
        # Require unanminous support to score SNP as different (missing data okay)        
        if allele_counts[0] == allele_counts[0] + allele_counts[1] and allele_counts[3] == allele_counts[2] + allele_counts[3]:
        	true_dif_snp = True # G0 ref, G1 alt case
        if allele_counts[1] == allele_counts[0] + allele_counts[1] and allele_counts[2] == allele_counts[2] + allele_counts[3]:
        	true_dif_snp = True # G1 ref, G0 alt case
        
        # Store the results
        if true_key not in divergent_counts:
        	divergent_counts[true_key] = 0
        if true_dif_snp:
        	divergent_counts[true_key] += 1
                
# Output region calls in organized manner here
for key in regions:
	output = regions[key]
	if key not in divergent_counts:
		snp_count = 0
	else:
		snp_count = divergent_counts[key]
	output.append(output[2] - output[1])
	output.append(snp_count)
	divergence = float(snp_count)/(output[2] - output[1])
	output.append(1 - divergence)
	output.append(divergence / age_calibration * 1000000)
	print '\t'.join([str(x) for x in output])