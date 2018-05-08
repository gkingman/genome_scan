# This program takes as input:
#[1] a VCF file with variants called on all populations of interest and 
#[2] a config file detailing which populations belong to which groups
#[3] an indication for which group to use to generate the custom ref

import sys
from scipy.stats import fisher_exact
from math import log
from math import factorial

def choose(n, x):
	return factorial(n) / factorial (x) / factorial(n - x)

def explicit_p_calculation(hr, het, ha, g0nnsc, g0r):
	# Takes as input:
	# [1] Number homozygous reference samples
	# [2] Number heterozygous samples
	# [3] Number homozygous alternate samples
	# [4] Total number of non-null samples in Group 0
	# [5] Total number of reference alleles called in group 0
	tot_calls = hr + het + ha # This is the total number of non-null samples

	# Enumerate the bounds of possibilities
	# Triangle inequality 1: G0 cannot have more than [2] + [3] calls not in [1]
	# Triangle inequality 2: G0 cannot have more than min([1], [4]) calls in [1]
	max_out = het + ha # We can pick at most het + ha samples that aren't reference
	min_in = max(g0nnsc - max_out, 0) # Therefore, we must have at least this many reference calls
	max_in = min(hr, g0nnsc) # We cannot have more g0 reference calls than total reference calls or g0 calls
	
	# i is the number of reference reads in group 0
	tot = 0.0
	probs = {}
	for i in xrange(min_in, max_in + 1):
		min_het = max(g0nnsc - ha - i, 0) # We must be able to fill the remaining calls with the ha calls
		max_het = min(het, g0nnsc - i) # We cannot have more g0 het calls than total het calls or remaining calls 
		for j in xrange(min_het, max_het + 1): # Consider all possibilities for het calls
			k = g0nnsc - i - j # This is the necessary amount of ha calls in g0
			
			p = choose(hr, i) * choose(het, j) * choose(ha, k) / float(choose(tot_calls, g0nnsc)) # Compute the probability of this i, j, and k		
			ref_count = str(2 * i + j) # Count total number of reference alleles in group 0 in this case
			if ref_count not in probs:
				probs[ref_count] = 0.
			probs[ref_count] += p
				
	counts = list(probs)
	counts.sort(key=int)
	ps = [0., 0., 0.] # Sum the probability of being greater, equal, or less than the observed g0 reference allele count
	for count in counts:
		if int(count) > g0r:
			ps[0] += probs[count]
		elif int(count) < g0r:
			ps[2] += probs[count]
		else:
			ps[1] += probs[count]
	
	# Calculate one-sided p-value
	one_sided_p = ps[1] + min(ps[0], ps[2])
	
	# Convert dict to list for easier operations
	list_probs = []
	for key in probs:
		list_probs.append([int(key), probs[key]])
	list_probs.sort()	
	
	# Now calculate the other tail
	max_score = max(list_probs, key=lambda x: x[1]) # What is the most likely reference allele count in g0? Use this to split the distribution in half
	tail2_p = 0.
	if max_score[0] > g0r: # Case 1: we are in the lower part of the curve and our other tail is the higher part
		for score in list_probs: 
			if score [0] >= max_score[0] and score[1] <= ps[1]: # Iterate over the higher part of the curve and if this position is more extreme than our observed statistic
				tail2_p += score[1] # Increase our sum
	else: # Case 2: we are in the higher part of the curve and our other tail is the lower part
		for score in list_probs: 
			if score [0] <= max_score[0] and score[1] <= ps[1]: # Iterate over the higher part of the curve and if this position is more extreme than our observed statistic
				tail2_p += score[1] # Increase our sum
	
	two_sided_p = min(one_sided_p + tail2_p, 1)
	
	return two_sided_p

# Import population groups
group0 = []
group1 = []
with open(sys.argv[2]) as f1:
    for line in f1:
        line_data = line.strip().split()
        if len(line_data) != 2:
            continue
        if line_data[1] == '0':
            group0.append(line_data[0])
        else:
            group1.append(line_data[0])

# Now go through the VCF file and analyze line by line
group0_indices = []
group1_indices = []
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
        
        # Initialize counts of each allele
        allele_counts = [0, 0, 0, 0] # [0] Group 0 Ref, [1] Group 0 Alt, [2] Group 1 Ref, [3] Group 1 Alt
        samples = [0, 0, 0] # Homozygous reference, heterozygous, homozygous alternate

        # Deal with special case, 'REFERENCE', if the genome itself is assigned to one group or the other
        if "REFERENCE" in group0:
            allele_counts[0] += 2
        elif "REFERENCE" in group1:
            allele_counts[2] += 2

        # Flag any funny cases where the genotype is not the first data field
        data_format = line_data[8] # This is where VCF stores the field format


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
                samples[0] += 1
            elif GT.count("0") == 1: # If call is heterozygous
                allele_counts[0] += 1
                allele_counts[1] += 1
                samples[1] += 1
            else: # If call is homozygous non-reference
                allele_counts[1] += 2
                samples[2] += 1

        # Now look at entries in group 1
        for entry in group1_indices:
            call = line_data[entry]
            if call == "./.":
                continue # Skip over if there's no data

            # Get the genotype call
            GT = call.split(":")[0]
            
            if GT.count("0") == 2: # If call is 2 reference alleles
                allele_counts[2] += 2
                samples[0] += 1
            elif GT.count("0") == 1: # If call is heterozygous
                allele_counts[2] += 1
                allele_counts[3] += 1
                samples[1] += 1
            else: # If call is homozygous non-reference
                allele_counts[3] += 2
                samples[2] += 1

        # Filter out SNPs that cannot be informative
        # First, a filter for statistical power - skip entries with fewer than 10% of calls for either allele, and with fewer than 10 calls total
        ref_counts = allele_counts[0] + allele_counts[2]
        alt_counts = allele_counts[1] + allele_counts[3]
        if ref_counts > alt_counts * 9 or alt_counts > ref_counts * 9:
            continue # This requirement is subject to additional consideration
        if sum(allele_counts) < 10:
        	continue # Require a total of at least 10 calls

        # Skip entries without at least three calls from each group - no conclusions can be made on so few calls
        if allele_counts[0] + allele_counts[1] < 4 or allele_counts[2] + allele_counts[3] < 4:
        	continue

        # For comparison with previous results, calculate chi_sq - this is not used functionally
        g0_total = allele_counts[0] + allele_counts[1]
        g1_total = allele_counts[2] + allele_counts[3]
        e0 = float(ref_counts) * float(g0_total)/float(g0_total + g1_total)
        e1 = float(alt_counts) * float(g0_total)/float(g0_total + g1_total)
        e2 = float(ref_counts) * float(g1_total)/float(g0_total + g1_total)
        e3 = float(alt_counts) * float(g1_total)/float(g0_total + g1_total)
        expected = [e0, e1, e2, e3]	
        chi_sq_val = 0
        for i in xrange(4):
        	chi_sq_val += (allele_counts[i] - expected[i]) ** 2 / expected[i]

        # Now calculate p-value
        p_value = explicit_p_calculation(samples[0], samples[1], samples[2], (allele_counts[0] + allele_counts[1])/2, allele_counts[0])
        output = line_data[:2]
        output.extend(line_data[3:5])
        output.extend([str(x) for x in allele_counts])
        output.extend([str(x) for x in samples])
        output.append(str(chi_sq_val))
        output.append(str(p_value))
        print "\t".join(output)        
