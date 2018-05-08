import sys
from operator import itemgetter
from math import log

data_file = sys.argv[1] # This should be final simulation output file, sorted in chr/pos order
min_peak = float(sys.argv[2]) # Minimum -log10 pvalue to consider calling a peak, recommend 5
shoulder_height = float(sys.argv[3]) # How many LOD can the values drop while we still call it a peak? Suggest 2
chr_distance = float(sys.argv[4]) # How far do we go from the last value to find the next? Recommend 20000?? Unclear

# Import data
snp_scores = {}
keys = {}
with open(data_file) as f1:
    for line in f1:
        line_data = line.strip().split()
        key = line_data[0] + line_data[1] + line_data[2] + '_' + line_data[3] # Unambiguous identifier
        chrom = line_data[0]
        pos = int(line_data[1])
        score = -log(float(line_data[-1]), 10)
        
        # If the score is too low for us to ever care about it, forget it
        if score < min_peak - shoulder_height:
        	continue

        # Now store the score
        if chrom not in snp_scores:
            snp_scores[chrom] = []
            keys[chrom] = {}
        snp_scores[chrom].append([pos, key, score])
        keys[chrom][key] = line_data

# Now go through and start calling peaks
for chromosome in snp_scores:
    	    
    windows = [] # 0 - center, 1 - current lower bound, 2 - current upper bound
    chr_data = snp_scores[chromosome]
    inpeaks = set()
    sorted_data = sorted(chr_data, key=itemgetter(2), reverse=True) # Sort by score

	# Start at the highest scoring SNP, work way down to seed new peaks
    for entry in sorted_data:
		position, key, score = entry[0], entry[1], entry[2]
	
		# Is this point contained in another peak window already?
		if key in inpeaks:
			continue	
		for window in windows:
			if window[1] <= position and position <= window[2]:
				inpeaks.add(key)
				continue
		if key in inpeaks:
			continue
	
		# If we're here, we need to start a new window, if the score is high enough
		if score < min_peak:
			break # The list is sorted, so if we fail once, we'll fail the rest of the time too		
		#windows.append([key, position, position, score]) # [Center key, lower bound, upper bound, chisq score]
		new_window = [key, position, position + 1, score]
		encroaching = False
		new_min_score = score - shoulder_height # We must be within 2 LOD of the peak

		# Extend newly seeded window to right as much as possible
		plus_sorted = sorted(chr_data)
		minus_sorted = sorted(chr_data, reverse=True)
		in_range = False
		for new_entry in plus_sorted:
			if new_entry[1] == key:
				in_range = True # Note when we pass our center point, it is the signal to start extending the window
				continue
			if in_range and new_entry[0] - new_window[2] <= chr_distance: # If we're to the right of the peak and close enough...				
				if new_entry[2] >= new_min_score: # And above the minimum score...
					new_window[2] = new_entry[0] + 1 # Move the edge					
					if new_entry[1] in inpeaks:
						encroaching = True # If we're encroaching somewhere, make a note
					inpeaks.add(new_entry[1]) # This point is now in a peak
					
			elif in_range:
				break # If not, the peak is over
		
		# Extend newly seeded window to the left as much as possible
		in_range = False
		for new_entry in minus_sorted:
			if new_entry[1] == key:
				in_range = True # Note when we pass our center point, it is the signal to start extending the window
				continue
			if in_range and new_window[1] - new_entry[0]  <= chr_distance: # If we're to the left of the peak and close enough...
				if new_entry[2] >= new_min_score: # And above the minimum score...

					new_window[1] = new_entry[0] # Move the edge					
					if new_entry[1] in inpeaks:
						encroaching = True						
					inpeaks.add(new_entry[1])  # This point is now in a peak
					
			elif in_range:
				break # If not, the peak is over
					
		# If we haven't encroached on another window, save this window		
		if not encroaching:
			windows.append(new_window)
			inpeaks.add(new_window[0])
			

    # Output data for this chromosome
    for peak_window in windows:    	

    	peak_window.append(peak_window[-1])
    	peak_window[-2] = peak_window[2] - peak_window[1] # Insert peak width before maximum score
    	print chromosome + '\t' + str(peak_window[1]) + '\t' + str(peak_window[2]) + '\t' + str(peak_window[3]) + '\t' + '\t'.join([str(x) for x in keys[chromosome][peak_window[0]][1:]]) 	

	

