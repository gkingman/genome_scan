import sys

new_min_width = int(sys.argv[2])

with open(sys.argv[1]) as f1:
    for line in f1:
        ld = line.strip().split()

        chrom, start, stop = ld[0], int(ld[1]), int(ld[2])
        width = stop - start
        if width < new_min_width:
            to_widen = (new_min_width - width) / 2
            start = start - to_widen
            stop = stop + to_widen
        print '\t'.join([str(x) for x in [chrom, start, stop]])
            
