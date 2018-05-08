import sys

with open(sys.argv[1]) as f1:
    for line in f1:
        ld = line.strip().split()
        chrom, start, stop = ld[0], int(ld[1]) - int(sys.argv[2]), int(ld[1]) + int(sys.argv[2])
        print '\t'.join([str(x) for x in [chrom, start, stop]])
