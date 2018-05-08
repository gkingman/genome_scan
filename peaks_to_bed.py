import sys
from math import log

with open(sys.argv[1]) as f1:
    for line in f1:
        ld = line.strip().split()
        chrom = ld[0]
        center = ld[4]
        score = -log(float(ld[-1]), 10)
        score = str(round(score, 2))
        outkey = chrom + ':' + center + '|' + score
        output = ld[:3]
        output.append(outkey)
        print '\t'.join(output)
