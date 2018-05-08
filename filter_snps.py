import sys
from math import log

with open(sys.argv[1]) as f1:
    for line in f1:
        ld = line.strip().split()
        score = -log(float(ld[-1]), 10)
        if score >= float(sys.argv[2]):
            print line.strip()
