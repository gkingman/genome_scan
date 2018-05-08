import sys

# Takes as input:
# [1] the peak file to select from
# [2] a size (ie, 10000)
# [3] either '+' or '-' - do we want regions larger or smaller?

size = int(sys.argv[2])

with open(sys.argv[1]) as f1:
    for line in f1:
        ld = line.strip().split()
        reg_size = int(ld[3])
        if sys.argv[3] == '-' and reg_size <= size:
            print line.strip()
        if sys.argv[3] == '+' and reg_size >= size:
            print line.strip()
