import random
import sys

random.seed()

n = int(sys.argv[1])
l = int(sys.argv[2])

for i in range(n):
    print(">%d" % (i + 1))
    print(''.join([random.choice('ACGT') for x in range(random.randint(l-10,l+10))]))
