#!/usr/bin/env python2

import sys
import os
import numpy
import StringIO
import numpy as np
import pprint as pprint
total_str=''
offsets_str=''

with open('triplets', 'r') as file:
    total_str= file.read()

triplets = np.genfromtxt(StringIO.StringIO(total_str), delimiter=' ', comments="#").astype(int)
base = triplets[:,0]
middle = triplets[:, 1]
head = triplets[:,2]

with open('offsets', 'r') as file:
    offsets_str = file.read()

offsets = np.genfromtxt(StringIO.StringIO(offsets_str), delimiter=' ', comments="#").astype(int)
offsets = list(offsets)
offsets += [len(base)]
print offsets
pairs = []
ranges1 = [(offsets[i], offsets[i+1]) for i in range(len(offsets)-1)]
ranges2 = [[ranges1[i], ranges1[j]] for i in range(len(ranges1) - 1) for j in range(i + 1, len(ranges1))]
print ranges1
pprint.pprint(ranges2)
pairs2 = []

for slbase in range(len(offsets) - 1):
    for slnext in range(slbase + 1, len(offsets) - 1):
        print offsets[slbase], offsets[slbase+1], offsets[slnext], offsets[slnext+1]
        print len(pairs2)
        for i in range(offsets[slbase], offsets[slbase+1]):
            for j in range(offsets[slnext], offsets[slnext+1]):
                if head[i] == base[j]:
                    pairs2 += [(i, j)]

#pprint.pprint(pairs2)
for i in range(0, offsets[-2]):
    for j in range(i + 1, len(base)):
        if head[i] == base[j]:
            pairs += [(i, j)]
print(pairs == pairs2)
print(len(pairs))
