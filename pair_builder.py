import numpy as np
connectivity = [5, 7, 9]
prefixsum = list(np.cumsum(connectivity))
prefixsum = [0] + prefixsum[:-1] + [prefixsum[-2] + connectivity[-1]]

base = [0 for i in range(1024)]
print connectivity
print prefixsum
for gid in range(1024):
    where_do_I_belong = 0
    for i in range(len(prefixsum) - 1):
        where_do_I_belong += i * (prefixsum[i] <=gid) * (prefixsum[i+1] > gid)
    base[gid] = where_do_I_belong
print base
