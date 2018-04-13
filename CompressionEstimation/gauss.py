import numpy as np
from scipy.stats import norm
import sys
#from statistics import mode

data1 = []
with open(sys.argv[1]) as f:
    for line in f:
        data1.append(int(line))

print('data len', len(data1))

print(len(set(data1)))

mu, std = norm.fit(data1)
print('center, std fit', mu, std)
#print 'center point', mode(data1)
mu, std = np.mean(data1), np.std(data1)
print('mu, std', mu, std)

data = np.random.normal(loc=mu, scale=std, size=int(sys.argv[2]))
#print data
mu, std = norm.fit(data)
print(mu, std)

s = set()
for i in data:
    if int(i) in s:
        continue
    else:
        s.add(int(i))

print('s', len(s))
