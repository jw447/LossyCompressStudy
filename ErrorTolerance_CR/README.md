# Cluster 18'

## Insight of Error Tolerance and Compression Ratio

We think the error tolerance plays an role in compression ratio for lossy compression, higher error tolerance(larger error bound) should yield higher compression ratio. But according to experiment, the increasing of __error bound__ doesn't linearly conclude increasing of __compression ratio__. Thus we want to find out the relation between error tolerance and compression ratio, and further more explain the relation by the algorithm of SZ.

### Experiment Plan

Typical SZ compression uses error bound:

absErrBound = 0.000001

relBoundRatio = 1E-5

pw_realBoundRatio = 1E-5 //what we are using

1. We start with point-wise relative Bound Ratio. From 1E-10 to 1E-0, we collect the compression matrics(which includes the compression ratio).
