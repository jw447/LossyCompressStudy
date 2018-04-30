# Understanding and Extrapolating Compression Ratios of HPC Scientific Data

Scientific simulations on high-performance computing (HPC) systems generate vast amounts of floating-point data that need to be reduced in order to lower the I/O cost. Lossy compressors trade accuracy for performance, and were demonstrated to be effective in reducing data volume. However, a key hurdle to a wide adoption of lossy compressors is that the reduction performance is not well understood, in particular, the trade-off between accuracy and the associated performance. The consequence is that domain scientists often need to exhaust all error bounds before they can figure the most cost-effective setup. As a result, the current practice of using lossy compressions to compress datasets is through trial and error, which is not efficient for large datasets.

This paper aims to analyze, understand, and estimate the reduction performance. In particular, we aim to extrapolate the compression ratios of SZ, a state-of-the-art lossy compressor that achieves superior performance, at various error bounds, based upon the compressor-internal metrics collected under a given error bound (i.e., base error bound). We propose a SZ compression ratio extraplolation schme based on the similarity among quantization factor distributions across a range of error bounds, and the fact that they approximately can be characterized using Gaussian distribution. We evaluate the performance extrapolation scheme using nine real HPC datasets, and the results confirm the efficacy of our approach.

## Getting Started

```
Please make sure SZ is configured properly and Python3 is installed for data analysis.
```

# Compression Ratio Estimation

1. baisc SZ Compression

./run.sh -c sz -i testdouble\_8\_8\_128.dat

2. Compression at certain error bounds

```
cd ErrorTolerance_CR

# Dataset names and error bounds can be set in this script.
./er_tolerance_expr.sh

```
make sure following data is captured:

1) original data
2) compression metrics listed in paper
3) curve-fitting prediction error

3. Result analysis

```
cd output/cr_estimation_cfg/cr_estimation

```

  - HitRatio_Estimate.py: function for hit ratio estimation
  - NodeCount_Estimation.py: function for nodecount estimation
  - CompressionRatio_esti_cfg.ipynb: conduct compression ratio estimation.
