Compression Estimation
======================

INTRODUCTION
------------

This work focuses on the compression ratio estimation of SZ compressor. We first try to reproduce the figure 10 and TABLE III in IPDPS paper. Then we reproduce the SZ compression ratio estimation based on sample data compression metrics and quantization factors distribution.
 
We validate the compression ratio method using various datasets: Astro, Blast2_p, Bump, Dpot, Eddy, Fish, Sedov_p, YF17_p, YF17_t.

SETTINGS
--------

We try to compare our estimation method with Tao's method, where Tao used same quantization factors for each dataset between sample and full dataset and we used optimized quantization factors for both sample and full data.

RUN EXPERIMENTS
---------------

In order to run the experiments, one can simply run the command in terminals.

./cr_estimation.sh

In cr_estimation.sh you can decide which datasets shall be used in the experiment(variable _fname_). In the experiment, an point-wise relative error bound shall be decided explicitly (variable _pwr_).

Under the same directory there should be a inputdata directory where experiment data is listed.