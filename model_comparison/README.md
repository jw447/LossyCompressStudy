MODEL COMPARISON
=================

INTRODUCTION
-----------------

This work focuses on the compression ratio estimation of SZ compressor. We first try to reproduce the figure 10 and TABLE III in IPDPS paper. Then we reproduce the SZ compression ratio estimation based on sample data compression metrics and quantization factors distribution.
 
We validate the compression ratio method using various datasets: fish, eddy, yf17....

SETTINGS
--------

__Quantization intervals__ (See paper for other datasets)

Fish: 10000

Eddy: 65535

YF17: 8192

__errorBoundMode__

PW_REL

__Sampling method__ 

 _chunk interval sampling_ and _chunk random sampling_ that can preserve the bounded locality which enables better compression estimation of full data.
 