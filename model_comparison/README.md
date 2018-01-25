MODEL COMPARISON
=================

INTRODUCTION
-----------------

This folder contains the work to predict SZ compression rate on different data. Data comes from three original dataset: Fish, Eddy and YF17. Each dataset is sampled using _block interval sampling_, _block random sampling_, _chunk interval sampling_, _chunk random sampling_, _point interval sampling_, _point random sampling_, _prefix sampling_. We calculate and print out the quantizaton factor for all data under all settings. Hope we can find something useful here.

SETTINGS
--------

__Quantization_intervals__

Fish: 10000
Eddy: 65535
YF17: 8192


__errorBoundMode__

PW_REL


INPUTDATA
---------

small: contains eddy, fish and yf17

full: 

original: uncompressed data



