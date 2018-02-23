#!/bin/bash

fname="astro_o.dat"

List="1e-10"

../model_comparison/sz.sh -c sz -i "../model_comparison/inputdata/full/$fname" > "../ErrorTolerance_CR/expr/$fname-$List.log"
