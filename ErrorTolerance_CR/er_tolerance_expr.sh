#!/bin/bash

fname="astro blast2_p bump dpot eddy fish sedov_p yf17_p yf17_t"
#fname="astro bump yf17_p yf17_t"
# fname="blast2_p astro"

for f in $fname; do
  echo file is $f;
  mkdir $f;

  List="1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e-0"
  #List="1e-10 1e-5 1e-1"
  # List="1e-5"
  for i in $List; do
    echo The current point-wise errorbound is $i
    ../CompressionEstimation/pwr_expr_sz.sh -c sz -e $i -i "../CompressionEstimation/inputdata/sample_1/$f.dat" > "./$f/$f-s1-$i.log"
    ../CompressionEstimation/pwr_expr_sz.sh -c sz -e $i -i "../CompressionEstimation/inputdata/sample_10/$f.dat" > "./$f/$f-s10-$i.log"
    ../CompressionEstimation/pwr_expr_sz.sh -c sz -e $i -i "../CompressionEstimation/inputdata/original/$f.dat" > "./$f/$f-f-$i.log"
  done
done
