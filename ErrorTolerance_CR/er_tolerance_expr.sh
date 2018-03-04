#!/bin/bash

# fname="astro blast2_p bump dpot eddy fish sedov_p yf17_p yf17_t"
#fname="astro bump yf17_p yf17_t"
fname="blast2_p astro"

for f in $fname; do
  echo file is $f;
  # mkdir $f;

  # List="1e-11 5e-11 1e-10 5e-10 1e-9 5e-9 1e-8 5e-8 1e-7 5e-7 1e-6 5e-6 1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 1e-1 5e-1 1e-0 5e-0"
  #List="1e-10 1e-5 1e-1"
  List="1e-5"
  for i in $List; do
    echo The current point-wise errorbound is $i
    ../model_comparison/pwr_expr_sz.sh -c sz -e $i -i "../model_comparison/inputdata/sample_1/$f.dat" > "./$f/$f-s1-$i.log"
    ../model_comparison/pwr_expr_sz.sh -c sz -e $i -i "../model_comparison/inputdata/sample_10/$f.dat" > "./$f/$f-s10-$i.log"
    ../model_comparison/pwr_expr_sz.sh -c sz -e $i -i "../model_comparison/inputdata/original/$f.dat" > "./$f/$f-f-$i.log"
  done
done
