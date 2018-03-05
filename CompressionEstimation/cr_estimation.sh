#!/bin/bash

# fname="astro blast2_p bump dpot eddy fish sedov_p yf17_p yf17_t"
#fname="astro bump yf17_p yf17_t"
fname="astro blast2_p bump dpot eddy fish sedov_p yf17_p yf17_t"

for f in $fname; do
  echo file is $f;
  mkdir $f;

  pwr="1e-5"
  for i in $pwr; do
    echo The current point-wise errorbound is $i
    ./pwr_expr_sz.sh -c sz -e $i -i "./inputdata/sample_1/$f.dat" > "./$f/$f-s1-$i.log"
    ./pwr_expr_sz.sh -c sz -e $i -i "./inputdata/sample_10/$f.dat" > "./$f/$f-s10-$i.log"
    ./pwr_expr_sz.sh -c sz -e $i -i "./inputdata/original/$f.dat" > "./$f/$f-f-$i.log"
  done
done
