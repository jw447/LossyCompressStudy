#!/bin/bash

# fname="astro blast2_p bump dpot eddy fish sedov_p yf17_p yf17_t"
#fname="astro bump yf17_p yf17_t"
# fname="astro blast2_p bump dpot eddy fish sedov_p yf17_p yf17_t"
fname="fish"
for f in $fname; do
  echo file is $f;
  # mkdir $f;

  pwr="1e-5"
  for i in $pwr; do
    echo The current point-wise errorbound is $i
    # ./pwr_expr_sz.sh -c sz -e $i -i "./inputdata/sample_1/$f.dat" > "./$f/$f-s1.log"
    ./pwr_expr_sz.sh -c sz -e $i -i "./inputdata/sample_10/$f.dat" > "$f-s10.log"
    ./pwr_expr_sz.sh -c sz -e $i -i "./inputdata/original/$f.dat" > "$f-f.log"
  done
done
