#!/bin/bash
List="dpot.dat astro.dat fish.dat sedov_p.dat blast2_p.dat eddy.dat yf17_p.dat yf17_t.dat bump.dat"

while getopts ':m:' flag; do
	case ${flag} in
		m) sample_mode=$OPTARG;;
	esac
done

echo "Sample mode is $sample_mode"

if [ "$sample_mode" == "1" ];
then
	file_path="sample_1/"
elif [ "$sample_mode" == "10" ];
then 
	file_path="sample_10/"
elif [ "$sample_mode" == "100" ];  
then
	file_path="original/"
fi

set -- $List

for i;do
	echo "$i"
		./sz.sh -c sz -i "./inputdata/$file_path/$i" > "./outputdata/$i.log"
done
