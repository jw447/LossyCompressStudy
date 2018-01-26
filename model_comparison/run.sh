#!/bin/bash
search_dir="./inputdata/full"

for i in $(ls $search_dir); do

	if [[ $i == "astro_"*".dat" ]]; then
		echo "$i"
		./sz.sh -c sz -i "$search_dir/$i" > "./outputdata/$i.log"
	fi

done
