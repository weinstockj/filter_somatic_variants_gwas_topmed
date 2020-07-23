#!/bin/bash
set -eou pipefail

gwas_result=$1
filter=$2
false_germline_variants=$3
ac_threshold=$4 # recommend 150
out=$5

echo "filtering now"
go run filter_gwas_results.go $1 $2 $3 | 
    cut -f1,2,4,5,6,7,10,11,13 | 
    awk -v threshold=$ac_threshold '{ if ($5 >= threshold) {print} }' | 
    bgzip -c > $5

echo "tabixing now"
tabix -S1 -s1 -b2 -e2 $5
