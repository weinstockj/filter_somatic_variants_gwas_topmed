#!/bin/bash
set -eou pipefail

gwas_result=$1
filter=$2
false_germline_variants=$3
ac_threshold=$4 # recommend 150
out=$5

echo "filtering now"
# go run filter_gwas_results.go $1 $2 $3 | 
binary=/net/topmed2/working/jweinstk/count_singletons/new_drivers/filter_somatic_variants_from_gwas/filter_gwas_results
$binary $1 $2 $3 | 
    cut -f1,2,4,5,6,7,9,10,11,13 | 
    awk -v threshold=$ac_threshold 'NR==1;{ if ($5 >= threshold && ($7-$5) >= threshold) {print} }' | 
    bgzip -c > $5

echo "tabixing now"
tabix -S1 -s1 -b2 -e2 $5
