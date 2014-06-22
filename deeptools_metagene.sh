#!/usr/bin/env bash

cd ~/GR_binding_site_analysis/hg19/;

for cov in ~/GR_18hr_ChIP/dex_GR_18hr_T47D_hg19.bw ~/GR_1hr_ChIP/GR_1hr_Dex_hg19_sort_noOverlap.bw ~/GR_8hr_Reddy/GR_8h_Dex_coverage_VehSub_hg19.bw ~/GR_8hr_Reddy/GR_A549_ChIP_Dex_8hr_hg19_N1_merge.bw ~/GR_8hr_Reddy/GR_A549_ChIP_Veh_8hr_hg19_N1_merge.bw
	do

	coverage=${cov##*/}  ## Better than using an external binary like basename. removes path greedy match for all characters up to the last forward slash
	dir=${coverage%%.*}  ## Remove extension.

	mkdir -p "$dir" || exit 1;
	
	for i in ls *.bed;
		do 
		
		computeMatrix prefernce-point -R $i -S $cov -referencePoint center -b 5000 -a 5000 -out $dir/${i%%.*}.matrix -p 4;
		profiler -m ${i%%.*}.matrix --plotType se --plotWidth 16 --plotHeight 10 --refPointLabel GR --regionLabels   -out ${i%%.*}.png
		done;
	done;