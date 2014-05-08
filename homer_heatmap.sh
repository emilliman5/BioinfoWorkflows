for file in ~/GR_1hr_CHIP/dexGRVsuntrGR_1hr_peak_locs.bed ~/GR_8hr_Reddy/reddy_8hrGR.bed ~/GR_18hr_ChIP/Vdavis_GR_18hrDex_hg18
	do
	filename=$(basename $file .bed)
	for dir in ~/HomerTagDirectories/hg18/*/
		do
		dirname=$(basename $dir)
		annotatePeaks.pl "$file" hg18 -size 10000 -hist 10 -ghist -d "$dir" > output/heatmaps/${filename}_${dirname}_heatmap_10000_tss.txt 
	
	done
done