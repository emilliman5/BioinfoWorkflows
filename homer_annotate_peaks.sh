for file in *.txt
do
filename=$(basename $file .txt)
	for dir in ~/HomerTagDirectories/hg18/*/
		do
		dirname=$(basename $dir)
		annotatePeaks.pl tss hg18 -list "$file" -size 10000 -hist 10 -ghist -d "$dir" > output/heatmaps/${filename}_${dirname}_heatmap_10000_tss.txt 
	
	done
done
