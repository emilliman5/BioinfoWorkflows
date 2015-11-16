#
annotatePeaks.pl tss hg18 -size 10000 -noadj -hist 10 -d ~/HomerTagDirectories/hg18/* > hg18_tss_meta_5000.txt
annotatePeaks.pl tss hg18 -size 10000 -noadj -hist 10 -ghist -d ~/HomerTagDirectories/hg18/* > hg18_tss_meta__heatmap_5000.txt
annotatePeaks.pl tss hg18 -size 10000 -noadj -d ~/HomerTagDirectories/hg18/* > hg18_tss_homer_5000.txt
#
annotatePeaks.pl ~/GR_8hr_Reddy/reddy_etal_grpeakd.bed hg18 -size 10000 -noadj -hist 10 -d ~/HomerTagDirectories/hg18/* > hg18_reddyGR_meta_5000.txt
annotatePeaks.pl ~/GR_8hr_Reddy/reddy_etal_grpeakd.bed hg18 -size 10000 -noadj -hist 10 -ghist -d ~/HomerTagDirectories/hg18/* > hg18_reddyGR_meta_heatmap_5000.txt
#
annotatePeaks.pl ~/GR_18hr_ChIP/Vdavis_GR_18hrDex_hg18 hg18 -size 10000 -noadj -hist 10 -d ~/HomerTagDirectories/hg18/* > hg18_VCDGR18hr_meta_5000.txt
annotatePeaks.pl ~/GR_18hr_ChIP/Vdavis_GR_18hrDex_hg18 hg18 -size 10000 -noadj -hist 10 -ghist -d ~/HomerTagDirectories/hg18/* > hg18_VCDGR18hr_meta_heatmap_5000.txt
#
