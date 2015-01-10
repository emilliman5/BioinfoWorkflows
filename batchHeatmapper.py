#python wrapper for deeptools heatmapper and profiler
#
#inputs:	matrix files (req)
#	output file basename (req)
#	sorted regions list or sort func (min, max, sd, average, median, etc)
#	plot size H & W (req)
#
#pass through all other options
#	
#determine the Zmax value and apply to all plots
#sort all files with the same list if specified or specified func.
#determine the scaling for each plot based on matrix with largest number of regions (heatmaps only)
#	--long term: allow the passing of feature line width in px for heatmaps
#
#call heatmapper or profiler and render plots (parallelize)
#put all plots into a new subdirectory