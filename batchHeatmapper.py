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
#sort all files with the same regions file, same func. or sort one plot and sort all others based on first plot.
#determine the scaling for each plot based on matrix with largest number of regions (heatmaps only)
#	--long term: allow the passing of feature line width in px for heatmap rows
#
#call heatmapper or profiler and render plots (parallelize)
#put all plots into a new subdirectory

# 1.    Write a python script that executes heatmapper
# 2.    Write a python script that will execute heatmapper over a list of files
# 3.    Write a python script that will take arguments to pass to heatmapper
# 4.    Write a python script that will extract values from heatmapper before executing

    


import argparse

from deeptools import parserCommon ## contains all of the option flags for heatmapper and profiler
from deeptools import heatmapper    ## should be able to classes/methods to determine some parameters


