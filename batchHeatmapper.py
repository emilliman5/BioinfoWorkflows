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

# 2.    Write a python script that will execute heatmapper over a list of files
            
# 4.    Write a python script that will extract values from heatmapper before executing

    
import argparse
import imp
import os
from deeptools import parserCommon ## contains all of the option flags for heatmapper and profiler
from deeptools import heatmapper    ## should be able to classes/methods to determine some parameters
hmScript=imp.load_source('hmScript', '/home/millimanej/workspace/deepTools/bin/heatmapper')

parser=argparse.ArgumentParser(description="Process matrix files into heatmaps in batch")
parser.add_argument("-f","--files", nargs='+', help="list of matrices to be processed in batch")
parser.add_argument("-zMx", type=float, help="Maximum value for heatmap intensities (for all matrices)", default=0)
parser.add_argument("--ext", choices=["png","pdf","eps","svg","emf"], default="png")
parser.add_argument("--prefix", default="", help="prefix for outfiles")
parser.add_argument("--suffix", default="", help="suffix for outfiles (not the extension")
batch_args=parser.parse_args()

if batch_args.zMx == 0:
    #retrieve largest z-value (heatmap color intensity) from list of files and set zMax to that.
    #this will not be trivial... zMax is determined during the plot function call. heatmapper will
    #need to be tweaked to accomplish this
    batch_args.zMx=""   

print batch_args.files

for f in batch_args.files:
    
    outfile = os.path.splitext(f)[0]
    if batch_args.suffix:    
        outfile=outfile+"."+batch_args.suffix
    
    if batch_args.prefix:
        outfile=batch_args.prefix+outfile

    outfile=outfile+"."+batch_args.ext
    print "Matrix file: "+outfile
    
    args = hmScript.parseArguments(['-m', f,'-o', outfile])
    
    if batch_args.zMx:
        args.zMax=batch_args.zMx
   
    hmScript.main(args)




