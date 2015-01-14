#python wrapper for deeptools heatmapper and profiler
#
#inputs:	matrix files (req)
#	        output file basename (req)
#	        sort function (regions file not needed because the matrices should be all different regions)
#	        plot size H & W (req)
#
#               pass through all other options
#
#
# 4.    Write a python script that will extract intensity values from heatmapper before executing


import argparse
import imp
import os
import multiprocessing
from deeptools import parserCommon ## contains many of the option flags for heatmapper and profiler
hmScript=imp.load_source('hmScript', '/home/millimanej/workspace/deepTools/bin/heatmapper')

parser=argparse.ArgumentParser(description="Batch process matrix files into heatmaps using deepTools heatmapper")
parser.add_argument("-f","--files", nargs='+', help="list of matrices to be processed in batch")
parser.add_argument("-zMx", type=float, help="Maximum value for heatmap intensities (for all matrices)", default=0)
parser.add_argument("--ext", choices=["png","pdf","eps","svg","emf"], default="png", help="Specifies filetype for heatmap image")
parser.add_argument("--prefix", default="", help="string to prepend outfile filenames")
parser.add_argument("--suffix", default="", help="string to append to outfile names")
parser.add_argument("-hh", action="store_true", help="Equivalent to: heatmapper --help")
parser.add_argument("-p",type=int, help="Number of processers to disrisbute heatmapper across", default=1)
parser.add_argument("-R", "--regions", help="Bed file to sort each heatmap by")

batch_args=parser.parse_args()

if batch_args.hh:
    args = hmScript.parseArguments(['--help'])
    hmScript.main(args)

for f in batch_args.files:

    outfile = os.path.splitext(f)[0]
    if batch_args.suffix:
        outfile=outfile+"."+batch_args.suffix

    if batch_args.prefix:
        outfile=batch_args.prefix+outfile

    outfile=outfile+"."+batch_args.ext

    args = hmScript.parseArguments(['-m', f,'-o', outfile])

    if batch_args.zMx:
        args.zMax=batch_args.zMx
        args.yMax=batch_args.zMx ## the profile atop the heatmaps should be on the same scale as the heatmap intenisites

    if batch_args.regions:
        args.
    
    hmScript.main(args)




