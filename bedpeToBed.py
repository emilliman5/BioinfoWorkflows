#!/usr/bin/env python
import csv
import sys
import re

bedpe=open(sys.argv[1], "rb")
f = csv.reader(bedpe , dialect="excel-tab")
of = csv.writer(sys.stdout, dialect="excel-tab")

for line in f:
    #take care of the header
    if(line[0][0] == "#"):
        of.writerow(line)
        continue
        
    else :
        of.writerow([line[0],min(int(line[1]), int(line[2]),int(line[4]),int(line[5])), max(int(line[1]), int(line[2]),int(line[4]),int(line[5])), line[6], line[7], "+"])
        
bedpe.close()