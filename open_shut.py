#!/usr/bin/env python

import os,csv,gzip

# write a script to create new file called closed.txt

outfile="closed.txt"
# on the HPCC cluster you can use the file
file="/bigdata/gen220/shared/simple/title.basics.tsv.gz"
# or if you run this on your own computer use this
#file="title.basics.tsv.gz"


if not os.path.exists(file):
    os.system("curl -O https://datasets.imdbws.com/title.akas.tsv.gz")

with gzip.open(file,"r") as fh:
    # now use the csv reader to read in the file, delimiter is '\t'
