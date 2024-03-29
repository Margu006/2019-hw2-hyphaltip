#!/usr/bin/env python3

import os,csv,gzip,re

# write a script to create new file called closed.txt

outfile="closed.txt"
# on the HPCC cluster you can use the file
file="/bigdata/gen220/shared/simple/title.basics.tsv.gz"
# or if you run this on your own computer use this
#file="title.basics.tsv.gz"


if not os.path.exists(file):
    os.system("curl -O https://datasets.imdbws.com/title.basics.tsv.gz")

door_count = 0
door_word_count = 0
opened = 0
closed = 0
open_or_closed = 0
open_and_closed = 0

openre = re.compile(r'\sOpen\s')
closedre = re.compile(r'\sClosed\s')
#openclosedre = re.compile(r'\s(Open|Closed)\s')

with gzip.open(file,"r") as fh:
    # now use the csv reader to read in the file, delimiter is '\t'
    header = next(fh)
    for line in fh:
        row = line.decode('utf-8').strip().split("\t")
        title = row[2]

        if "door" in title or "Door" in title:
            door_count += 1

        if " door " in title or " Door " in title:
            door_word_count += 1

        o = openre.search(title)
        c = closedre.search(title)
        if o:
            opened +=1

        if c:
            closed +=1

        if o or c:
            open_or_closed += 1

        if o and c:
            open_and_closed += 1

print("There are %d titles with door or Door"%(door_count))
print("There are %d titles with [Ddoor] as word"%(door_word_count))
print("There are %d titles with Open as a word and %d with Closed"%(
    opened,closed))

print("There are %d titles with Open or Closed"%(
    open_or_closed))

print("There are %d titles with Open and Closed"%(
    open_and_closed))
