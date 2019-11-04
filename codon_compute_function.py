#!/usr/bin/env python

import os, gzip, itertools

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

# this is a function to summarize the data from each file
def codon_compute(fh):
    # initialize some variables we will use and also return values back
    codons = {}
    genes = 0
    gene_length = 0
    gc = 0
    # parse the file handle to get all the sequences back in a dictionary
    seqs = aspairs(fh)

    # iterate through all the sequences
    for seq in seqs:
        genes += 1 # everytime we see a sequence it is another gene, add up
        # the total length of gene sequences are the
        # sum of all the sequence lengths
        gene_length += len(seq[1])
        # to count the number of G and C just iterate through every
        # character in the sequence string which is stored in seq[1]
        # remember seq[0] is the sequence ID
        for char in seq[1]:
            # could use regex or can just do this simply with a
            # set of 4 tests
            if char == "G" or char == "C" or char == "g" or char == "c":
                gc += 1 # add to the GC counter if we see G or C

        # process each codon use range(0,end,3) to count by 3s to the end
        for i in range(0,len(seq[1]),3):
            # the codon is the start to start + 3
            codon = seq[1][i:i+3]
            # use this code to count the codons storing in a dictionary
            # if we have seen the codon before add one to the value
            if codon in codons:
                codons[codon] += 1
            else:
                # if we have not seen this codon before set the value to 1
                codons[codon]  = 1
    # turn gc into an average by dividing by the total length of sequence
    # we have examined
    gc /= gene_length

    # turn the codon counts into fraction by dividing by the
    # the total number of codons in this organism
    for codon in codons:
        codons[codon] /= (gene_length / 3)
    # return a list of 4 values:
    # the gene number, gene length, average GC fraction, and
    # the codon table dictionary
    return(genes,gene_length,gc,codons)


# now all the functions are defined we are going to start the main code

# these are the URLs of two organisms coding sequence files

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"


# download the data files
if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

# initialize some empty dictionaries for storing information
codon_tables = {}
codons  = {}
# turn this into a loop - for each of the files we provide
for file in (file1,file2):
    # turn the filename into a shorter name by removing everything after
    # the first '.' for simpler output
    name = file.split(".")[0]

    # open the file, which is gzip compressed, the "rt" is needed for
    # python 3 especially
    with gzip.open(file,"rt") as fh:
        # call the function codon_compute on this file handle
        # this returns a list which has 4 entries
        res = codon_compute(fh)

        # get the codon table from this result, the last itenm in the list
        codon_tables[name] = res[3]

        # we need to know what are all the codons by keeping track
        # of a unique list of all possible codons seen
        # just assign a value of 1 to all seen codons
        # when this is done will always have the complete set
        for codon in res[3]:
            codons[codon] = 1

        # print out the summary values for each species
        print("species is ",file.split(".")[0])
        print("gene count        = %d"%(res[0]))
        print("total gene length = %d"%(res[1]))
        print("GC                = %.1f%%"%(100 * res[2]))
        print("Codon count       = %d"%(res[1] / 3))
        print("==")

# now we want to print out a table of the results
# first we make a header row which will be tab separated
# Codon Sp1 Sp2
header = ["Codon"]
header.extend(list(codon_tables.keys()))
print("\t".join(header))

# now lets print out the codons, sort the codons alphabetically
for codon in sorted(codons.keys()):
    # make a row for printing out, first column will be the codon strin
    row = [codon]
    # iterate through each of codon tables in the data array we made above
    for r in codon_tables:
        # add to this row the % value ( by multiplying by 100)
        # use %.4f to print out the formatted floating point string
        # with 4 numbers after the period
        row.append("%.4f"%(100*codon_tables[r][codon]))
    # now print the row, using the 'join' function where the separator is
    # tab ("\t")
    # could have also done this with csvwriter.writerow(row)
    print("\t".join(row))
