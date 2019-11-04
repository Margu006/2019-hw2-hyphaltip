#!/usr/bin/env python3

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

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

codon_table_Salmonella = {}
codon_table_Mycobacterium = {}

gene_counts = 0
gene_lens   = 0
GC          = 0

with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)

    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        # count the number of genes by adding 1 when we see one
        gene_counts += 1
        # total length based on the added lengths of all sequences
        gene_lens += len(seqstring)

        # iterate through each character in the string to see if it
        # is a G or C to keep track of these
        for ch in seqstring:
            if ch == "G" or ch == "C" or ch == "c" or ch == "g":
                GC += 1
        # iterate through the codons by going by 3s across the
        # sequence starting at 0 and going up to but not including
        # the length of the string (len(seqstring)
        for i in range(0,len(seqstring),3):
            codon = seqstring[i:i+3] # a codon is 3 bases extracted from string
            # if we have seen this codon aleady, just add one to its value
            # in the dictionary
            if codon in codon_table_Salmonella:
                codon_table_Salmonella[codon] += 1
            else:
                # otherwise set the value to 1
                codon_table_Salmonella[codon] = 1

# after we have seen all the codons, normalize the frequency
# by dividing by the total number of codons observed
for codon in codon_table_Salmonella:
    codon_table_Salmonella[codon] /= (gene_lens / 3)

# GC is the average across the total length of sequences
# so divide G/C count by total length
GC /= gene_lens

print("%d genes in %s"%(gene_counts,file1))
print("%d total gene length"%(gene_lens))
print("%.2f%% GC"%(100 * GC))
print("%d codon count"%(gene_lens/3))

print("===")

gene_counts = 0
gene_lens   = 0
GC          = 0

with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    for seq in seqs:
        seqname  = seq[0]
        seqstring= seq[1]
        # count the number of genes by adding one when we see on1
        gene_counts += 1
        gene_lens += len(seqstring)
        for ch in seqstring:
            if ch == "G" or ch == "C" or ch == "c" or ch == "g":
                GC += 1
        for i in range(0,len(seqstring),3):
            codon = seqstring[i:i+3]
            if codon in codon_table_Mycobacterium:
                codon_table_Mycobacterium[codon] += 1
            else:
                codon_table_Mycobacterium[codon] = 1

# normalize the frequency of codons by the total number of codons
for codon in codon_table_Mycobacterium:
    codon_table_Mycobacterium[codon] /= (gene_lens / 3)

# Make GC and average across all gene lengths
GC /= gene_lens

print("%d genes in %s"%(gene_counts,file2))
print("%d total gene length"%(gene_lens))
print("%.2f%% GC"%(100 * GC))# convert to a percentage, x by 100
print("%d codon count"%(gene_lens/3))
print("===")

# turn the filenames into a shorter name by removing everything after
# the first '.' for simpler output

print("\t".join(["Codon",file1.split(".")[0],file2.split(".")[0]]))
for codon in sorted(codon_table_Salmonella.keys()):
    # convert to a percentage x by 100
    codonsp1 = "%.2f%%"%(100 * codon_table_Salmonella[codon])
    codonsp2 = "%.2f%%"%(100 * codon_table_Mycobacterium[codon])

    print("\t".join([codon,codonsp1,codonsp2]))
