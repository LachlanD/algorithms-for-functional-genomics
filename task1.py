#!/usr/bin/env python2

from Bio import SeqIO
import sys

try:
  kvalue = sys.argv[2]
  cutoff = sys.argv
except IOError as e:
  print "Reqires 3 arguments fastq file, kmer length, coverage cutoff"
  print e
  sys.exit(1)

try:
  for read in SeqIO.parse(sys.argv[1], "fastq"):
    print read  
except IOError as e:
  print "Could not read fastq input file (see below)!"
  print e
  sys.exit(1)


