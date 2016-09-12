#!/usr/bin/env python2

# Lachlan Dryburgh 188607
# COMP90014 - Assignment 2 - Task 1
# 13/09/2016

from Bio import SeqIO
import sys

# DeBruijn Graph class
class DeBruijn:
  def __init__(self, reads, klength, cutoff):
    self.k = klength
    self.cut = cutoff
    self.kmers = self.kDict(reads)
    self.graph = self.graph()
  def __repr__(self):
    g_str = "\n"
    for k,v in self.graph.iteritems():
      g_str = g_str + str(k) + "->" + str(v) + "\n"
    return g_str     

# Construct kmer dictionary 
# I have choosen to increment bot the kmer and its reverse complement as 
# keeping both is useful for ensuring the algorithm is working properly
  def kDict(self, reads):
    kmers = {}
    for r in reads:
      for i in xrange(len(r)-self.k+1):
        kmer = r.seq[i:i+self.k]
        rc = str(kmer.reverse_complement())
        kmer = str(kmer)
        kmers[kmer] = kmers.get(kmer, 0) + 1
        kmers[rc] = kmers.get(rc, 0) + 1
    return kmers

# Construct a graph by linking k-1 overlapping kmers
# Cutoff is only applied in part 2    
  def graph(self):
    g = {}
    for k,v in self.kmers.iteritems():
      g[k] = []
      for l in ['A', 'C', 'G', 'T']:
        if k[1:]+l in self.kmers:
          g[k].append(k[1:]+l)
    return g            

# Parser
try:
  kvalue = int(sys.argv[2])
  cutoff = int(sys.argv[3])
  reads = SeqIO.parse(sys.argv[1], "fastq")
  graph = DeBruijn(reads, kvalue, cutoff) 
except IOError as e:
  print "Reqires 3 arguments fastq file, kmer length, coverage cutoff"
  print e
  sys.exit(1)

# print the DeBruijn graph, can redirect to file
print graph

