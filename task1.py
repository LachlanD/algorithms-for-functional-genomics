#!/usr/bin/env python2

from Bio import SeqIO
import sys

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

  def kDict(self, reads):
    kmers = {}
    for r in reads:
      for i in xrange(len(r)-self.k+1):
        kmer = r.seq[i:i+self.k]
        rc = str(kmer.reverse_complement())
        kmer = str(kmer)
        if kmer not in kmers:
          kmers[kmer] = 1
        else:
          kmers[kmer] += 1
        if rc not in kmers:
          kmers[rc] = 1
        else:
          kmers[rc] += 1
    return kmers
    
  def graph(self):
    g = {}
    for k,v in self.kmers.iteritems():
      g[k] = []
      for l in ['A', 'C', 'G', 'T']:
        if k[1:]+l in self.kmers:
          g[k].append(k[1:]+l)
    return g            

try:
  kvalue = int(sys.argv[2])
  cutoff = int(sys.argv[3])
  reads = SeqIO.parse(sys.argv[1], "fastq")
  graph = DeBruijn(reads, kvalue, cutoff) 
except IOError as e:
  print "Reqires 3 arguments fastq file, kmer length, coverage cutoff"
  print e
  sys.exit(1)

print graph

