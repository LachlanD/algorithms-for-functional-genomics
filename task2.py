#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import networkx as nx
import matplotlib.pyplot as plt
import copy

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
      if v < self.cut:
        continue
      g[k] = []
      for l in ['A', 'C', 'G', 'T']:
        con = k[1:]+l 
        if con in self.kmers and self.kmers[con] >= self.cut:
          g[k].append(k[1:]+l)
    return g            
  
  def starts(self):
    starts = []
    graph = copy.deepcopy(self.graph)
    while len(graph) > 0:
      k,v = graph.popitem()    
      kmer = k
      links = v
      while len(links)==1:
        kmer = links[0]
        links = graph.pop(kmer, [])
    

      rc = str(Seq(kmer).reverse_complement())
      if not rc in starts:
        starts.append(rc)
    return starts

  def contigs(self):
    starts = self.starts()
    contigs = {}
    graph = copy.deepcopy(self.graph)
    for s in starts:
      contig = s
      links = graph.pop(s,[])
      while len(links)==1:
        contig += links[0][-1]
        links = graph.pop(links[0],[])
      contigs[s] = (contig,links)
    return contigs
    

try:
  kvalue = int(sys.argv[2])
  cutoff = int(sys.argv[3])
  reads = SeqIO.parse(sys.argv[1], "fastq")
  graph = DeBruijn(reads, kvalue, cutoff) 
except IOError as e:
  print "Reqires 3 arguments fastq file, kmer length, coverage cutoff"
  print e
  sys.exit(1)

f = open('reads.fa', 'w')
G = nx.DiGraph()
labels = {}
contigs = graph.contigs()
nkmers = len(graph.kmers)/2
clens = []
j = 0
for c,l in contigs.iteritems():
  labels[c] = l[0]
  for edge in l[1]:
    G.add_edge(c,edge)
  clen = len(l[0])
  if clen>=50:
    j += 1
    clens.append(clen)
    f.write('>c' + str(j) + '\n' + l[0] + '\n\n')

f.close()
 
clens.sort()

i = len(clens)/2
while sum(clens[:i])<((sum(clens)+1)/2):
  i += 1
          
print 'No. kmers: ' + str(nkmers)
print 'No. of nodes: ' + str(len(contigs))
print 'No of contigs: ' + str(len(clens))
print 'Combined contig length: ' + str(sum(clens))
print 'Average contig length: ' + str(sum(clens)/len(clens))
print 'Longest contig: ' + str(max(clens))
print 'N50: ' + str(i)

pos = nx.spring_layout(G)

nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_edges(G, pos, arrows = True)    
nx.draw_networkx_labels(G,pos, font_size=4)

plt.show()  
