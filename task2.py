#!/usr/bin/env python

# Lachlan Dryburgh 188607
# Comp90014 - Assignment 2 - Task 2 and 3
# 13/09/2016

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
    # construct the kmer dictionary
    kmers = {}
    for r in reads:
      for i in xrange(len(r)-self.k+1):
        kmer = r.seq[i:i+self.k]
        # Add both the kmer and its reverse complement to the dictionary
        rc = str(kmer.reverse_complement())
        kmer = str(kmer)
        kmers[kmer] = kmers.get(kmer, 0) + 1
        kmers[rc] = kmers.get(rc, 0) + 1
    return kmers
    
  def graph(self):
    g = {}
    for k,v in self.kmers.iteritems():
    # Apply the cutoff
      if v < self.cut:
        continue
      g[k] = []
      # Add links if k-1 overlapping kmer exist
      for l in ['A', 'C', 'G', 'T']:
        con = k[1:]+l
        # Shouldn't link to a kmer that is below the coverage cutoff
        if self.kmers.get(con, 0) >= self.cut:
          g[k].append(con)
    return g            

# Find kmers that are the start of contigs
  def starts(self):
    starts = set([])
    for k,v in self.graph.iteritems():
      if not len(v)==1:
        rc = str(Seq(k).reverse_complement())
        for s in v:
          #If a kmer is found that is linked to multiple kmers each of those links will start a new contig
          starts.add(s)
        if rc not in self.kmers:
          print rc + " wasn't found in kmers"
          sys.exit(1)
        #If there is any other number than 1 outgoing links then the rc is a start
        starts.add(rc)
    return starts

# Construct contigs
  def contigs(self):
    starts = self.starts()
    contigs = {}
    for s in starts:
      if not s in self.graph:
        print "start " + s + " not found"
        sys.exit(1)
      # Begin with the start kmers
      kmer = s
      contig = kmer
      links = self.graph.pop(kmer)
      # Extend the contig with linked kmers until the number of links isn't 1
      # or a new start kmer is reached 
      while len(links)==1 and not links[0] in starts:
        kmer = links[0]
        contig += kmer[-1]
        links = self.graph.pop(kmer)
      contigs[s] = (contig,links)      
    
    #Any kmers remaing in the graph must be part of unbranched loops
    while len(self.graph)>0:
      # start at an arbitary kmer
      start,link = self.graph.popitem()
      contig = start
      while len(link) == 1 and not link[0]==start:
        # traverse the looped contig until reaching the starting kmer again  
        contig += link[0][-1]
        link = self.graph.pop(link[0], [])
        
      contigs[start] = (contig,link)

    return contigs
    
# Parse file
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
totalkmers = sum(graph.kmers.values())
uniquekmers = len(graph.kmers)

clens = []
j = 0
p =[]
for c,l in contigs.iteritems():

# Make a pretty graph
  labels[c] = l[0]
  for edge in l[1]:
    G.add_edge(c,edge)
  clen = len(l[0])

# write contigs of length 50 or greater to file
# use p to prevent the printing and counting of reverse complements
  if clen>=50 and not str(Seq(l[0]).reverse_complement()) in p  :
    j += 1
    p.append(l[0])
    clens.append(clen)
     
    f.write('>c' + str(j) + '\n' + l[0] + '\n\n')

# Sanity check, should always find a reverse complement for a contig
  #rc = str(Seq(l[0]).reverse_complement())
  #if rc[:graph.k] in contigs:
  #  if not rc == contigs[rc[:graph.k]][0]:
  #    print l[0][:k] +  " didn't match the rc correctly"
  #else:
  #  print "No contig found starting with rc kmer " + rc
    
    
    
  #else:
  #  print 'match ' + l[0]

f.close()
 
clens.sort()


i = len(clens)/2
while sum(clens[:i])<((sum(clens)+1)/2):
  i += 1
          
print "Unique kmers " + str(uniquekmers)
print "No. of nodes: " + str(len(contigs))
print 'No of contigs: ' + str(len(clens)) + " (length 50 or greater)"
print 'Combined contig length: ' + str(sum(clens))
print 'Average contig length: ' + str(sum(clens)/len(clens))
print 'Longest contig: ' + str(max(clens))
print 'N50: ' + str(i) 

pos = nx.spring_layout(G)

nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_edges(G, pos, arrows = True)    
nx.draw_networkx_labels(G,pos, font_size=4)

plt.show()


