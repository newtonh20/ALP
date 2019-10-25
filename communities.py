__author__ = 'Newton Campbell'
from networkx import *
from community import *

G=nx.erdos_renyi_graph(100, 0.01)
dendo = generate_dendogram(G)
for level in range(len(dendo) - 1) :
    print "partition at level", level, "is", partition_at_level(dendo, level)