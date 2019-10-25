import igraph
import networkx as nx
import numpy as np

def igraph_to_networkx(igraph_g):
    A = igraph_g.get_adjacency()
    A = np.matrix(A.data)
    return nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
def networkx_to_igraph(G):
    new = igraph.Graph()
    new.add_vertices(G.nodes())
    edges = G.edges()
    new.add_edges(edges)

    weights = G.edges(data=True)
    print str(weights)
    if nx.get_edge_attributes(g,"weight"):
        new.es['weight'] = [edge[2]['weight'] for edge in weights]
    else:
        new.es['weight'] = [1 for edge in weights]

    return new
  
G = nx.DiGraph()
G.add_edges_from([('a', 'c', {'weight': 1}),
                  ('a', 'b', {'weight': 3}),
                  ('c', 'a', {'weight': 1}),
                  ('c', 'd', {'weight': 2}),
                  ('b', 'a', {'weight': 1}),
                  ('d', 'c', {'weight': 1})])

print G.edges(data=True)

g = nx.Graph()
g.add_edges_from(G.edges_iter(), weight=0)

print g.edges(data=True)

for u, v, d in G.edges_iter(data=True):
    g[u][v]['weight'] += d['weight']

print g.edges(data=True)
ig_g = networkx_to_igraph(g)
# calculate dendrogram
dendrogram = ig_g.community_edge_betweenness()

# convert it into a flat clustering
clusters = dendrogram.as_clustering()
# get the membership vector
membership = clusters.membership
print "Clusters", str(clusters)
print clusters[0]
print "Membership", str(membership)

clusterings = ig_g.community_leading_eigenvector(clusters=2)
print str(clusterings)