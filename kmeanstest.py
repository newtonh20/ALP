import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import KMeans
import urllib, zipfile, StringIO, community

G=nx.scale_free_graph(100)
D = nx.to_numpy_matrix(G)
print str(D)
S = cosine_similarity(D) # 100 x 100
S[S <= 0.49999] = 0 # using < 0.5 is not working well here for some reason
S[S != 0] = 1
print str(S)
# S as adjacency matrix
#G = nx.from_numpy_matrix(S)


pos = nx.spring_layout(G, k=0.05)
colors = 'bgrcmykw'

nx.draw_networkx_edges(G, pos, alpha=0.5);
k = 3 # use number of clusters found by MM (8)
kmeans = KMeans(k)
kmeans.fit(D)

cluster_nodes = {}
node_clusters = {}
for i in range(k):
    list_nodes = np.where(kmeans.labels_ == i)[0].tolist()
    #list of nodes in community = cluster_nodes[community#]
    #community for node = node_clusters[node]
    cluster_nodes[i] = list_nodes
    node_clusters.update({l:i for l in list_nodes})
print cluster_nodes
print node_clusters