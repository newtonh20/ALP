import networkx as nx
import random
from operator import itemgetter
from copy import deepcopy

def greedy_modularity_agglomoritive_partition(G,
                                              number_of_partitions,
                                              dendogram=False):

    """Use a greedy modularity maximization algorithm to
    combine communities to create the partition of a graph.

    This algorithm proceeds by iteratively computing the change
    in modularity from combining two communities who share an edge.
    Communities are choosen to combine based on which edge will
    increase the modularity of the partition the most. This
    algorithm can either return the partition itself or a
    dendogram of the process by which communities were combined.

    Parameters
    ----------
    G : NetworkX Graph, DiGraph or MultiGraph
      The Graph to be partitioned
    number_of_partitions : int
      Number of partitions to partition the graph
    dendogram : boolean, optional default=False
      Whether to return the dendogram of community combinations

    Returns
    -------
    C : list of sets
      Partition of the Graph
    or
    D : NetworkX Graph
      Dendogram of community combinations

    Raises
    ------
    NetworkXError
      If number of partitions is not in [1,n]

    Notes
    -----
    If a Dendogram is returned nodes have an attribute 'order' which indicates
    the order in which they were added to the dendogram. This algorithm likely
    over calculates the change in modularity, so further investigations are
    likely needed.

    Examples
    --------
    >>> G = nx.barbell_graph(3,0)
    >>> C = nx.greedy_modularity_agglomoritive_partition(G,2)
    >>> C
    [set([0,1,2]),set([3,4,5])]
    >>> D = nx.greedy_modularity_agglomoritive_partition(G,1,dendogram=True)
    >>> D.nodes(data=True)
    [(frozenset([0, 1, 2, 3, 4, 5]), {'order': 5}),
     (frozenset([0, 1, 2]), {'order': 4}),
     (frozenset([5]), {'order': 0}),
     (frozenset([3]), {'order': 0}),
     (frozenset([0]), {'order': 0}),
     (frozenset([3, 4, 5]), {'order': 2}),
     (frozenset([4]), {'order': 0}),
     (frozenset([1]), {'order': 0}),
     (frozenset([0, 1]), {'order': 3}),
     (frozenset([2]), {'order': 0}),
     (frozenset([4, 5]), {'order': 1})]

    References
    ----------
    ..[1] M.E.J. Newman 'Fast Algorithm for detecting community structure',
          Physical Review E 69 066133 2004"""
    
    if not 1 <= number_of_partitions <= G.order():
        raise nx.NetworkXError("Number of partitions must be in [0,n]")
    elif number_of_partitions == 1 and not dendogram:
        return [set(G)]
    elif number_of_partitions == n:
        if dendogram:
            node_order = 0
            D = nx.Graph()
            D.add_nodes_from(map(frozenset,[[n] for n in G]),order=node_order)
            return D
        else:
            return map(set,[[n] for n in G])
    
    C = map(set,[[n] for n in G])
    Q = nx.modularity(G,C)
    if dendogram:
        node_order = 0
        D = nx.DiGraph()
        D.add_nodes_from(map(frozenset,C),order=node_order,Q=Q)
        node_order += 1
    affil = dict(zip(G,range(len(G))))
    m = float(G.number_of_edges())
    edges = G.edges()
    while len(filter(None,C)) > number_of_partitions:
        del_q = []
        for (u,v) in edges:
            if not affil[u] == affil[v]:
                e_ij = 0.0
                for i in C[affil[u]]:
                    e_ij += len(set(G.neighbors(i)).\
                                intersection(C[affil[v]]))/(2.0*m)
                ai_aj = sum(G.degree(C[affil[u]]).values())*\
                        sum(G.degree(C[affil[v]]).values())/(4.0*(m**2.0))
                del_q.append((2.0*(e_ij - ai_aj),(u,v)))
        q,(max_u,max_v) = max(del_q,key=itemgetter(0))
        Q += q
        aff_max_v = affil[max_v]
        edges.remove((max_u,max_v))
        Cu_new = C[affil[max_u]].union(C[affil[max_v]])
        if dendogram:
            D.add_node(frozenset(Cu_new),order=node_order,Q=Q)
            D.add_edge(frozenset(C[affil[max_u]]),frozenset(Cu_new))
            D.add_edge(frozenset(C[affil[max_v]]),frozenset(Cu_new))
            node_order += 1
        C[affil[max_u]] = Cu_new
        for j in C[affil[max_v]]:
            affil[j] = affil[max_u]
        C[aff_max_v] = None
    if dendogram:
        return D
    else:
        return filter(None,C)


def greedy_modularity_dendrogram(G):
    T=nx.DiGraph()
    level = 0
    T.add_nodes_from(G, level=level)
    
    n = G.number_of_nodes()
    m = float(G.number_of_edges())
    # assign each node to a unique group 0,..,n-1
    aff = dict(zip(G, range(n)))

    edges = G.edges()
    q = 0.0
    while len(node) < n:
        for (u,v) in ( (s,t) in G.edges_iter() if aff[u] != aff[v]) ):
                e_ij = 0.0
                for i in C[affil[u]]:
                    e_ij += len(set(G.neighbors(i)).\
                                intersection(C[affil[v]]))/(2.0*m)
                ai_aj = sum(G.degree(C[affil[u]]).values())*\
                        sum(G.degree(C[affil[v]]).values())/(4.0*(m**2.0))
                del_q.append((2.0*(e_ij - ai_aj),(u,v)))
        q,(max_u,max_v) = max(del_q)
        aff_max_v = affil[max_v]
        edges.remove((max_u,max_v))
  
        D.add_node(frozenset(Cu_new),order=node_order)
        D.add_edge(frozenset(C[affil[max_u]]),frozenset(Cu_new))
        D.add_edge(frozenset(C[affil[max_v]]),frozenset(Cu_new))
        node_order += 1

        affil[j] = affil[max_u]

    return T
