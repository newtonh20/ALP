import networkx as nx
import sys
def get_graph_features(h):
    if h.is_directed():
        g = nx.Graph(h)
    else:
        g = h

    def runf(func):
        try:
            print "Running", str(func)
            val = func(g)
            print "Value:", str(val)
            sys.stdout.flush()
            return val
        except:
            return -1

    a = 1 if h.is_directed() else 0
    b = runf(nx.number_of_nodes)
    c = runf(nx.number_of_edges)
    d = runf(nx.estrada_index) if b < 1000 and not g.is_directed() else -1
    e = 1 if not g.is_directed() and runf(nx.is_chordal) else 0
    f = runf(nx.graph_clique_number)
    h = runf(nx.graph_number_of_cliques)
    i = runf(nx.transitivity)
    j = runf(nx.average_clustering)
    k = runf(nx.average_node_connectivity) if b < 1000 else -1
    l = runf(nx.edge_connectivity) if b < 1000 else -1
    m = runf(nx.node_connectivity) if b < 1000 else -1
    n = runf(nx.diameter) if b < 10000 else -1
    try:
        o = len(nx.periphery(g)) if b < 1000 else -1
    except:
        o = -1
    p = 1 if runf(nx.is_eulerian) else 0
    q = runf(nx.average_shortest_path_length) if b < 5000 else -1
    try:
        r = nx.connected_double_edge_swap(g, nx.number_of_edges(g)) if b < 5000 else -1
    except:
        r = -1

    s = 1 if runf(nx.is_tree) else 0
    t = runf(nx.density)
    return (a, b, c, d, e, f, h, i, j, k, l, m, n, o, p, q, r, s, t)
