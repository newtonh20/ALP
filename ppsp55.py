"""
This module implements community detection.
"""
__all__ = ["partition_at_level", "modularity", "generate_dendrogram", "generate_dendogram", "induced_graph", "greedy_modularity_agglomoritive_partition"]
__author__ = """Thomas Aynaud (thomas.aynaud@lip6.fr)"""
#    Copyright (C) 2009 by
#    Thomas Aynaud <thomas.aynaud@lip6.fr>
#    All rights reserved.
#    BSD license.

__PASS_MAX = -1
__MIN = 0.0000001
import networkx as nx
import sys
import types
import array
import threading
from contextlib import closing
import MySQLdb
import math
import random as r
import gc
import pprint
from random import choice, sample
import time
from memory_profiler import profile
import os
import operator
from multiprocessing import Pool
import itertools
import collections
import heapq
from networkx.utils import generate_unique_node
from networkx import NetworkXError
from itertools import count
import numpy as np
import random
from operator import itemgetter
from copy import deepcopy
import igraph

def igraph_to_networkx(igraph_g):
    A = igraph_g.get_adjacency()
    A = np.matrix(A.data)
    return nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
def networkx_to_igraph(G):
    #new = igraph.Graph(directed=nx.is_directed(G))
    new = igraph.Graph(directed=nx.is_directed(G))
    new.add_vertices(G.nodes())
    edges = G.edges()
    new.add_edges(edges)

    weights = G.edges(data=True)
    #print str(weights)
    if nx.get_edge_attributes(g,"weight"):
        new.es['weight'] = [edge[2]['weight'] for edge in weights]
    else:
        new.es['weight'] = [1 for edge in weights]
    return new
    

def partition_at_level(dendrogram, level) :
    """Return the partition of the nodes at the given level

    A dendrogram is a tree and each level is a partition of the graph nodes.
    Level 0 is the first partition, which contains the smallest communities, and the best is len(dendrogram) - 1.
    The higher the level is, the bigger are the communities

    Parameters
    ----------
    dendrogram : list of dict
       a list of partitions, ie dictionnaries where keys of the i+1 are the values of the i.
    level : int
       the level which belongs to [0..len(dendrogram)-1]

    Returns
    -------
    partition : dictionnary
       A dictionary where keys are the nodes and the values are the set it belongs to

    Raises
    ------
    KeyError
       If the dendrogram is not well formed or the level is too high

    See Also
    --------
    best_partition which directly combines partition_at_level and generate_dendrogram to obtain the partition of highest modularity

    Examples
    --------
    >>> G=nx.erdos_renyi_graph(100, 0.01)
    >>> dendo = generate_dendrogram(G)
    >>> for level in range(len(dendo) - 1) :
    >>>     print "partition at level", level, "is", partition_at_level(dendo, level)
    """
    partition = dict()
    partition = dendrogram[0].copy()
    for index in range(1, level + 1) :
        for node, community in partition.items() :
            partition[node] = dendrogram[index][community]
    return partition


def modularity(partition, graph) :
    """Compute the modularity of a partition of a graph

    Parameters
    ----------
    partition : dict
       the partition of the nodes, i.e a dictionary where keys are their nodes and values the communities
    graph : networkx.Graph
       the networkx graph which is decomposed

    Returns
    -------
    modularity : float
       The modularity

    Raises
    ------
    KeyError
       If the partition is not a partition of all graph nodes
    ValueError
        If the graph has no link
    TypeError
        If graph is not a networkx.Graph

    References
    ----------
    .. 1. Newman, M.E.J. & Girvan, M. Finding and evaluating community structure in networks. Physical Review E 69, 26113(2004).

    Examples
    --------
    >>> G=nx.erdos_renyi_graph(100, 0.01)
    >>> part = best_partition(G)
    >>> modularity(part, G)
    """
    if type(graph) != nx.Graph :
        raise TypeError("Bad graph type, use only non directed graph")
    print "Computing modularity"
    inc = dict([])
    deg = dict([])
    links = graph.size(weight='weight')
    if links == 0 :
        raise ValueError("A graph without link has an undefined modularity")

    for node in graph :
        com = partition[node]
        deg[com] = deg.get(com, 0.) + graph.degree(node, weight = 'weight')
        for neighbor, datas in graph[node].items() :
            weight = datas.get("weight", 1)
            if partition[neighbor] == com :
                if neighbor == node :
                    inc[com] = inc.get(com, 0.) + float(weight)
                else :
                    inc[com] = inc.get(com, 0.) + float(weight) / 2.

    res = 0.
    for com in set(partition.values()) :
        res += (inc.get(com, 0.) / links) - (deg.get(com, 0.) / (2.*links))**2
    return res

def best_partition(graph, partition = None) :
    """Compute the partition of the graph nodes which maximises the modularity
    (or try..) using the Louvain heuristices

    This is the partition of highest modularity, i.e. the highest partition of the dendrogram
    generated by the Louvain algorithm.

    Parameters
    ----------
    graph : networkx.Graph
       the networkx graph which is decomposed
    partition : dict, optionnal
       the algorithm will start using this partition of the nodes. It's a dictionary where keys are their nodes and values the communities

    Returns
    -------
    partition : dictionnary
       The partition, with communities numbered from 0 to number of communities

    Raises
    ------
    NetworkXError
       If the graph is not Eulerian.

    See Also
    --------
    generate_dendrogram to obtain all the decompositions levels

    Notes
    -----
    Uses Louvain algorithm

    References
    ----------
    .. 1. Blondel, V.D. et al. Fast unfolding of communities in large networks. J. Stat. Mech 10008, 1-12(2008).

    Examples
    --------
    >>>  #Basic usage
    >>> G=nx.erdos_renyi_graph(100, 0.01)
    >>> part = best_partition(G)

    >>> #other example to display a graph with its community :
    >>> #better with karate_graph() as defined in networkx examples
    >>> #erdos renyi don't have true community structure
    >>> G = nx.erdos_renyi_graph(30, 0.05)
    >>> #first compute the best partition
    >>> partition = community.best_partition(G)
    >>>  #drawing
    >>> size = float(len(set(partition.values())))
    >>> pos = nx.spring_layout(G)
    >>> count = 0.
    >>> for com in set(partition.values()) :
    >>>     count = count + 1.
    >>>     list_nodes = [nodes for nodes in partition.keys()
    >>>                                 if partition[nodes] == com]
    >>>     nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 20,
                                    node_color = str(count / size))
    >>> nx.draw_networkx_edges(G,pos, alpha=0.5)
    >>> plt.show()
    """
    dendo = generate_dendrogram(graph, partition)
    return partition_at_level(dendo, len(dendo) - 1)

def generate_dendogram(graph, part_init = None) :
    """Deprecated, use generate_dendrogram"""
    return generate_dendrogram(graph, part_init)
    
    
def generate_dendrogram(graph, part_init = None) :
    """Find communities in the graph and return the associated dendrogram

    A dendrogram is a tree and each level is a partition of the graph nodes.  Level 0 is the first partition, which contains the smallest communities, and the best is len(dendrogram) - 1. The higher the level is, the bigger are the communities


    Parameters
    ----------
    graph : networkx.Graph
        the networkx graph which will be decomposed
    part_init : dict, optionnal
        the algorithm will start using this partition of the nodes. It's a dictionary where keys are their nodes and values the communities

    Returns
    -------
    dendrogram : list of dictionaries
        a list of partitions, ie dictionnaries where keys of the i+1 are the values of the i. and where keys of the first are the nodes of graph

    Raises
    ------
    TypeError
        If the graph is not a networkx.Graph

    See Also
    --------
    best_partition

    Notes
    -----
    Uses Louvain algorithm

    References
    ----------
    .. 1. Blondel, V.D. et al. Fast unfolding of communities in large networks. J. Stat. Mech 10008, 1-12(2008).

    Examples
    --------
    >>> G=nx.erdos_renyi_graph(100, 0.01)
    >>> dendo = generate_dendrogram(G)
    >>> for level in range(len(dendo) - 1) :
    >>>     print "partition at level", level, "is", partition_at_level(dendo, level)
    """
    if type(graph) != nx.Graph :
        raise TypeError("Bad graph type, use only non directed graph")

    #special case, when there is no link
    #the best partition is everyone in its community
    if graph.number_of_edges() == 0 :
        part = dict([])
        for node in graph.nodes() :
            part[node] = node
        return part

    current_graph = graph.copy()
    status = Status()
    status.init(current_graph, part_init)
    mod = __modularity(status)
    status_list = list()
    __one_level(current_graph, status)
    new_mod = __modularity(status)
    partition = __renumber(status.node2com)
    status_list.append(partition)
    mod = new_mod
    current_graph = induced_graph(partition, current_graph)
    status.init(current_graph)

    while True :
        __one_level(current_graph, status)
        new_mod = __modularity(status)
        if new_mod - mod < __MIN :
            break
        partition = __renumber(status.node2com)
        status_list.append(partition)
        mod = new_mod
        current_graph = induced_graph(partition, current_graph)
        status.init(current_graph)
    return status_list[:]


def induced_graph(partition, graph) :
    """Produce the graph where nodes are the communities

    there is a link of weight w between communities if the sum of the weights of the links between their elements is w

    Parameters
    ----------
    partition : dict
       a dictionary where keys are graph nodes and  values the part the node belongs to
    graph : networkx.Graph
        the initial graph

    Returns
    -------
    g : networkx.Graph
       a networkx graph where nodes are the parts

    Examples
    --------
    >>> n = 5
    >>> g = nx.complete_graph(2*n)
    >>> part = dict([])
    >>> for node in g.nodes() :
    >>>     part[node] = node % 2
    >>> ind = induced_graph(part, g)
    >>> goal = nx.Graph()
    >>> goal.add_weighted_edges_from([(0,1,n*n),(0,0,n*(n-1)/2), (1, 1, n*(n-1)/2)])
    >>> nx.is_isomorphic(int, goal)
    True
    """
    ret = nx.Graph()
    ret.add_nodes_from(partition.values())

    for node1, node2, datas in graph.edges_iter(data = True) :
        weight = datas.get("weight", 1)
        com1 = partition[node1]
        com2 = partition[node2]
        w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
        ret.add_edge(com1, com2, weight = w_prec + weight)

    return ret


def __renumber(dictionary) :
    """Renumber the values of the dictionary from 0 to n
    """
    count = 0
    ret = dictionary.copy()
    new_values = dict([])

    for key in dictionary.keys() :
        value = dictionary[key]
        new_value = new_values.get(value, -1)
        if new_value == -1 :
            new_values[value] = count
            new_value = count
            count = count + 1
        ret[key] = new_value

    return ret


def __load_binary(data) :
    """Load binary graph as used by the cpp implementation of this algorithm
    """
    data = open(data, "rb")

    reader = array.array("I")
    reader.fromfile(data, 1)
    num_nodes = reader.pop()
    reader = array.array("I")
    reader.fromfile(data, num_nodes)
    cum_deg = reader.tolist()
    num_links = reader.pop()
    reader = array.array("I")
    reader.fromfile(data, num_links)
    links = reader.tolist()
    graph = nx.Graph()
    graph.add_nodes_from(range(num_nodes))
    prec_deg = 0

    for index in range(num_nodes) :
        last_deg = cum_deg[index]
        neighbors = links[prec_deg:last_deg]
        graph.add_edges_from([(index, int(neigh)) for neigh in neighbors])
        prec_deg = last_deg

    return graph


def __one_level(graph, status) :
    """Compute one level of communities
    """
    modif = True
    nb_pass_done = 0
    cur_mod = __modularity(status)
    new_mod = cur_mod

    while modif  and nb_pass_done != __PASS_MAX :
        cur_mod = new_mod
        modif = False
        nb_pass_done += 1

        for node in graph.nodes() :
            com_node = status.node2com[node]
            degc_totw = status.gdegrees.get(node, 0.) / (status.total_weight*2.)
            neigh_communities = __neighcom(node, graph, status)
            __remove(node, com_node,
                    neigh_communities.get(com_node, 0.), status)
            best_com = com_node
            best_increase = 0
            for com, dnc in neigh_communities.items() :
                incr =  dnc  - status.degrees.get(com, 0.) * degc_totw
                if incr > best_increase :
                    best_increase = incr
                    best_com = com
            __insert(node, best_com,
                    neigh_communities.get(best_com, 0.), status)
            if best_com != com_node :
                modif = True
        new_mod = __modularity(status)
        if new_mod - cur_mod < __MIN :
            break


class Status :
    """
    To handle several data in one struct.

    Could be replaced by named tuple, but don't want to depend on python 2.6
    """
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}

    def __init__(self) :
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])

    def __str__(self) :
        return ("node2com : " + str(self.node2com) + " degrees : "
            + str(self.degrees) + " internals : " + str(self.internals)
            + " total_weight : " + str(self.total_weight))

    def copy(self) :
        """Perform a deep copy of status"""
        new_status = Status()
        new_status.node2com = self.node2com.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.total_weight = self.total_weight

    def init(self, graph, part = None) :
        """Initialize the status of a graph with every node in one community"""
        count = 0
        self.node2com = dict([])
        self.total_weight = 0
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.total_weight = graph.size(weight = 'weight')
        if part == None :
            for node in graph.nodes() :
                self.node2com[node] = count
                deg = float(graph.degree(node, weight = 'weight'))
                if deg < 0 :
                    raise ValueError("Bad graph type, use positive weights")
                self.degrees[count] = deg
                self.gdegrees[node] = deg
                self.loops[node] = float(graph.get_edge_data(node, node,
                                                 {"weight":0}).get("weight", 1))
                self.internals[count] = self.loops[node]
                count = count + 1
        else :
            for node in graph.nodes() :
                com = part[node]
                self.node2com[node] = com
                deg = float(graph.degree(node, weight = 'weight'))
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.gdegrees[node] = deg
                inc = 0.
                for neighbor, datas in graph[node].items() :
                    weight = datas.get("weight", 1)
                    if weight <= 0 :
                        raise ValueError("Bad graph type, use positive weights")
                    if part[neighbor] == com :
                        if neighbor == node :
                            inc += float(weight)
                        else :
                            inc += float(weight) / 2.
                self.internals[com] = self.internals.get(com, 0) + inc



def __neighcom(node, graph, status) :
    """
    Compute the communities in the neighborood of node in the graph given
    with the decomposition node2com
    """
    weights = {}
    for neighbor, datas in graph[node].items() :
        if neighbor != node :
            weight = datas.get("weight", 1)
            neighborcom = status.node2com[neighbor]
            weights[neighborcom] = weights.get(neighborcom, 0) + weight

    return weights


def __remove(node, com, weight, status) :
    """ Remove node from community com and modify status"""
    status.degrees[com] = ( status.degrees.get(com, 0.)
                                    - status.gdegrees.get(node, 0.) )
    status.internals[com] = float( status.internals.get(com, 0.) -
                weight - status.loops.get(node, 0.) )
    status.node2com[node] = -1


def __insert(node, com, weight, status) :
    """ Insert node into community and modify status"""
    status.node2com[node] = com
    status.degrees[com] = ( status.degrees.get(com, 0.) +
                                status.gdegrees.get(node, 0.) )
    status.internals[com] = float( status.internals.get(com, 0.) +
                        weight + status.loops.get(node, 0.) )


def __modularity(status) :
    """
    Compute the modularity of the partition of the graph faslty using status precomputed
    """
    links = float(status.total_weight)
    result = 0.
    for community in set(status.node2com.values()) :
        in_degree = status.internals.get(community, 0.)
        degree = status.degrees.get(community, 0.)
        if links > 0 :
            result = result + in_degree / links - ((degree / (2.*links))**2)
    return result

def girvan_newman (G):
    if len(G.nodes()) == 1:
        return [G.nodes()]

    def find_best_edge(G0):
        """
        Networkx implementation of edge_betweenness
        returns a dictionary. Make this into a list,
        sort it and return the edge with highest betweenness.
        """
        eb = nx.edge_betweenness_centrality(G0)
        eb_il = eb.items()
        eb_il.sort(key=lambda x: x[1], reverse=True)
        return eb_il[0][0]

    components = nx.connected_component_subgraphs(G)

    while len(components) == 1:
        G.remove_edge(*find_best_edge(G))
        components = nx.connected_component_subgraphs(G)

    result = [c.nodes() for c in components]

    for c in components:
        result.extend(girvan_newman(c))

    return result


############## DBAdapter ##########################################
__author__ = 'Newton'
#
# Database Adapter for measurements
#
test = False
DEBUG = False
debug = False
pickled = False
g = None
landmarks = []
alp_landmarks = []
num_alt_calculations = 0
num_alp_calculations = 0
num_alt_estimates = 0
cluster_nodes = None
num_alp_estimates = 0
dendo = None
tree_level = 1
num_alt_landmarks = 0
gid = 0
tid = 0
eid = 0

host = "192.168.1.111"
user = "newtonh20"
password = "Newt13hc!"
db = "alp"

test=False
def query_insert(query):
    try:
        con = MySQLdb.connect(host, user, password, db)
        with closing(con.cursor()) as cur:
            try:
                cur.execute(query)
                con.commit()
                return cur.lastrowid
            except Exception as e:
                print(str(e))
                print "Could not commit", query
                con.rollback()
                import time
                time.sleep(5)
                return query_insert(query)

        con.close()
    except:
        print "Error with query:", query
        print "Could not read/write from database."

def query(sql):
    if test:
        print "Testing,", sql
        return None
    try:
        con = MySQLdb.connect(host, user, password, db)
        with closing(con.cursor()) as cur:
            cur = con.cursor()
            cur.execute(sql)
            result = cur.fetchall()
        con.close()
    except:
        print "Error with query:", query
        print "Could not read/write from database."
        result = None

    return result


############## LANDMARKS ##########################################

#import matplotlib.pyplot as plt

#TODO: Implement this in C using http://www.cs.yale.edu/homes/aspnes/pinewiki/C(2f)Graphs.html

pos = None
def mode(d):
    c = collections.Counter(d.values())
    most_common_value = c.most_common(1)[0][0]
    for name, val in d.iteritems():
        if val == most_common_value:
            return name

def get_path_size(G, path,weight='weight'):
    return sum(G[u][v].get(weight, 1) for u, v in zip(path[:-1], path[1:]))
    '''if not nx.get_edge_attributes(G,"weight"):
        return len(path)
    tot_size = 0
    for i in range(len(path)-1):
        tot_size += G[path[i]][path[i+1]]['weight']

    return tot_size'''

#random-p
def alp_random_ls_opt(g, subgraphs, trials_to_nodes_in_graph=.1, queries_per_trial=25):
    if debug:
        print "Running ALP random landmark selection with", str(trials_to_nodes_in_graph), \
            "trial ratio and", str(queries_per_trial), "queries per trial."
    #perform the landmark selection within each subgraph
    if g.is_directed():
        #print subgraphs.values()
        #perform on the connected subgraph (what can be reached by dijkstra's)
        #landmarks = [random_ls_opt(g, 1, trials_to_nodes_in_graph, queries_per_trial, candidate_nodes=sub) if len(sub) > 1 else sub[0] for sub in subgraphs.values()]
        #landmarks = [random_ls_opt(g.subgraph(single_source_dijkstra_path_length(g, sub[0], target_cutoffs=sub, only_targets=False).keys()), 1, trials_to_nodes_in_graph, queries_per_trial, candidate_nodes=sub) if len(sub) > 1 else sub[0] for sub in subgraphs.values()]
        landmarks = [random_ls_opt(list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0], 1, trials_to_nodes_in_graph, queries_per_trial) for sub in subgraphs.values()]
    else:
        landmarks = [random_ls_opt(list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0], 1, trials_to_nodes_in_graph, queries_per_trial) for sub in subgraphs.values()]
    return [l[0] for l in landmarks if l]

#performs random landmark selection for ALP without trials
def alp_random_ls_non_opt(g,subgraphs):
    #perform the landmark selection within each subgraph
    if g.is_directed():
        return [sample(sub, 1)[0] for sub in subgraphs.values()]
        #landmarks = [random_ls_non_opt(list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0], 1) for sub in subgraphs.values()]
    else:
        landmarks = [random_ls_non_opt(list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0], 1) for sub in subgraphs.values()]
    return [l[0] for l in landmarks if l]



#performs random landmark selection without trials
def random_ls_non_opt(g, k):
    #[choice(g.nodes()) for i in range(k)]
    return sample(g.nodes(), k)

# ALT random landmark selection
# Based on the number of vertices in the graph, k vertices are chosen at random to serve as landmarks.
# A series of sample queries are run with each landmark to determine the best set. This is a brute
# force method of performing landmark selection for ALT. However, in terms of lower bounds, random
# landmarks demonstrate better performance than any of the following methods of landmark selection
#candidate nodes parameter implies we have a set we would like to choose from
def random_ls_opt(g, k, trials_to_nodes_in_graph=.1, queries_per_trial=25, candidate_nodes=[]):
    if k < 1:
        raise ValueError("Cannot specify " + str(k) + " as number of landmarks.")

    temp_ls = []
    n_len = nx.number_of_nodes(g)

    if n_len < 3 and k == 1:
        return [g.nodes()[0]] if not candidate_nodes else [candidate_nodes[0]]

    def temp_ALT(v, end):
        return max([abs(g.node[v][str(l)] - g.node[end][str(l)]) for l in temp_ls])

    min_ls = []
    num_trials = int(n_len * trials_to_nodes_in_graph)
    #try at least two trials
    num_trials = 2 if num_trials <= 2 else num_trials
    if candidate_nodes and num_trials > len(candidate_nodes):
        num_trials = len(candidate_nodes)
    '''if debug:
        print "Running", str(num_trials), "trials for random landmark selection."'''
    #grab the one that has the smallest avg search space
    min_search_space = n_len

    #establish a set of test nodes such that we can vet each landmark set
    test_nodes = []
    for j in range(queries_per_trial):
        samps = sample(g.nodes(), 2)
        test_nodes.append((samps[0], samps[1]))
    #igraph_g = networkx_to_igraph(g)
    #for each landmark configuration
    for i in range(num_trials):
        if debug:
            print "Running random landmark selection, trial", str(i)
        #choose a set of landmarks
        temp_ls = sample(g.nodes(), k) if not candidate_nodes else sample(candidate_nodes, k)
        temp_ls = list(set(temp_ls))
        #[choice(g.nodes()) for i in range(k)]
        for l in temp_ls:
            nx.set_node_attributes(g, str(l), nx.single_source_dijkstra_path_length(g, l))
        
            

        #issue a series of arbitrary shortest path queries using the temp_ALT algorithm as the heuristic
        #compare
        total_search_space = 0
        for j in range(queries_per_trial):
            avg_search_space = n_len
            try:
                for (s,t) in test_nodes:
                    path, astar_size = astar_path(g,s,t,temp_ALT, search_space_size=True)
                    total_search_space += astar_size
                avg_search_space = total_search_space/queries_per_trial
            except nx.NetworkXNoPath:
                pass
            if avg_search_space < min_search_space:
                min_search_space = avg_search_space
                min_ls = temp_ls

        #clean up graph labeling
        for l in temp_ls:
            for n in g.nodes():
                del g.node[n][str(l)]
    #print "Subgraph size:", str(n_len), ", Landmark chosen (Vertex,Search Space): ", str((min_ls, min_search_space))
    return min_ls

def alp_pagerank(g, subgraphs, func=max, custom_alpha=0.9, custom_max_iter=10000):
    import operator
    landmarks = []
    for sub in subgraphs.values():
        if g.is_directed():
            g_s = list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0]
        else:
            g_s = list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0]

        n_len = nx.number_of_nodes(g_s)

        if n_len < 3:
            landmarks.append(g_s.nodes()[0])
        else:
            try:
                pr = nx.pagerank(g_s, alpha=custom_alpha, max_iter=custom_max_iter)
            except:
                print g_s.nodes()
                exit(1)

            if func == mode:
                landmarks.append(mode(pr))
            else:
                landmarks.append(func(pr.iteritems(), key=operator.itemgetter(1))[0])

    return landmarks

#remove all but the keys
def remove_all_but(keys, the_dict):
    return {k:the_dict[k] for k in keys}

#eccentricity-based farthest
#ALP farthest concentrates on finding landmarks in each cluster that have the highest eccentricity in the overall graph
#Identifying a new landmark by finding the furthest landmark from the rest in the set would result, for each clus
def alp_farthest_ecc(g, subgraphs, opt = False, trials=10, dijkstra_sample_nodes=100):
    if not g:
        raise ValueError("Cannot specify graph as null.")
    elif not subgraphs:
        raise ValueError("Cannot specify subgraphs as null.")

    ls = []
    #dictionary of eccentricities
    e_dict = {}
    reverse_graph = g.reverse(True)
    d_samples = sample(g.nodes(), dijkstra_sample_nodes)
    #get the distance from other nodes to the landmark
    
    def run_dij(samp):
        e_dict[s] = max(single_source_dijkstra_path_length(reverse_graph, samp, target_cutoffs=d_samples[:], only_targets=True).values())
    
    for g_p in subgraphs.values():
        samples = sample(g_p, trials) if trials <= len(g_p) else g_p
        for s in samples:
            threading.Thread(target=run_dij, args=(s,)).start()
        
        for t in threading.enumerate():
            main_thread = threading.currentThread()
            if t is main_thread:
                continue
            t.join()
            
    #e = nx.eccentricity(g)
    i = 0
    for g_p in subgraphs.values():
        max_nodes = []
        max_ecc = 0
        #first, get the nodes with max eccentricity in the cluster
        for v in [k for k in e_dict.keys() if k in g_p]:
            if e_dict[v] > max_ecc:
                max_nodes = [v]
                max_ecc = e_dict[v]
            elif e_dict[v] == max_ecc:
                max_nodes.append(v)
        ls.append(choice(max_nodes))

    return ls

#Gets the index of maximum element in a list. If a conflict occurs, the index of the last largest is returned
def maxl(l): 
    return l.index(reduce(lambda x,y: max(x,y), l))

def lsum(first,second):
    return [x + y for x, y in zip(first, second)]

#cluster-based farthest
#ALP farthest concentrates on finding landmarks in each cluster that are farthest from all other landmarks in the graph
#Identifying a new landmark by finding the furthest landmark from the rest in the set would result, for each clus
def alp_farthest(g, subgraphs, opt = False):
    if not g:
        raise ValueError("Cannot specify graph as null.")
    elif not subgraphs:
        raise ValueError("Cannot specify subgraphs as null.")

    ls = []

    #e = nx.eccentricity(g)
    i = 0
    igraph_g = networkx_to_igraph(g)
    for g_p in subgraphs.values():
        i += 1
        #print str(i)
        #print g_p
        
        if not ls:
            ls.append(choice(g_p))
        else:
            paths = igraph_g.shortest_paths(ls,g_p[:],'weight',igraph.OUT)
            psum = paths[0]
            for i in range(1,len(paths)):
                psum = lsum(psum, paths[i])
                
            ls.append(g_p[maxl(psum)])
            '''ls_dists = []
            for l in ls:
                if nx.get_edge_attributes(g,"weight") and not opt:
                    ls_dists.append(remove_all_but(g_p, single_source_dijkstra_path_length(g, l, target_cutoffs=g_p[:])))
                else:
                    ls_dists.append(remove_all_but(g_p, single_source_shortest_path_length(g, l, target_cutoffs=g_p[:])))


            #print "Ran dijkstra:", ls_dists
            #find the node furthest from current landmarks
            max_d = 0

            max_node = g_p[0]
            #identify node with max eccentricity
            for node in g_p:
                total = 0
                for dists in ls_dists:
                    total += dists[node]

                avg_dist = total/len(ls_dists)
                if avg_dist > max_d:
                    max_d = avg_dist
                    max_node = node

            #print "Max dist=", str(max_d)
            #print "Max node=", str(max_node)
            #add to set
            ls.append(max_node)
            #print ls'''
    return ls

def farthest(g, k, farthest_cutoff=0):
    if k < 1:
        raise ValueError("Cannot specify " + str(k) + " as number of landmarks.")
    elif farthest_cutoff < 0:
        raise ValueError("Cannot specify " + str(farthest_cutoff) + " as cutoff.")

    #landmark set
    ls = []
    #we must store the tree lengths to compute the best upper bounds later
    ls_dists = []

    #identify a start vertex v in V
    ls.append(choice(g.nodes()))

    #add v to the set of landmarks
    while(len(ls) != k):
        v_p = ls[-1]

        #identify a vertex v'' farthest from v' (which is the most recent in the set)
        cf = nx.number_of_nodes(g) if farthest_cutoff == 0 else farthest_cutoff
        depths = single_source_dijkstra_path_length(g, v_p, target_cutoffs=cf)

        #add the lengths dictionary to the landmark dists set
        ls_dists.append(depths)

        #get the distance that maximizes the upper bounds
        max_furth = v_p
        max_furth_dist =0

        #now go through v' nodes for depths and identify max
        for node in depths.keys():
            upper = 0
            #uppder_bound = d(v1,v2) + d(v2,v3) + ... + d(v_n-1, v_n)
            for dist in ls_dists:
                upper += dist[node]

            #now see if this is the highest upper bound
            if upper > max_furth_dist:
                max_furth = node
                max_furth_dist = upper

        #add v_p to the set of landmarks
        ls.append(max_furth)

    return ls


''' Implementation of ALT farthest landmark selection that leverages the upper
    bound of the generalized polygon inequality to
    quickly determine the farthest landmark from the remainder of the set'''
def farthest_d(g, k, farthest_cutoff=0):
    if k < 1:
        raise ValueError("Cannot specify " + str(k) + " as number of landmarks.")
    elif farthest_cutoff < 0:
        raise ValueError("Cannot specify " + str(farthest_cutoff) + " as cutoff.")

    #landmark set
    ls = []
    #we must store the tree lengths to compute the best upper bounds later
    ls_dists = []

    #identify a start vertex v in V
    ls.append(choice(g.nodes()))

    #add v to the set of landmarks
    while(len(ls) < k):
        v_p = ls[-1]

        #identify a vertex v'' farthest from v' (which is the most recent in the set)
        cf = nx.number_of_nodes(g) if farthest_cutoff == 0 else farthest_cutoff
        depths = single_source_shortest_path_length(g, v_p, cutoff=cf, target_cutoffs=g.nodes())

        #add the lengths dictionary to the landmark dists set
        ls_dists.append(depths)

        #get the distance that maximizes the upper bounds
        max_furth = v_p
        max_furth_dist =0

        #now go through v' nodes for depths and identify max
        for node in depths.keys():
            upper = 0
            #uppder_bound = d(v1,v2) + d(v2,v3) + ... + d(v_n-1, v_n)
            for dist in ls_dists:
                upper += dist[node]

            #now see if this is the highest upper bound
            if upper > max_furth_dist:
                max_furth = node
                max_furth_dist = upper

        #add v_p to the set of landmarks
        ls.append(max_furth)

    return ls

def planar(g, k):
    global pos
    if not pos:
        pos=nx.spectral_layout(g)
    #nx.draw(g, pos=pos)
    #plt.savefig(g.graph['name']+'.png')
    pos_planar = {key:[value[0]-.5,value[1]-.5] for (key,value) in pos.iteritems()}
    ls = []
    sector = 0

    # we label the sectors as follows
    #            |
    #       0    |     1
    #       -,+  |     +,+
    #   --------------------
    #       -,-  |    +,-
    #       2    |     3
    #
    '''zeros = {key:value for (key, value) in pos_planar.iteritems() if value[0] < 0 and value[1] > 0}
    ones =  {key:value for (key, value) in pos_planar.iteritems() if value[0] >= 0 and value[1] > 0}
    twos =  {key:value for (key, value) in pos_planar.iteritems() if value[0] < 0 and value[1] <= 0}
    threes =  {key:value for (key, value) in pos_planar.iteritems() if value[0] > 0 and value[1] < 0}'''

    while(len(ls) < k):
        if sector == 0:
            zeros = {key:value for (key, value) in pos_planar.iteritems() if value[0] < 0 and value[1] > 0}
            if zeros:
                #find the farthest point (most negative/positive point) in the set
                lcand = zeros.keys()[0]
                min_x = zeros[lcand][0]
                max_y = zeros[lcand][1]
                for key, value in zeros.iteritems():
                    if DEBUG:
                        print value
                    if value[0] < min_x and value[1] > max_y:
                        lcand = key
                        min_x = value[0]
                        max_y = value[1]

                ls.append(lcand)
                del pos_planar[lcand]
            sector += 1
        elif sector == 1:
            ones =  {key:value for (key, value) in pos_planar.iteritems() if value[0] >= 0 and value[1] > 0}
            if ones:
                #find the farthest point (most positive/positive point) in the set
                lcand = ones.keys()[0]
                max_x = ones[lcand][0]
                max_y = ones[lcand][1]
                for key, value in ones.iteritems():
                    if DEBUG:
                        print value
                    if value[0] >= max_x and value[1] > max_y:
                        lcand = key
                        max_x = value[0]
                        max_y = value[1]

                ls.append(lcand)
                del pos_planar[lcand]
            sector += 1
        elif sector == 2:
            twos =  {key:value for (key, value) in pos_planar.iteritems() if value[0] < 0 and value[1] <= 0}
            if twos:
                #find the farthest point (most negative/negative point) in the set
                lcand = twos.keys()[0]
                min_x = twos[lcand][0]
                min_y = twos[lcand][1]
                for key, value in twos.iteritems():
                    if DEBUG:
                        print value
                    if value[0] < min_x and value[1] <= min_y:
                        lcand = key
                        min_x = value[0]
                        min_y = value[1]

                ls.append(lcand)
                del pos_planar[lcand]
            sector += 1
        elif sector == 3:
            threes =  {key:value for (key, value) in pos_planar.iteritems() if value[0] > 0 and value[1] < 0}
            if threes:
                #find the farthest point (most positive/negative point) in the set
                lcand = threes.keys()[0]
                max_x = threes[lcand][0]
                min_y = threes[lcand][1]
                for key, value in threes.iteritems():
                    if DEBUG:
                        print value
                    if value[0] > max_x and value[1] < min_y:
                        lcand = key
                        max_x = value[0]
                        min_y = value[1]

                ls.append(lcand)
                del pos_planar[lcand]
            sector = 0
    return ls


def anti_planar(g, subgraphs):
    if not g:
        raise ValueError("Cannot specify graph as null.")
    elif not subgraphs:
        raise ValueError("Cannot specify subgraphs as null.")

    ls = []
    import operator
    for g_p in subgraphs.values():
        if g.is_directed():
            subg = list(nx.strongly_connected_component_subgraphs(g.subgraph(g_p)))[0]
        else:
            subg = list(nx.connected_component_subgraphs(g.subgraph(g_p), False))[0]
        e = nx.eccentricity(subg)
        ls.append(min(e, key=e.get))
        
    return ls   

def alp_planar(g, subgraphs):
    if not g:
        raise ValueError("Cannot specify graph as null.")
    elif not subgraphs:
        raise ValueError("Cannot specify subgraphs as null.")

    ls = []
    import operator
    for g_p in subgraphs.values():
        if g.is_directed():
            subg = list(nx.strongly_connected_component_subgraphs(g.subgraph(g_p)))[0]
        else:
            subg = list(nx.connected_component_subgraphs(g.subgraph(g_p), False))[0]

        ls.append(nx.periphery(subg)[0])
        
    return ls

def alp_betweenness(g, subgraphs, func=max):
    import operator
    landmarks = []
    for sub in subgraphs.values():
        if g.is_directed():
            g_s = list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0]
        else:
            g_s = list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0]

        n_len = nx.number_of_nodes(g_s)

        if n_len < 3:
            landmarks.append(g_s.nodes()[0])
        else:
            try:
                bc = nx.betweenness_centrality(g_s, normalized=True)
            except:
                print "Betweenness Centrality returned nothing."
                return []


            if func == mode:
                landmarks.append(mode(bc))
            else:
                landmarks.append(func(bc.iteritems(), key=operator.itemgetter(1))[0])

    return landmarks

def alp_closeness(g, subgraphs, func=max):
    import operator
    landmarks = []
    for sub in subgraphs.values():
        if g.is_directed():
            g_s = list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0]
        else:
            g_s = list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0]

        n_len = nx.number_of_nodes(g_s)

        if n_len < 3:
            landmarks.append(g_s.nodes()[0])
        else:
            try:
                bc = nx.closeness_centrality(g_s, normalized=True)
            except:
                print "Closeness Centrality returned nothing."
                return []


            if func == mode:
                landmarks.append(mode(bc))
            else:
                landmarks.append(func(bc.iteritems(), key=operator.itemgetter(1))[0])

    return landmarks

def alp_katz(g, subgraphs, func=max):
    import operator
    landmarks = []
    for sub in subgraphs.values():
        if g.is_directed():
            g_s = list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0]
        else:
            g_s = list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0]

        n_len = nx.number_of_nodes(g_s)

        if n_len < 10:
            landmarks.append(g_s.nodes()[0])
        else:
            try:
                bc = nx.katz_centrality(g_s, max_iter=100000, normalized=True)
            except:
                try:
                    print "Running Katz centrality again for tolerance = ", str(.001)
                    bc = nx.katz_centrality(g_s, max_iter=100000, normalized=True, tol=.001)
                except:
                    print "Katz Centrality returned nothing."
                    return []


            if func == mode:
                landmarks.append(mode(bc))
            else:
                landmarks.append(func(bc.iteritems(), key=operator.itemgetter(1))[0])

    return landmarks

def alp_load(g, subgraphs, func=max):
    import operator
    landmarks = []
    for sub in subgraphs.values():
        if g.is_directed():
            g_s = list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0]
        else:
            g_s = list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0]

        n_len = nx.number_of_nodes(g_s)

        if n_len < 3:
            landmarks.append(g_s.nodes()[0])
        else:
            try:
                bc = nx.load_centrality(g_s, normalized=True)
            except:
                print "Load Centrality returned nothing."
                return []


            if func == mode:
                landmarks.append(mode(bc))
            else:
                landmarks.append(func(bc.iteritems(), key=operator.itemgetter(1))[0])

    return landmarks

################NETWORKK PATHPLANNING################################
__author__ = 'Acer'



#TODO: Implement this in C using http://www.cs.yale.edu/homes/aspnes/pinewiki/C(2f)Graphs.html
# -*- coding: utf-8 -*-
"""
Shortest path algorithms for weighed graphs.
"""
#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.



def single_source_shortest_path_length(G,source,cutoff=None,target_cutoffs=[]):
    """Compute the shortest path lengths from source to all reachable nodes.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    cutoff : integer, optional
        Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    lengths : dictionary
        Dictionary of shortest path lengths keyed by target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length=nx.single_source_shortest_path_length(G,0)
    >>> length[4]
    4
    >>> print(length)
    {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}

    See Also
    --------
    shortest_path_length
    """
    tc = []
    if target_cutoffs:
        tc = set(target_cutoffs[:])

    seen=set()                  # level (number of hops) when seen in BFS
    seen_dict = {}
    level=0                  # the current level
    nextlevel={source:1}  # dict of nodes to check at next level
    while nextlevel:
        thislevel=nextlevel  # advance to next level
        nextlevel={}         # and start a new list (fringe)

        for v in thislevel:
            if v not in seen:
                seen.add(v)
                nextlevel.update(G[v]) # add neighbors of v
                if v in tc:
                    seen_dict[v] = level
                    tc.remove(v)
                    if not tc:
                        return seen_dict

        level=level+1
    return seen_dict  # return all path lengths as dictionary

def dijkstra_path(G, source, target, weight='weight'):
    """Returns the shortest path from source to target in a weighted graph G.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node

    target : node
       Ending node

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    Returns
    -------
    path : list
       List of nodes in a shortest path.

    Raises
    ------
    NetworkXNoPath
       If no path exists between source and target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> print(nx.dijkstra_path(G,0,4))
    [0, 1, 2, 3, 4]

    Notes
    ------
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    See Also
    --------
    bidirectional_dijkstra()
    """
    (length,path)=single_source_dijkstra(G, source, target=target,
                                         weight=weight)
    try:
        #print "Length to", str(target) + ":", str(path[target])
        return path[target]
    except KeyError:
        raise nx.NetworkXNoPath("node %s not reachable from %s"%(source,target))


def dijkstra_path_length(G, source, target, weight='weight'):
    """Returns the shortest path length from source to target
    in a weighted graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       starting node for path

    target : node label
       ending node for path

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    Returns
    -------
    length : number
        Shortest path length.

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> print(nx.dijkstra_path_length(G,0,4))
    4

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    See Also
    --------
    bidirectional_dijkstra()
    """
    length=single_source_dijkstra_path_length(G, source, weight=weight)
    try:
        return length[target]
    except KeyError:
        raise nx.NetworkXNoPath("node %s not reachable from %s"%(source,target))


def single_source_dijkstra_path(G,source, cutoff=None, weight='weight'):
    """Compute shortest path between source and all other reachable
    nodes for a weighted graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path.

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    paths : dictionary
       Dictionary of shortest path lengths keyed by target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> path=nx.single_source_dijkstra_path(G,0)
    >>> path[4]
    [0, 1, 2, 3, 4]

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    See Also
    --------
    single_source_dijkstra()

    """
    (length,path)=single_source_dijkstra(G,source, weight = weight)
    return path



def single_source_dijkstra_path_length(G, source, weight= 'weight', target_cutoffs=[], only_targets=False):
    #print "Target cutoffs:", target_cutoffs
    """Compute the shortest path length between source and all other
    reachable nodes for a weighted graph.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight.

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    target_cutoffs: list, optional
        if specified, stop the search and return the discovered lengths

    Returns
    -------
    length : dictionary
       Dictionary of shortest lengths keyed by target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length=nx.single_source_dijkstra_path_length(G,0)
    >>> length[4]
    4
    >>> print(length)
    {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    See Also
    --------
    single_source_dijkstra()

    """
    push = heapq.heappush
    pop = heapq.heappop
    dist = {}  # dictionary of final distances
    final_dist ={}
    target_cutoffs = set(target_cutoffs)

    if source in target_cutoffs:
        target_cutoffs.remove(source)
    seen = {source:0}
    c = count()
    fringe=[] # use heapq with (distance,label) tuples
    push(fringe, (0, next(c), source))
    while fringe:
        (d, _, v) = pop(fringe)

        if v in dist:
            continue # already searched this node.

        dist[v] = d
        if v in target_cutoffs:
            target_cutoffs.remove(v)
            final_dist[v] = d
            if not target_cutoffs:
                #print dist
                return final_dist if only_targets else dist

        #for ignore,w,edgedata in G.edges_iter(v,data=True):
        #is about 30% slower than the following
        edata=iter(G[v].items())

        for w,edgedata in iter(G[v].items()):
            vw_dist = dist[v] + edgedata.get(weight,1)

            if w not in seen or vw_dist < seen[w]:
                seen[w] = vw_dist
                push(fringe, (vw_dist, next(c), w))
                #push(fringe,(vw_dist,w))
    if target_cutoffs:
        raise ValueError("There are still target cutoffs:", str(target_cutoffs))
    return dist


def single_source_dijkstra(G,source,target=None,cutoff=None,weight='weight'):


    """Compute shortest paths and lengths in a weighted graph G.

    Uses Dijkstra's algorithm for shortest paths.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    target : node label, optional
       Ending node for path

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    distance,path : dictionaries
       Returns a tuple of two dictionaries keyed by node.
       The first dictionary stores distance from the source.
       The second stores the path from the source to that node.


    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length,path=nx.single_source_dijkstra(G,0)
    >>> print(length[4])
    4
    >>> print(length)
    {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}
    >>> path[4]
    [0, 1, 2, 3, 4]

    Notes
    ---------
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    Based on the Python cookbook recipe (119466) at
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/119466

    This algorithm is not guaranteed to work if edge weights
    are negative or are floating point numbers
    (overflows and roundoff errors can cause problems).

    See Also
    --------
    single_source_dijkstra_path()
    single_source_dijkstra_path_length()
    """
    if source==target:
        return (0, [source])
    dist = {}  # dictionary of final distances
    paths = {source:[source]}  # dictionary of paths
    seen = {source:0}
    fringe=[] # use heapq with (distance,label) tuples
    heapq.heappush(fringe,(0,source))
    while fringe:
        (d,v)=heapq.heappop(fringe)
        if v in dist:
            continue # already searched this node.
        dist[v] = d
        if v == target:
            break
        #for ignore,w,edgedata in G.edges_iter(v,data=True):
        #is about 30% slower than the following
        if G.is_multigraph():
            edata=[]
            for w,keydata in G[v].items():
                minweight=min((dd.get(weight,1)
                               for k,dd in keydata.items()))
                edata.append((w,{weight:minweight}))
        else:
            edata=iter(G[v].items())

        for w,edgedata in edata:
            vw_dist = dist[v] + edgedata.get(weight,1)

            if w in dist:
                if vw_dist < dist[w]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif w not in seen or vw_dist < seen[w]:
                seen[w] = vw_dist
                heapq.heappush(fringe,(vw_dist,w))
                paths[w] = paths[v]+[w]
    return (dist,paths)

# -*- coding: utf-8 -*-
"""Shortest paths and path lengths using A* ("A star") algorithm.
"""

#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
def astar_path(G, source, target, heuristic=None, weight='weight', search_space_nodes=False, search_space_size=False):
    """Return a list of nodes in a shortest path between source and target
    using the A* ("A-star") algorithm.

    There may be more than one shortest path.  This returns only one.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path

    heuristic : function
       A function to evaluate the estimate of the distance
       from the a node to the target.  The function takes
       two nodes arguments and must return a number.

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight.

    search_space_nodes: boolean, optional (default=False)
        Keep track of visited nodes and return them when returning the path

    search_space_size: boolean, optional (default=False)
        Count the number of visited nodes and return them when returning the path

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> print(nx.astar_path(G,0,4))
    [0, 1, 2, 3, 4]
    >>> G=nx.grid_graph(dim=[3,3])  # nodes are two-tuples (x,y)
    >>> def dist(a, b):
    ...    (x1, y1) = a
    ...    (x2, y2) = b
    ...    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    >>> print(nx.astar_path(G,(0,0),(2,2),dist))
    [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)]


    See Also
    --------
    shortest_path, dijkstra_path

    """
    if G.is_multigraph():
        raise NetworkXError("astar_path() not implemented for Multi(Di)Graphs")

    if not heuristic:
        # The default heuristic is h=0 - same as Dijkstra's algorithm
        def heuristic(u, v):
            return 0

    push = heapq.heappush
    pop = heapq.heappop

    # The queue stores priority, node, cost to reach, and parent.
    # Uses Python heapq to keep in priority order.
    # Add a counter to the queue to prevent the underlying heap from
    # attempting to compare the nodes themselves. The hash breaks ties in the
    # priority and is guarenteed unique for all nodes in the graph.
    c = heapq.count()
    queue = [(0, next(c), source, 0, None)]

    # Maps enqueued nodes to distance of discovered paths and the
    # computed heuristics to target. We avoid computing the heuristics
    # more than once and inserting the node into the queue too many times.
    enqueued = {}
    # Maps explored nodes to parent closest to the source.
    explored = {}
    i = 0
    while queue:
        i = i + 1
        # Pop the smallest item from queue.
        _, __, curnode, dist, parent = pop(queue)

        if curnode == target:
            path = [curnode]
            node = parent
            while node:
                path.append(node)
                node = explored[node]
            path.reverse()
            if search_space_nodes:
                return path, explored.keys()
            elif search_space_size:
                return path, len(explored.keys())
            elif h:
                return path, i #len(explored.keys())

            return path

        if curnode in explored:
            continue

        explored[curnode] = parent

        for neighbor, w in G[curnode].items():
            if neighbor in explored:
                continue
            ncost = dist + w.get(weight, 1)
            if neighbor in enqueued:
                qcost, h = enqueued[neighbor]
                # if qcost < ncost, a longer path to neighbor remains
                # enqueued. Removing it would need to filter the whole
                # queue, it's better just to leave it there and ignore
                # it when we visit the node a second time.
                if qcost <= ncost:
                    continue
            else:
                h = heuristic(neighbor, target)
            enqueued[neighbor] = ncost, h
            push(queue, (ncost + h, next(c), neighbor, ncost, curnode))

    raise nx.NetworkXNoPath("Node %s not reachable from %s" % (source, target))

#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
def astar_path_pathmax(G, source, target, heuristic=None, weight='weight', search_space_nodes=False, search_space_size=False):
    """Return a list of nodes in a shortest path between source and target
    using the A* ("A-star") algorithm.

    There may be more than one shortest path.  This returns only one.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path

    heuristic : function
       A function to evaluate the estimate of the distance
       from the a node to the target.  The function takes
       two nodes arguments and must return a number.

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight.

    search_space_nodes: boolean, optional (default=False)
        Keep track of visited nodes and return them when returning the path

    search_space_size: boolean, optional (default=False)
        Count the number of visited nodes and return them when returning the path

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> print(nx.astar_path(G,0,4))
    [0, 1, 2, 3, 4]
    >>> G=nx.grid_graph(dim=[3,3])  # nodes are two-tuples (x,y)
    >>> def dist(a, b):
    ...    (x1, y1) = a
    ...    (x2, y2) = b
    ...    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    >>> print(nx.astar_path(G,(0,0),(2,2),dist))
    [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)]


    See Also
    --------
    shortest_path, dijkstra_path

    """
    if G.is_multigraph():
        raise NetworkXError("astar_path() not implemented for Multi(Di)Graphs")

    if not heuristic:
        # The default heuristic is h=0 - same as Dijkstra's algorithm
        def heuristic(u, v):
            return 0

    push = heapq.heappush
    pop = heapq.heappop

    # The queue stores priority, node, cost to reach, and parent.
    # Uses Python heapq to keep in priority order.
    # Add a counter to the queue to prevent the underlying heap from
    # attempting to compare the nodes themselves. The hash breaks ties in the
    # priority and is guarenteed unique for all nodes in the graph.
    c = heapq.count()
    queue = [(0, next(c), source, 0, None)]

    # Maps enqueued nodes to distance of discovered paths and the
    # computed heuristics to target. We avoid computing the heuristics
    # more than once and inserting the node into the queue too many times.
    enqueued = {}
    # Maps explored nodes to parent closest to the source.
    explored = {}
    i = 0
    # For consistency, store the estimates
    estimates = {source:heuristic(source,target)}
    first = False
    while queue:
        i = i + 1
        # Pop the smallest item from queue.
        _, __, curnode, dist, parent = pop(queue)

        if curnode == target:
            path = [curnode]
            node = parent
            while node:
                path.append(node)
                node = explored[node]
            path.reverse()
            if search_space_nodes:
                return path, explored.keys()
            elif search_space_size:
                return path, len(explored.keys())+1
            elif h:
                return path, i #len(explored.keys())

            return path

        if curnode in explored:
            continue

        explored[curnode] = parent

        for neighbor, w in G[curnode].items():
            if neighbor in explored:
                continue
            ncost = dist + w.get(weight, 1)
            if neighbor in enqueued:
                qcost, h = enqueued[neighbor]
                # if qcost < ncost, a longer path to neighbor remains
                # enqueued. Removing it would need to filter the whole
                # queue, it's better just to leave it there and ignore
                # it when we visit the node a second time.
                if qcost <= ncost:
                    continue
            else:
                h = max(heuristic(neighbor, target), estimates[curnode]-w.get(weight,1))
                estimates[neighbor] = h
            enqueued[neighbor] = ncost, h
            push(queue, (ncost + h, next(c), neighbor, ncost, curnode))

    raise nx.NetworkXNoPath("Node %s not reachable from %s" % (source, target))

def astar_path_length_pathmax(G, source, target, heuristic=None, weight='weight'):
    """Return the length of the shortest path between source and target using
    the A* ("A-star") algorithm.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path

    heuristic : function
       A function to evaluate the estimate of the distance
       from the a node to the target.  The function takes
       two nodes arguments and must return a number.

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    See Also
    --------
    astar_path

    """
    path = astar_path_pathmax(G, source, target, heuristic, weight)
    return sum(G[u][v].get(weight, 1) for u, v in zip(path[:-1], path[1:]))

    

def astar_path_length(G, source, target, heuristic=None, weight='weight'):
    """Return the length of the shortest path between source and target using
    the A* ("A-star") algorithm.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node for path

    target : node
       Ending node for path

    heuristic : function
       A function to evaluate the estimate of the distance
       from the a node to the target.  The function takes
       two nodes arguments and must return a number.

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    See Also
    --------
    astar_path

    """
    path = astar_path(G, source, target, heuristic, weight)
    return sum(G[u][v].get(weight, 1) for u, v in zip(path[:-1], path[1:]))


def dijkstra_predecessor_and_distance(G,source, cutoff=None, weight='weight'):
    """Compute shortest path length and predecessors on shortest paths
    in weighted graphs.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    pred,distance : dictionaries
       Returns two dictionaries representing a list of predecessors
       of a node and the distance to each node.

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    The list of predecessors contains more than one element only when
    there are more than one shortest paths to the key node.
    """
    push=heapq.heappush
    pop=heapq.heappop
    dist = {}  # dictionary of final distances
    pred = {source:[]}  # dictionary of predecessors
    seen = {source:0}
    fringe=[] # use heapq with (distance,label) tuples
    push(fringe,(0,source))
    while fringe:
        (d,v)=pop(fringe)
        if v in dist: continue # already searched this node.
        dist[v] = d
        if G.is_multigraph():
            edata=[]
            for w,keydata in G[v].items():
                minweight=min((dd.get(weight,1)
                               for k,dd in keydata.items()))
                edata.append((w,{weight:minweight}))
        else:
            edata=iter(G[v].items())
        for w,edgedata in edata:
            vw_dist = dist[v] + edgedata.get(weight,1)
            if cutoff:
                if vw_dist>cutoff:
                    continue
            if w in dist:
                if vw_dist < dist[w]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif w not in seen or vw_dist < seen[w]:
                seen[w] = vw_dist
                push(fringe,(vw_dist,w))
                pred[w] = [v]
            elif vw_dist==seen[w]:
                pred[w].append(v)
    return (pred,dist)

##########ppsp45#####################################
__author__ = 'ncampbel'
'''Optimal ALT vs ALP Real Graph Trials Over Different Landmarks'''

'''Dijkstra(G =(V,E))'''
def enum(**enums):
    return type('Enum', (object,), enums)

heuristics = enum(DIJKSTRA=1, ALT=2, ALP=3, PCD=4, ALT_OPT=5, ALP_OPT=6,ALT_ND_PATHMAX=11,ALP_ND_PATHMAX=12, ALT_ND_PATHMAX_KMEANS=13,ALP_ND_PATHMAX_KMEANS=14, ALT_ND_PATHMAX_EDGE_BETWEENNESS=15, ALP_ND_PATHMAX_EDGE_BETWEENNESS=16, ALT_ND_PATHMAX_LV_EB=17, ALP_ND_PATHMAX_LV_EB=18)
embedding_methods = enum(RANDOM=1, FARTHEST_D=2, PLANAR=3, BETWEENNESS=4, PAGERANK=5, PAGERANK_MODE=6, PAGERANK_MIN=7, CLOSENESS=8, EDGEBETWEENNESS=9, KATZ=10, LOAD=11, FARTHEST_ECC=13)
em = 0

#returns most granular louvain clustering
# - a pair representing two dictionaries {community_id: [nodes]}, {node: community_id}
def louvain_clustering(g, k):
    global dendo
    #don't recluster for change in embedding method
    if not dendo:
        dendo = generate_dendogram(g)
        print "Crap. Okay, we'll generate the dendogram."
    else:
        print "Oh good. I do not have to generate the dendogram"
    #this was a parameter
    if k > len(dendo):
        return None, None
    #k = len(dendo)
    #comms = partition_at_level(dendo, len(dendo) - int(3*len(dendo)/4))
    comms = partition_at_level(dendo, len(dendo) - k)
    ret_dict = dict()
    #rdk = set(ret_dict.keys())
    for k in comms.keys():
        if comms[k] in ret_dict:
            ret_dict[comms[k]].append(k)
        else:
            ret_dict[comms[k]] = [k]
    #print "Comms:", comms
    #print "Ret_dict:", ret_dict
    return ret_dict, comms

#performs k means clustering on adjacency matrix
# - a pair representing two dictionaries {community_id: [nodes]}, {node: community_id}
def kmeans_clustering(g, k):
    #from sklearn.metrics.pairwise import cosine_similarity
    from sklearn.cluster import KMeans
    ret_dict = {}
    comms = {}
    nodes=g.nodes()
    dict_nodes={nodes[i]:i for i in range(len(nodes))}

    degrees=g.degree()
    X=[]    
    for node in nodes:
        neighbors_node=g.neighbors(node)
        node_nb_vector=[0]*len(nodes)
        for node_nb in neighbors_node:
            node_nb_vector[dict_nodes[node_nb]]=float(1)/(degrees[node_nb])       
        X.append(node_nb_vector)


    X=np.array(X)
    print "Converting to numpy matrix"

    print "Forming cosine similarity matrix"
    #S = cosine_similarity(D) # 100 x 100
    #S[S <= 0.49999] = 0 # using < 0.5 is not working well here for some reason
    #S[S != 0] = 1
    kmeans = KMeans(k)
    print "Executing K-Means on similarity matrix"
    kmeans.fit(X)
    if not kmeans:
        return None, None
    print "Forming partition structure"
    node_list = g.nodes()
    '''for k in comms.keys():
        if comms[k] in ret_dict:
            ret_dict[comms[k]].append(k)
        else:
            ret_dict[comms[k]] = [k]'''
    for i in range(k):
        list_nodes = np.where(kmeans.labels_ == i)[0].tolist()
        #list of nodes in community = ret_dict[community#]
        #community for node = comms[node]
        ret_dict[i] = [node_list[ln] for ln in list_nodes]
        comms.update({node_list[l]:i for l in list_nodes})
    #print ret_dict, comms
    #ret_dict2 = dict((k, v) for k, v in ret_dict.iteritems() if v)
    file = open("kmeans_test.txt", "w")
    file.write(str(ret_dict))
    file.write(str(comms))
    file.close()
    return ret_dict, comms

#performs k means clustering on adjacency matrix
# - a pair representing two dictionaries {community_id: [nodes]}, {node: community_id}
def edge_betweenness_clustering(g, k):
    ret_dict = {}
    comms = {}
    #print "Converting to IGraph"
    ig_g = networkx_to_igraph(g.to_undirected())
    # calculate dendrogram
    #print "Clustering with edge betweenness "
    #clusters = ig_g.community_leading_eigenvector(clusters=k, weights='weight')
    #communities = ig_g.community_edge_betweenness(clusters=k)
    #communities = ig_g.community_edge_betweenness(clusters=k)
    #csteps = int(math.log10(len(g.nodes())))**2 if len(g.nodes()) > 100 else 4
    communities = ig_g.community_fastgreedy()
    #communities = ig_g.community_walktrap(steps=csteps)
    clusters = communities.as_clustering(k)
    node_list = g.nodes()
    for i in range(len(clusters)):
        list_nodes = list(clusters[i])
        #list of nodes in community = ret_dict[community#]
        #community for node = comms[node]
        ret_dict[i] = [node_list[ln] for ln in list_nodes]
        comms.update({node_list[l]:i for l in list_nodes})
    return ret_dict, comms  

def kmean_community_detection_neighbor_information(g, k):
    ret_dict = {}
    comms = {}
    #import Pycluster
    from sklearn.cluster import KMeans
    nodes=g.nodes()
    dict_nodes={nodes[i]:i for i in range(len(nodes))}

    degrees=g.degree()
    X=[]
    def yield_node_nb_vector():
        for node in nodes:
            neighbors_node=g.neighbors(node)
            node_nb_vector=[0]*len(nodes)
            for node_nb in neighbors_node:
                node_nb_vector[dict_nodes[node_nb]]=float(1)/(degrees[node_nb])       
            yield node_nb_vector


    X=np.array([y for y in yield_node_nb_vector()])
    print "Vector space assembled. Executing k-means"
    #labels, error, nfound = Pycluster.kcluster(X, k)
    kmeans = KMeans(k)
    print "Executing K-Means on similarity matrix"
    kmeans.fit(X)
    if not kmeans:
        return None, None
    print "Forming partition structure"
    node_list = g.nodes()
    '''for k in comms.keys():
        if comms[k] in ret_dict:
            ret_dict[comms[k]].append(k)
        else:
            ret_dict[comms[k]] = [k]'''
    for i in range(k):
        list_nodes = np.where(kmeans.labels_ == i)[0].tolist()
        #list of nodes in community = ret_dict[community#]
        #community for node = comms[node]
        ret_dict[i] = [node_list[ln] for ln in list_nodes]
        comms.update({node_list[l]:i for l in list_nodes})
    #print ret_dict, comms
    #ret_dict2 = dict((k, v) for k, v in ret_dict.iteritems() if v)
    file = open("kmeans_test.txt", "w")
    file.write(str(ret_dict))
    file.write(str(comms))
    file.close()
    return ret_dict, comms
partition = {}

def louvain_kmeans(g, k):
    global partition
    print "Executing Louvain K-Means"
    under = g.to_undirected()
    cluster_nodes = {}
    node_clusters = {}
    cnodes = {}
    nclusters = {}
    len_nodes = len(g)
    print "Getting Louvain partition"
    partition = best_partition(under)
    num_coms = max(set(partition.values()))+1
    sub_coms = int(k/num_coms)
    print "Best partition produced", str(num_coms), "partitions. Identifying", str(sub_coms), "partitions within each using K-Means."
    if sub_coms > 2:
        #divide k by number of clusters and use as k for k-means
        delta = 0
        for com in set(partition.values()):
            ret_dict, comms = kmean_community_detection_neighbor_information(g.subgraph([nodes for nodes in partition.keys() if partition[nodes] == com]), sub_coms)
            #list of nodes in community = ret_dict[community#]
            #community for node = comms[node]
            ret_dict2 = dict((int(k)+delta, v) for k, v in ret_dict.iteritems())
            comms2 = dict((k, int(v)+delta) for k, v in comms.iteritems())
            delta = delta + len(ret_dict.keys())
            #get the highest k and set it to delta
            cluster_nodes.update(ret_dict2)
            node_clusters.update(comms2)
        return cluster_nodes, node_clusters
        
    #do louvain regular louvain instead
    for m in range(1, len_nodes):
        cluster_nodes, node_clusters = louvain_clustering(under, m)
        if not cluster_nodes:
            break
        elif len(cluster_nodes.keys()) < k:
            num_landmarks = len(cluster_nodes.keys())
            cnodes = cluster_nodes
            nclusters = node_clusters
    return cnodes, nclusters
    
def louvain_edge_betweenness(g, k):
    return edge_betweenness_clustering(g, k)
    '''global partition
    print "Executing Louvain K-Means"
    under = g.to_undirected()
    cluster_nodes = {}
    node_clusters = {}
    cnodes = {}
    nclusters = {}
    len_nodes = len(g)
    print "Getting Louvain partition"
    partition = best_partition(under)
    num_coms = max(set(partition.values()))+1
    sub_coms = int(k/num_coms)
    
    print "Best partition produced", str(num_coms), "partitions. Identifying", str(sub_coms), "partitions within each using K-Means."
    if sub_coms > 2:
        rem_coms = k % num_coms
        #divide k by number of clusters and use as k for k-means
        delta = 0
        first_k = True
        for com in set(partition.values()):
            if first_k:
                ret_dict, comms = edge_betweenness_clustering(g.subgraph([nodes for nodes in partition.keys() if partition[nodes] == com]), sub_coms+rem_coms)
                first_k = False
            else:
                ret_dict, comms = edge_betweenness_clustering(g.subgraph([nodes for nodes in partition.keys() if partition[nodes] == com]), sub_coms)
            #list of nodes in community = ret_dict[community#]
            #community for node = comms[node]
            ret_dict2 = dict((int(k)+delta, v) for k, v in ret_dict.iteritems())
            comms2 = dict((k, int(v)+delta) for k, v in comms.iteritems())
            delta = delta + len(ret_dict.keys())
            #get the highest k and set it to delta
            cluster_nodes.update(ret_dict2)
            node_clusters.update(comms2)
        #file = open("leb.txt", "w")
        #file.write(str(cluster_nodes))
        #file.write(str(node_clusters))
        #file.close()
        return cluster_nodes, node_clusters
        
    #do louvain regular louvain instead
    for m in range(1, len_nodes):
        cluster_nodes, node_clusters = louvain_clustering(under, m)
        if not cluster_nodes:
            break
        elif len(cluster_nodes.keys()) < k:
            num_landmarks = len(cluster_nodes.keys())
            cnodes = cluster_nodes
            nclusters = node_clusters
    return cnodes, nclusters'''

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
            if val == np.nan:
                return -1
            sys.stdout.flush()
            return val
        except:
            return -1

    a = 1 if h.is_directed() else 0
    b = nx.number_of_nodes(h)
    c = nx.number_of_edges(h)
    d = runf(nx.estrada_index) if b < 10000 and not h.is_directed() else -1
    e = 1 if not g.is_directed() and runf(nx.is_chordal) else 0
    f = runf(nx.graph_clique_number)
    h = runf(nx.graph_number_of_cliques)
    i = runf(nx.transitivity)
    j = runf(nx.average_clustering)
    k = runf(nx.average_node_connectivity) if c < 1000 else -1
    l = runf(nx.edge_connectivity) if c < 1000 else -1
    m = runf(nx.node_connectivity) if c < 1000 else -1
    n = runf(nx.diameter) if c < 1000 else -1
    try:
        o = len(nx.periphery(g)) if b < 1000 else -1
    except:
        o = -1
    p = 1 if runf(nx.is_eulerian) else 0
    q = runf(nx.average_shortest_path_length) if b < 3000 else -1
    try:
        r = nx.connected_double_edge_swap(g, nx.number_of_edges(g)) if b < 5000 else -1
    except:
        r = -1

    s = 1 if runf(nx.is_tree) else 0
    t = runf(nx.density)
    return (a, b, c, d, e, f, h, i, j, k, l, m, n, o, p, q, r, s, t)

#stores the distances from landmarks to subgraph nodes and other landmarks
alt_landmark_dists = {}
#stores the distance from subgraph nodes to their landmarks
alt_landmark_froms = {}
def ALT(s, t):
    global g, alp_landmarks, cluster_nodes, num_alt_landmarks, landmarks,num_alt_calculations, num_alt_estimates, tid,gid, em
    if not landmarks:
        print "Preprocessing ALT"

        pstart = time.time()
        if em == embedding_methods.RANDOM:
            print "ALT Random Landmark Optimized"
            landmarks = random_ls_opt(g, num_alt_landmarks, .00001, 25)
            #landmarks = random_ls_non_opt(g, num_alt_landmarks)
        elif em == embedding_methods.FARTHEST_D or em == embedding_methods.FARTHEST_ECC:
            print "ALT Farthest-D"
            landmarks = farthest_d(g, num_alt_landmarks)
        elif em == embedding_methods.PLANAR:
            print "ALT Planar"
            landmarks = planar(g, num_alt_landmarks)
        elif em == embedding_methods.BETWEENNESS:
            print "Betweenness landmark selection"
            landmarks = sample(alp_betweenness(g, cluster_nodes), num_alt_landmarks)
        elif em == embedding_methods.PAGERANK:
            print "PageRank landmark selection"
            landmarks = sample(alp_pagerank(g, cluster_nodes), num_alt_landmarks)
        elif em == embedding_methods.PAGERANK_MODE:
            print "PageRank Mode landmark selection"
            landmarks = sample(alp_pagerank(g, cluster_nodes, func=mode), num_alt_landmarks)
        elif em == embedding_methods.PAGERANK_MIN:
            print "PageRank Min landmark selection"
            landmarks = sample(alp_pagerank(g, cluster_nodes, func=min), num_alt_landmarks)
        elif em == embedding_methods.CLOSENESS:
            print "ALT Closeness landmark selection"
            landmarks = sample(alp_closeness(g, cluster_nodes), num_alt_landmarks)
        elif em == embedding_methods.KATZ:
            print "Katz landmark selection"
            landmarks = sample(alp_katz(g, cluster_nodes), num_alt_landmarks)
        elif em == embedding_methods.LOAD:
            print "Load centrality landmark selection"
            landmarks = sample(alp_load(g, cluster_nodes), num_alt_landmarks)
        elif em == embedding_methods.FARTHEST_ECC:
            print "Farthest-ECC landmark selection"
            landmarks = farthest_d(g, num_alt_landmarks)
        #num_alt_landmarks = len(landmarks)
        #landmarks = alp_landmarks
        print "ALT Landmarks chosen:", str(landmarks)
        if DEBUG:
            print "Landmarks:", landmarks
        reverse_graph = g.reverse(True)
        def get_landmark_paths(l):
            alt_landmark_dists[l] = single_source_dijkstra_path_length(reverse_graph, l, target_cutoffs=g.nodes())
            alt_landmark_froms[l] = single_source_dijkstra_path_length(g, l, target_cutoffs=g.nodes())

        #landmark_dists is now a matrix for the landmark set
        from multiprocessing.pool import ThreadPool
        pool = ThreadPool(processes=min(int(len(landmarks)/8)+8,12))

        #do 8 landmarks at a time
        pool.map_async(get_landmark_paths, landmarks)
        pool.close()
        pool.join()

        pend = time.time()
        prep_time = pend - pstart

        prep_sql = "INSERT INTO preprocessing VALUES (NULL," + str(tid) + "," + str(heuristics.ALT_ND_PATHMAX_LV_EB) + "," + str(gid) + "," + str(prep_time) + ")"
        if DEBUG:
            print "Prep SQL:" , prep_sql
        query_insert(prep_sql)

        if DEBUG:
            print "Landmarks:", landmarks

    if s == t:
        return 0
        
    num_alt_calculations += len(landmarks)
    num_alt_estimates += 1

    #triangle inequality heuristic
    return max([abs(alt_landmark_froms[l][s] - alt_landmark_dists[l][t]) for l in landmarks])



#landmark ownership dictionary
alp_landmark_dict = {}
#stores the distances from landmarks to subgraph nodes and other landmarks
alp_landmark_dists = {}
#stores the distance from subgraph nodes to their landmarks
alp_landmark_froms = {}
reverse_graph = None
    
def ALP(v, end):
    global g, alp_landmarks, cluster_nodes, node_clusters, num_alp_calculations, num_alp_estimates, tid, gid, em,reverse_graph
    if not alp_landmarks:
        print "Preprocessing ALP"
        pstart = time.time()
        #NOT USING EDGEBETWEENNESS OR DISPERSION
        #embedding_methods = enum(RANDOM=1, FARTHEST_D=2, PLANAR=3, BETWEENNESS=4, PAGERANK=5, PAGERANK_MODE=6,
        # PAGERANK_MIN=7, CLOSENESS=8, KATZ=10, LOAD=11, FARTHEST_DC=13)
        #
        if em == embedding_methods.RANDOM:
            print "Random landmark selection"
            #alp_landmarks = alp_random_ls_non_opt(g, cluster_nodes)
            alp_landmarks = alp_random_ls_opt(g, cluster_nodes, .001, 25)
            #alp_landmarks = alp_random_ls_non_opt(g, cluster_nodes)
        elif em == embedding_methods.FARTHEST_D:
            print "Farthest-D landmark selection"
            alp_landmarks = alp_farthest(g, cluster_nodes)
        elif em == embedding_methods.PLANAR:
            print "Planar landmark selection"
            alp_landmarks = alp_planar(g, cluster_nodes)
        elif em == embedding_methods.BETWEENNESS:
            print "Betweenness landmark selection"
            alp_landmarks = alp_betweenness(g, cluster_nodes)
        elif em == embedding_methods.PAGERANK:
            print "PageRank landmark selection"
            alp_landmarks = alp_pagerank(g, cluster_nodes)
        elif em == embedding_methods.PAGERANK_MODE:
            print "PageRank Mode landmark selection"
            alp_landmarks = alp_pagerank(g, cluster_nodes, func=mode)
        elif em == embedding_methods.PAGERANK_MIN:
            print "PageRank Min landmark selection"
            alp_landmarks = alp_pagerank(g, cluster_nodes, func=min)
        elif em == embedding_methods.CLOSENESS:
            print "ALP Closeness landmark selection"
            alp_landmarks = alp_closeness(g, cluster_nodes)
        elif em == embedding_methods.KATZ:
            print "Katz landmark selection"
            alp_landmarks = alp_katz(g, cluster_nodes)
        elif em == embedding_methods.LOAD:
            print "Load centrality landmark selection"
            alp_landmarks = alp_load(g, cluster_nodes)
        elif em == embedding_methods.FARTHEST_ECC:
            print "Farthest-ECC landmark selection"
            alp_landmarks = alp_farthest_ecc(g, cluster_nodes)
        if not alp_landmarks:
            raise NetworkXError("Could not establish landmarks for this embedding method.")
        if DEBUG:
            print "Landmarks chosen:", alp_landmarks
        else:
            print "Landmarks chosen:", str(len(alp_landmarks))


        landmark_set = alp_landmarks[:]
        #print "Landmarks:", landmarks
        #One final shortest path computation amongst the landmark set
        #landmarks_dists = nx.floyd_warshall_numpy(g, alp_landmarks[:])
        reverse_graph = g.reverse(True)
        
        def get_landmark_paths(l):
            sub = cluster_nodes[node_clusters[l]][:]
            #alp_label = 'ALP_'+str(l)
            alp_landmark_dict.update({s:l for s in sub})
            #distances to the landmark from alp_landmarks+sub because it is on the reverse graph
            alp_landmark_dists[l] = single_source_dijkstra_path_length(reverse_graph, l, target_cutoffs=alp_landmarks[:]+list(sub), only_targets=True)
            alp_landmark_dists[l].update({l:0})
            #distances from the landmark to nodes of its subgraph
            alp_landmark_froms.update(single_source_dijkstra_path_length(g, l, target_cutoffs=sub, only_targets=True))
            alp_landmark_froms.update({l:0})
            #return ({s:l for s in sub},(l,single_source_dijkstra_path_length(rgraph, l, target_cutoffs=alp_landmarks[:]+sub[:], only_targets=True)),single_source_dijkstra_path_length(g, l, target_cutoffs=sub, only_targets=True))
        
        
        #do 10 landmarks at a time
        if len(alp_landmarks) > 2:       
            max_t = min(int(len(alp_landmarks)/8)+8,20)
            #from joblib import Parallel, delayed
            #results = Parallel(n_jobs=max_t,backend="threading")(delayed(get_landmark_paths)(l, reverse_graph) for l in alp_landmarks)

            t_count = 0
            for l in alp_landmarks:
                threading.Thread(target = get_landmark_paths, args=(l,)).start()
                t_count = t_count + 1
                if threading.active_count() > max_t:
                    print "On thread", str(t_count), ". Waiting for ", str(threading.active_count()), " threads to finish before starting next.", str(len(alp_landmarks)-t_count), "remaining."
                    main_thread = threading.currentThread()
                    for t in threading.enumerate():
                        if t is main_thread:
                            continue
                        t.join()
                    collected = gc.collect()
                    #print "Garbage collector: collected %d objects." % (collected)
                    time.sleep(1)
            
            
        else:
            map(get_landmark_paths, alp_landmarks)

        print "Landmark distances calculated. Storing..."

        pend = time.time()
        prep_time = pend - pstart
        print "Landmark trees constructed. Querying graph."
        prep_sql = "INSERT INTO preprocessing VALUES (NULL," + str(tid) + "," + str(heuristics.ALP_ND_PATHMAX_LV_EB) + "," + str(gid) + "," + str(prep_time) + ")"
        query_insert(prep_sql)
        #time to wrap up remaining threads
        
        for t in threading.enumerate():
            main_thread = threading.currentThread()
            if t is main_thread:
                continue
            t.join()
        file = open("landmarks.txt", "w")
    
        for l in alp_landmarks:
            file.write(str(l) + ": " + str(alp_landmark_dists[l]) + "\n")
            #print str(l), ":", str(alp_landmark_dists[l])
        file.close()
    
    if v == end:
        return 0
    '''ALP Heuristics for dual landmark'''
    '''if DEBUG:
        print "V:", str(v)
        print "Node:", str(g.node[v])
        print "End:", str(end)
        print "Node:", str(g.node[end])'''

    #get the landmark for s

    l_1 = alp_landmark_dict[v]
    #get the landmark for t
    l_2 = alp_landmark_dict[end]
        
    l1_v = alp_landmark_froms[v]
    t_l2 = alp_landmark_dists[l_2][end]
    if l_1 == l_2:
        num_alp_calculations += 1
        return abs(l1_v-alp_landmark_dists[l_1][end])

    l2_l1 = alp_landmark_dists[l_1][l_2]

    estimates = [min([g[v][u].get('weight', 1) for u in g[v]]), abs(l1_v-l2_l1)-t_l2, abs(l1_v-t_l2)-l2_l1, abs(l2_l1-t_l2)-l1_v]

    #num_alp_calculations += 6
    estimates.append((abs(l1_v - l2_l1) * abs(l2_l1 - t_l2) - l1_v * t_l2) / l2_l1)
    #if alp_landmark_dists[l_2][l_1] > alp_landmark_dists[l_1][v] + t_l2:
    #v_l1 = alp_landmark_dists[l_1][v]
    #l2_t = alp_landmark_dists[l_2][end]
    #l1_l2 = alp_landmark_dists[l_2][l_1]
    
    #if l2_l1 > v_l1 + t_l2:
    if l2_l1 > l1_v + t_l2:
        num_alp_calculations += 4
        estimates += [abs(l1_v-l2_l1)+abs(t_l2-l2_l1)-l2_l1]
        #estimates += [abs(v_l1-l1_l2)+abs(t_l2-l2_l1)-l2_l1]
    '''else:'''
    '''num_alp_calculations += 6
    if l1_l2 != 0:
        estimates.append((abs(v_l1 - l1_l2) * abs(l1_l2 - l2_t) - v_l1 * l2_t) / l1_l2)'''

    num_alp_estimates += 1
    return max(estimates)


last_num_alp_landmarks = 12
def test(n=1, k=1, debug=False, eid=0, tid=0, gid=0,st_list=[], max_labels=0):
    global cluster_nodes, node_clusters, g,num_alp_calculations,num_alp_estimates,em,landmarks,num_alt_calculations, num_alt_estimates,em, num_alt_landmarks, last_num_alp_landmarks, alp_landmarks
    
    len_nodes = len(g)
    #print "Cluster Nodes:", str(len(node_clusters))
    #get ALT number of landmarks
    print "Getting ALT # Landmarks"
    
    '''if max_labels > 0:
        for k in reversed(range(1, int(len_nodes))):
            if 2*(k * len_nodes) < max_labels:
                num_alt_landmarks = k
                break'''
    num_alt_landmarks = 2**(2+k)
    num_landmarks  = 0
    cnodes = dict()
    nclusters = dict()
    #num_alt_landmarks = 2
    #get ALP number of landmarks
    #print "Creating undirected graph for clustering"
    #under = g.to_undirected()
    #print "Running Edge Betweenness clustering algorithm"
    print "Running K-Means clustering algorithm"
    #k = 10
    '''if max_labels - (2*len_nodes) <= 0:
        return None
    k = int(math.sqrt(max_labels-(2*len_nodes)))
    if k > int(math.sqrt(len_nodes**2-len_nodes)):
        print "Stopping"
        return None'''
    print "Creating", str(k), "clusters"
    cluster_nodes, node_clusters = louvain_clustering(g.to_undirected(), k)
    #cluster_nodes, node_clusters = edge_betweenness_clustering(g,k)
    #cluster_nodes, node_clusters = louvain_edge_betweenness(g, k)
    if not (cluster_nodes and node_clusters):
        return None
    num_landmarks = len(cluster_nodes.keys())
    

    '''for k in range(1, len_nodes):
        cluster_nodes, node_clusters = louvain_clustering(under, k)
        if not cluster_nodes:
            break
        elif len(cluster_nodes.keys())*len(cluster_nodes.keys())+2*len_nodes < max_labels:
            num_landmarks = len(cluster_nodes.keys())
            cnodes = cluster_nodes
            nclusters = node_clusters'''
    
    if num_landmarks == 0 or num_alt_landmarks > num_landmarks:
        print "Num Landmarks:", str(num_landmarks), "Num ALT Landmarks:", str(num_alt_landmarks)
        return
        
    #ALP will have this many landmarks as well
    print "Cluster Nodes:", str(num_landmarks)


    import calendar
    #time since epoch
    print "Running queries..."
    cur_time = calendar.timegm(time.gmtime())
    #sys.stdout.flush()
    def query_graph(i):
        alp_time = 0
        alt_time = 0
        empty_a_star_time = 0
        if not landmarks:
            first = True
        else:
            first = False
        num_alp_calculations = 0
        num_alp_estimates = 0
        s = i[0]
        t = i[1]
        '''if debug:
            print "S:", str(s)
            print "T:", str(t)'''

        '''if debug:
            print "Running A* with no heuristic (Dijkstra's)"'''
        try:
            start = time.time()

            path, astar_size = astar_path(g,s,t, search_space_size=True)
            '''if debug:
                print "A* (no heuristic): ", str(path), "Size:", astar_size'''
            end = time.time()
            empty_a_star_time = end - start
            dij_path_weight = get_path_size(g,path)
            #print str(end - start), "seconds."
            #ALP Query and Error
            num_alp_calculations = 0
            num_alp_estimates = 0
            start = time.time()
            path, alp_size = astar_path_pathmax(g,s,t,ALP, search_space_size=True)
            '''if debug:
                print "A* (ALP)): ",  str(path), "Size:", alp_size, "Number of ALP calculations:", str(num_alp_calculations), "Numbber of ALP Estimates:", str(num_alp_estimates)'''
            end = time.time()
            alp_time = end - start

            alp_path_weight = get_path_size(g,path)
            alp_error = float(abs(alp_path_weight-ALP(s,t)))/float(alp_path_weight)
            '''if debug:
                print(str(empty_a_star_time) + ',' + str(astar_size) + ',' + str(num_alt_calculations) + ',' + str(num_alt_estimates) + ',' + str(alp_time) + ',' + str(alp_size) + ',' + str(num_alp_calculations) + ',' + str(num_alp_estimates) + ',' + str(len(path)) + ',' + str(num_landmarks) + '\n')'''

            num_alt_calculations = 0
            num_alt_estimates = 0
            start = time.time()
            path, alt_size = astar_path_pathmax(g,s,t,ALT, search_space_size=True)
            '''if debug:
                print "A* (ALT)): ", str(path), "Size:", alt_size, "Number of ALT calculations:", str(num_alt_calculations), "Numbber of ALT Estimates:", str(num_alt_estimates)'''
            end = time.time()
            alt_time = end - start

            #print str(end - start), "seconds."
            alt_path_weight = get_path_size(g,path)
            alt_error =  float(abs(alt_path_weight-ALT(s,t)))/float(alt_path_weight)
            '''if debug:
                print(str(empty_a_star_time) + ',' + str(astar_size) + ',' + str(alt_time) + ',' + str(alt_size) + ',' + str(num_alt_calculations) + ',' + str(num_alt_estimates) + ',' + str(alp_time) + ',' + str(alp_size) + ',' + str(num_alp_calculations) + ',' + str(num_alp_estimates) + ',' + str(len(path)) + ',' + str(num_landmarks) + '\n')'''

            if first:
                #otherwise, time results are skewed
                '''if landmarks and alp_landmarks:
                    cluster_nodes = {}
                    node_clusters = {}'''
                first = False
            else:
                a_star_sql = "(NULL," + str(tid) + "," + str(heuristics.DIJKSTRA) + ", NULL," + str(s) + "," + str(t) + "," + str(len(path)) + ", NULL," + str(empty_a_star_time) + "," + str(astar_size) + ",0,0," +str(dij_path_weight) + ")"
                #_star_sql = "(NULL," + str(tid) + "," + str(heuristics.DIJKSTRA) + ", NULL," + str(s) + "," + str(t) + "," + str(len(path)) + ", NULL," + str(empty_a_star_time) + "," + str(astar_size) + ",0,0)"
                alt_sql = "(NULL," + str(tid) + "," + str(heuristics.ALT_ND_PATHMAX_LV_EB) + "," + str(em) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(len(landmarks)) + "," + str(alt_time) + "," + str(alt_size) + "," + str(num_alt_calculations) + "," + str(num_alt_estimates) + "," +str(alt_path_weight) + ")"
                #alp_sql = "(NULL," + str(tid) + "," + str(heuristics.ALP_ND_PATHMAX_LV_EB) + "," + str(em) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(num_landmarks) + "," + str(alp_time) + "," + str(alp_size) + "," + str(num_alp_calculations) + "," + str(num_alp_estimates) + ")"
                #ALP Query and Error
                
                alp_sql = "(NULL," + str(tid) + "," + str(heuristics.ALP_ND_PATHMAX_LV_EB) + "," + str(em) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(num_landmarks) + "," + str(alp_time) + "," + str(alp_size) + "," + str(num_alp_calculations) + "," + str(num_alp_estimates) + "," +str(alp_path_weight) + ")"
                #Dijkstra Query
                a_query_sql = "INSERT INTO query VALUES " + a_star_sql
                query_insert(a_query_sql)

                #ALT Query and Error
                alt_query_sql = "INSERT INTO query VALUES " + alt_sql
                query_id = query_insert(alt_query_sql)
                alt_error_sql = "INSERT INTO error VALUES ("+ str(query_id) + ", " + str(alt_error) + ")"
                query_insert(alt_error_sql)

                #ALP Query and Error
                alp_query_sql = "INSERT INTO query VALUES " + alp_sql
                query_id = query_insert(alp_query_sql)
                alp_error_sql = "INSERT INTO error VALUES (" + str(query_id) + ", " + str(alp_error) + ")"
                query_insert(alp_error_sql)
        except ValueError as e:
            print "Error: Could not run query (S:", str(s), "T:", str(t), "). Selection failed on embedding method", str(em)
            print str(e)
            tend_time = time.time()
            for l in alp_landmarks:
                for v in [n for n in g.nodes() if 'ALP_'+str(l) in g.node[n]]:
                    del g.node[v]['ALP_'+str(l)]
            return
            #exit(1)'''
    
    #landmark_dists is now a matrix for the landmark set
    #setup landmark configuration
    for st in st_list:
        query_graph(st)
    '''query_graph(st_list[0])
    from multiprocessing.pool import ThreadPool
    pool = ThreadPool(processes=8)

    #do 8 landmarks at a time
    pool.map_async(query_graph, st_list[1:])
    pool.close()
    pool.join()'''
    #insert trial
    tend_time = time.time()
    #trial has ended. wipe the graph

    #first ALT
    for l in landmarks:
        for v in [n for n in g.nodes() if 'ALT_'+str(l) in g.node[n]]:
            del g.node[v]['ALT_' + str(l)]

    #now ALP
    for l in alp_landmarks:
        for v in [n for n in g.nodes() if 'ALP_'+str(l) in g.node[n]]:
            del g.node[v]['ALP_' + str(l)]


def get_files_by_file_size(dirname, reverse=False):
    """ Return list of file paths in directory sorted by file size """

    # Get list of files
    filepaths = []
    for basename in os.listdir(dirname):
        filename = os.path.join(dirname, basename)
        if os.path.isfile(filename):
            filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i in xrange(len(filepaths)):
        filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

    # Sort list by file size
    # If reverse=True sort from largest to smallest
    # If reverse=False sort from smallest to largest
    filepaths.sort(key=lambda filename: filename[1], reverse=reverse)

    # Re-populate list with just filenames
    for i in xrange(len(filepaths)):
        filepaths[i] = filepaths[i][0]

    return filepaths

if __name__ == "__main__":
    experiment_id = query_insert("INSERT INTO experiments (description, start_time) values ('ALT vs ALP PathMax (4-points) Arbitrated Landmark Selection', '" + time.strftime('%Y-%m-%d %H:%M:%S') + "')")
    if experiment_id:
        first = False
        second = False
        num_n = 1000

        for i in range(10):
            #choose different graph size (x2)
            #try:
            #choose different graph
            print "Running trial", str(i)

            g_name = "NULL"
            #data_path = r'C:\Users\Acer\Google Drive\Dissertation\conference\STOC_2015\simulation_environment\stoc1\data\exp3'
            #data_path = os.path.join("C:\\","Users","Newton","Google Drive","Dissertation","conference","STOC_2015","simulation_environment","stoc1","data","exp3")
            data_path = "/home/newtonh20/google_drive/Dissertation/conference/STOC_2015/simulation_environment/stoc1/data/large"

            #data_path = "/home/newtonh20/weighted"
            #files = os.listdir(data_path)
            #files.sort()
            files = get_files_by_file_size(data_path)

            #######WINDOWS##############
            #import glob
            #files = glob.glob(data_path + "\\*.*")
            #print files
            ###########################

            try:
                for g_name in files:
                    print "Analyzing graph", g_name
                    #sys.stdout.flush()

                    myg = nx.read_weighted_edgelist(g_name,create_using=nx.DiGraph())
                    g_name = os.path.basename(os.path.normpath(g_name))
                    gs = []
                    if nx.is_strongly_connected(myg):
                        gs = [myg]
                    else:
                        gs = list(nx.strongly_connected_component_subgraphs(myg))


                        #g = max(nx.connected_component_subgraphs(g, False), key=len)
                    for gr in gs:
                        
                        g = gr.copy()
                        try:
                            if len(g.nodes()) < 1000:
                                continue
                            cluster_nodes = None                            
                            dendo = None
                            graph_sql = "SELECT graph_id FROM graphs WHERE source='" + g_name + "' and num_nodes=" + str(len(g.nodes())) + " and num_edges=" + str(len(g.edges()))
                            print nx.info(g)
                            #print graph_sql
                            graph_id = None
                            result = query(graph_sql)
                            print result
                            already = False
                            if result:
                                print "Graph features already analyzed"
                                graph_id = result[0][0]
                                already = True
                            else:
                                print "Analyzing graph features"
                                graph_insert = "INSERT INTO graphs values " + str(get_graph_features(g) + (g_name,))
                                graph_insert = graph_insert.replace("values (", "values (NULL, ")
                                print "Inserting new graph into database:", graph_insert
                                graph_id = query_insert(graph_insert)

                            print "Running trial on", g_name, "graph."
                            if graph_id:
                                g.graph['name'] = graph_id
                                #setup a num_n queries of the same queries
                                st_l = []
                                for c in range(num_n):
                                    s,t = sample(g.nodes(),2)
                                    st_l.append((s,t))
                                o = 1
                                
                                #for o in range(first_em,14):
                                if first:
                                    first = False
                                    kstart = 2
                                else:
                                    kstart = 1
                                for k in range(kstart,5):
                                    if second:
                                        second = False
                                        start = 11
                                    else:
                                        start = 1
                                    #NOT USING EDGEBETWEENNESS OR DISPERSION for o in [1,4,6,7,8,11,13]:
                                    #embedding_methods = enum(RANDOM=1, FARTHEST_D=2, PLANAR=3, BETWEENNESS=4, PAGERANK=5, PAGERANK_MODE=6,
                                    # PAGERANK_MIN=7, CLOSENESS=8, KATZ=10, LOAD=11)
                                    
                                    for o in range(start,14):
                                        if o in [3,9,10,12]:
                                            continue
                                        else:
                                            em = o

                                            try:
                                                trial_sql = "INSERT INTO alt_alp_comparison_trials VALUES (NULL," + str(experiment_id) + "," + str(graph_id) + ")"
                                                print "Trial SQL:", trial_sql
                                                trial_id = query_insert(trial_sql)

                                                if trial_id:
                                                    print "Trial ID:", str(trial_id)
                                                    #choose different landmarks
                                                    #sys.stdout.flush()
                                                    eid = experiment_id
                                                    tid = trial_id
                                                    gid = graph_id
                                                    if eid and tid and gid:
                                                        test(num_n, k, st_list=st_l, eid=experiment_id, tid=trial_id, gid=graph_id)
                                                    else:
                                                        print "Could not run test: ", str(eid), str(tid), str(gid)
                                                    ALP_landmark_dict = dict()
                                                    landmarks = []
                                                    alp_landmarks = []
                                                    num_alp_calculations = 0
                                                    num_alp_estimates = 0
                                                    alp_landmark_dict = {}
                                                    alp_landmark_dists = {}
                                                    alp_landmark_froms = {}
                                                    alt_landmark_dists = {}
                                                    alt_landmark_froms = {}
                                                    sys.exc_clear()
                                                    sys.exc_traceback = sys.last_traceback = None
                                                    
                                                    g.clear()
                                                    del g
                                                    g = gr.copy()
                                                else:
                                                    print "SQL Error: Trouble retrieving last trial query"
                                            except KeyError as e:
                                                print str(e)
                                last_num_alp_landmarks = 0
                            else:
                                print "SQL Error: Trouble retrieving last graph query"
                        except TypeError as e:
                            print str(e)
            except ValueError:
                print "Error during or after graph generation", sys.exc_info()[0]
        #help(r.randint)
