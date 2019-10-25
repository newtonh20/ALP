__author__ = 'ncampbel'
'''ALP OVER THE SAME SET OF LANDMARKS THAN ALT using different types of graphs recorded in MySQL from files'''
import networkx as nx
from community import *
import random as r
#import matplotlib.pyplot as plt
import heapq
import heapq
import networkx as nx
from networkx.utils import generate_unique_node
from networkx import NetworkXError
from random import choice
import time
import DBAdapter as dba
from memory_profiler import profile
import os
import sys

'''Dijkstra(G =(V,E))'''
def enum(**enums):
    return type('Enum', (object,), enums)

heuristics = enum(DIJKSTRA=1, ALT=2, ALP=3, PCD=4)
embedding_methods = enum(RANDOM=1)
#Warning: Continues forever if no path exists
'''def dijkstra(g, start, end):
    def flatten(L):
        while len(L) > 0:
            yield L[0]
        L = L[1]

    q = [(0, start, ())]
    visited = set()
    while True:
        (cost, v1, path) = heapq.heappop(q)
        if v1 not in visited:
            visited.add(v1)
            if v1 == end:
                return list(flatten(path))[::-1] + [v1]
            path = (v1, path)
            print str(g[v1])
            for (v2, cost2) in g[v1].iteritems():
                if v2 not in visited:
                    heapq.heappush(q, (cost + cost2, v2, path))  # test with n trials'''
# -*- coding: utf-8 -*-
"""
Shortest path algorithms for weighed graphs.
"""
__author__ = """\n""".join(['Aric Hagberg <hagberg@lanl.gov>',
                            'Loic Seguin-C. <loicseguin@gmail.com>',
                            'Dan Schult <dschult@colgate.edu>'])
#    Copyright (C) 2004-2011 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

__all__ = ['dijkstra_path',
           'dijkstra_path_length',
           'bidirectional_dijkstra',
           'single_source_dijkstra',
           'single_source_dijkstra_path',
           'single_source_dijkstra_path_length',
           'all_pairs_dijkstra_path',
           'all_pairs_dijkstra_path_length',
           'dijkstra_predecessor_and_distance',
           'bellman_ford','negative_edge_cycle']

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


def single_source_dijkstra_path_length(G, source, cutoff= None,
                                       weight= 'weight', target_cutoffs=[]):
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
    dist = {}  # dictionary of final distances
    seen = {source:0}
    fringe=[] # use heapq with (distance,label) tuples
    heapq.heappush(fringe,(0,source))
    while fringe:
        (d,v)=heapq.heappop(fringe)
        if v in dist:
            continue # already searched this node.
        dist[v] = d

        if v in target_cutoffs:
            if v in target_cutoffs:
                target_cutoffs.remove(v)

            if not target_cutoffs:
                return dist

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
            if cutoff is not None:
                if vw_dist>cutoff:
                    continue

            if w in dist:
                if vw_dist < dist[w]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?')
            elif w not in seen or vw_dist < seen[w]:
                seen[w] = vw_dist
                heapq.heappush(fringe,(vw_dist,w))

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
            if cutoff is not None:
                if vw_dist>cutoff:
                    continue
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

    if heuristic is None:
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
            while node is not None:
                path.append(node)
                node = explored[node]
            path.reverse()
            if search_space_nodes:
                return path, explored.keys()
            elif search_space_size:
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
            if cutoff is not None:
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


def all_pairs_dijkstra_path_length(G, cutoff=None, weight='weight'):
    """ Compute shortest path lengths between all nodes in a weighted graph.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    distance : dictionary
       Dictionary, keyed by source and target, of shortest path lengths.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length=nx.all_pairs_dijkstra_path_length(G)
    >>> print(length[1][4])
    3
    >>> length[1]
    {0: 1, 1: 0, 2: 1, 3: 2, 4: 3}

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    The dictionary returned only has keys for reachable node pairs.
    """
    paths={}
    for n in G:
        paths[n]=single_source_dijkstra_path_length(G,n, cutoff=cutoff,
                                                    weight=weight)
    return paths

def all_pairs_dijkstra_path(G, cutoff=None, weight='weight'):
    """ Compute shortest paths between all nodes in a weighted graph.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    cutoff : integer or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    distance : dictionary
       Dictionary, keyed by source and target, of shortest paths.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> path=nx.all_pairs_dijkstra_path(G)
    >>> print(path[0][4])
    [0, 1, 2, 3, 4]

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    See Also
    --------
    floyd_warshall()

    """
    paths={}
    for n in G:
        paths[n]=single_source_dijkstra_path(G, n, cutoff=cutoff,
                                             weight=weight)
    return paths

def bellman_ford(G, source, weight = 'weight'):
    """Compute shortest path lengths and predecessors on shortest paths
    in weighted graphs.

    The algorithm has a running time of O(mn) where n is the number of
    nodes and m is the number of edges.  It is slower than Dijkstra but
    can handle negative edge weights.

    Parameters
    ----------
    G : NetworkX graph
       The algorithm works for all types of graphs, including directed
       graphs and multigraphs.

    source: node label
       Starting node for path

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    Returns
    -------
    pred, dist : dictionaries
       Returns two dictionaries keyed by node to predecessor in the
       path and to the distance from the source respectively.

    Raises
    ------
    NetworkXUnbounded
       If the (di)graph contains a negative cost (di)cycle, the
       algorithm raises an exception to indicate the presence of the
       negative cost (di)cycle.  Note: any negative weight edge in an
       undirected graph is a negative cost cycle.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(5, create_using = nx.DiGraph())
    >>> pred, dist = nx.bellman_ford(G, 0)
    >>> pred
    {0: None, 1: 0, 2: 1, 3: 2, 4: 3}
    >>> dist
    {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}

    >>> from nose.tools import assert_raises
    >>> G = nx.cycle_graph(5, create_using = nx.DiGraph())
    >>> G[1][2]['weight'] = -7
    >>> assert_raises(nx.NetworkXUnbounded, nx.bellman_ford, G, 0)

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    The dictionaries returned only have keys for nodes reachable from
    the source.

    In the case where the (di)graph is not connected, if a component
    not containing the source contains a negative cost (di)cycle, it
    will not be detected.

    """
    if source not in G:
        raise KeyError("Node %s is not found in the graph"%source)
    numb_nodes = len(G)

    dist = {source: 0}
    pred = {source: None}

    if numb_nodes == 1:
       return pred, dist

    if G.is_multigraph():
        def get_weight(edge_dict):
            return min([eattr.get(weight,1) for eattr in edge_dict.values()])
    else:
        def get_weight(edge_dict):
            return edge_dict.get(weight,1)

    for i in range(numb_nodes):
        no_changes=True
        # Only need edges from nodes in dist b/c all others have dist==inf
        for u, dist_u in list(dist.items()): # get all edges from nodes in dist
            for v, edict in G[u].items():  # double loop handles undirected too
                dist_v = dist_u + get_weight(edict)
                if v not in dist or dist[v] > dist_v:
                    dist[v] = dist_v
                    pred[v] = u
                    no_changes = False
        if no_changes:
            break
    else:
        raise nx.NetworkXUnbounded("Negative cost cycle detected.")
    return pred, dist

def negative_edge_cycle(G, weight = 'weight'):
    """Return True if there exists a negative edge cycle anywhere in G.

    Parameters
    ----------
    G : NetworkX graph

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    Returns
    -------
    negative_cycle : bool
        True if a negative edge cycle exists, otherwise False.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.cycle_graph(5, create_using = nx.DiGraph())
    >>> print(nx.negative_edge_cycle(G))
    False
    >>> G[1][2]['weight'] = -7
    >>> print(nx.negative_edge_cycle(G))
    True

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    This algorithm uses bellman_ford() but finds negative cycles
    on any component by first adding a new node connected to
    every node, and starting bellman_ford on that node.  It then
    removes that extra node.
    """
    newnode = generate_unique_node()
    G.add_edges_from([ (newnode,n) for n in G])

    try:
        bellman_ford(G, newnode, weight)
    except nx.NetworkXUnbounded:
        G.remove_node(newnode)
        return True
    G.remove_node(newnode)
    return False


def bidirectional_dijkstra(G, source, target, weight = 'weight'):
    """Dijkstra's algorithm for shortest paths using bidirectional search.

    Parameters
    ----------
    G : NetworkX graph

    source : node
       Starting node.

    target : node
       Ending node.

    weight: string, optional (default='weight')
       Edge data key corresponding to the edge weight

    Returns
    -------
    length : number
        Shortest path length.

    Returns a tuple of two dictionaries keyed by node.
    The first dictionary stores distance from the source.
    The second stores the path from the source to that node.

    Raises
    ------
    NetworkXNoPath
        If no path exists between source and target.

    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length,path=nx.bidirectional_dijkstra(G,0,4)
    >>> print(length)
    4
    >>> print(path)
    [0, 1, 2, 3, 4]

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    In practice  bidirectional Dijkstra is much more than twice as fast as
    ordinary Dijkstra.

    Ordinary Dijkstra expands nodes in a sphere-like manner from the
    source. The radius of this sphere will eventually be the length
    of the shortest path. Bidirectional Dijkstra will expand nodes
    from both the source and the target, making two spheres of half
    this radius. Volume of the first sphere is pi*r*r while the
    others are 2*pi*r/2*r/2, making up half the volume.

    This algorithm is not guaranteed to work if edge weights
    are negative or are floating point numbers
    (overflows and roundoff errors can cause problems).

    See Also
    --------
    shortest_path
    shortest_path_length
    """
    if source is None or target is None:
        raise nx.NetworkXException(
            "Bidirectional Dijkstra called with no source or target")
    if source == target: return (0, [source])
    #Init:   Forward             Backward
    dists =  [{},                {}]# dictionary of final distances
    paths =  [{source:[source]}, {target:[target]}] # dictionary of paths
    fringe = [[],                []] #heap of (distance, node) tuples for extracting next node to expand
    seen =   [{source:0},        {target:0} ]#dictionary of distances to nodes seen
    #initialize fringe heap
    heapq.heappush(fringe[0], (0, source))
    heapq.heappush(fringe[1], (0, target))
    #neighs for extracting correct neighbor information
    if G.is_directed():
        neighs = [G.successors_iter, G.predecessors_iter]
    else:
        neighs = [G.neighbors_iter, G.neighbors_iter]
    #variables to hold shortest discovered path
    #finaldist = 1e30000
    finalpath = []
    dir = 1
    while fringe[0] and fringe[1]:
        # choose direction
        # dir == 0 is forward direction and dir == 1 is back
        dir = 1-dir
        # extract closest to expand
        (dist, v )= heapq.heappop(fringe[dir])
        if v in dists[dir]:
            # Shortest path to v has already been found
            continue
        # update distance
        dists[dir][v] = dist #equal to seen[dir][v]
        if v in dists[1-dir]:
            # if we have scanned v in both directions we are done
            # we have now discovered the shortest path
            return (finaldist,finalpath)

        for w in neighs[dir](v):
            if(dir==0): #forward
                if G.is_multigraph():
                    minweight=min((dd.get(weight,1)
                               for k,dd in G[v][w].items()))
                else:
                    minweight=G[v][w].get(weight,1)
                vwLength = dists[dir][v] + minweight #G[v][w].get(weight,1)
            else: #back, must remember to change v,w->w,v
                if G.is_multigraph():
                    minweight=min((dd.get(weight,1)
                               for k,dd in G[w][v].items()))
                else:
                    minweight=G[w][v].get(weight,1)
                vwLength = dists[dir][v] + minweight #G[w][v].get(weight,1)

            if w in dists[dir]:
                if vwLength < dists[dir][w]:
                    raise ValueError("Contradictory paths found: negative weights?")
            elif w not in seen[dir] or vwLength < seen[dir][w]:
                # relaxing
                seen[dir][w] = vwLength
                heapq.heappush(fringe[dir], (vwLength,w))
                paths[dir][w] = paths[dir][v]+[w]
                if w in seen[0] and w in seen[1]:
                    #see if this path is better than than the already
                    #discovered shortest path
                    totaldist = seen[0][w] + seen[1][w]
                    if finalpath == [] or finaldist > totaldist:
                        finaldist = totaldist
                        revpath = paths[1][w][:]
                        revpath.reverse()
                        finalpath = paths[0][w] + revpath[1:]
    raise nx.NetworkXNoPath("No path between %s and %s." % (source, target))

def random_ls(g, k):
    return [choice(g.nodes()) for i in range(k)]

#random landmark selection for ALP
#g - graph
#c - clusters
def alp_random_ls(subgraphs):
    return [choice(subgraphs[c]) for c in subgraphs.keys()]

def independent_set_ls(g,k):
    return nx.maximal_independent_set(g)

#returns most granular louvain clustering
# - a pair representing two dictionaries {community_id: [nodes]}, {node: community_id}
def louvain_clustering(g):
    dendo = generate_dendogram(g)
    comms = partition_at_level(dendo, len(dendo) - 1)
    ret_dict = dict()
    for k in comms.keys():
        if comms[k] in ret_dict.keys():
            ret_dict[comms[k]].append(k)
        else:
            ret_dict[comms[k]] = [k]
    #print "Comms:", comms
    #print "Ret_dict:", ret_dict
    return ret_dict, comms

#g - graph
#l - a set of landmarks
#k - distance to step out
#returns a set of overlapping clusters with radius k
# - a pair representing two dictionaries {community_id: [nodes]}, {node: community_id}
def k_clustering(g, l, k=1):
    def get_k_neighbors(n, i=0, li=[]):
        if i == k:
            return n
        else:
            i += 1
            li.extend([get_k_neighbors(ne, i, li) for ne in g.neighbors(n)])

    ret_dict = dict()
    for j in range(len(l)):
        ret_dict[j] = get_k_neighbors(l[j])

    nodes = []
    for v in ret_dict.values():
        nodes.extend(v)

    if not set(g.nodes).difference(set(nodes)):
        return ret_dict
    else:
        return k_clustering(g, l, k+1)

def dist(a, b):
    return 1
    '''(x1, y1) = a
    (x2, y2) = b
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5'''

def get_graph_features(g):
    def runf(func):
        try:
            print "Running", str(func)
            val = func(g)
            print "Value:", str(val)
            sys.stdout.flush()
            return val
        except:
            return -1

    a = 1 if runf(nx.is_directed) else 0
    b = runf(nx.number_of_nodes)
    c = runf(nx.number_of_edges)
    d = runf(nx.estrada_index) if b < 1000 else -1
    e = 1 if runf(nx.is_chordal) else 0
    f = runf(nx.graph_clique_number)
    h = runf(nx.graph_number_of_cliques)
    i = runf(nx.transitivity)
    j = runf(nx.average_clustering)
    k = runf(nx.average_node_connectivity) if b < 10000 else -1
    l = runf(nx.edge_connectivity) if b < 10000 else -1
    m = runf(nx.node_connectivity) if b < 10000 else -1
    n = runf(nx.diameter) if b < 1000 else -1
    try:
        o = len(nx.periphery(g)) if b < 1000 else -1
    except:
        o = -1
    p = 1 if runf(nx.is_eulerian) else 0
    q = runf(nx.average_shortest_path_length) if b < 10000 else -1
    try:
        r = nx.connected_double_edge_swap(g, nx.number_of_edges(g)) if b < 10000 else -1
    except:
        r = -1

    s = 1 if runf(nx.is_tree) else 0
    t = runf(nx.density)
    return (a, b, c, d, e, f, h, i, j, k, l, m, n, o, p, q, r, s, t)

#def dijkstra(G,source,target):
#    return bidirectional_dijkstra(G,source,target)
ALP_landmark_dict = dict()
landmarks = []
alp_landmarks = []
num_alt_calculations = 0
num_alp_calculations = 0
num_alt_estimates = 0
num_alp_estimates = 0
def test(n=1, num_nodes=500000, g=None, debug=True, eid=0, tid=0, gid=0):
    if not g:
        g = nx.barabasi_albert_graph(num_nodes, 5)
    #g = nx.erdos_renyi_graph(1000000, .1)
    #g = nx.random_lobster(100000, .4, .4)
    #g = nx.random_powerlaw_tree(10000)
    #g = nx.powerlaw_cluster_graph(100000, 10, .2)




    #print str(g.adj)
    #nx.draw(g)
    #plt.show()

    #num_landmarks = 20

    print "Running Louvain algorithm"
    cluster_nodes, node_clusters = louvain_clustering(g)
    #print "Cluster Nodes:", str(len(node_clusters))

    #ALP will have this many landmarks as well
    num_landmarks = len(cluster_nodes.keys())
    print "Cluster Nodes:", str(num_landmarks)
    global num_alt_calculations
    global num_alp_calculations
    global num_alt_estimates
    global num_alp_estimates

    def ALT(s, t):
        global landmarks
        global num_alt_calculations
        global num_alt_estimates
        if not landmarks:
            print "Preprocessing"

            pstart = time.time()
            #run preprocessing algorithm
            #if alp_landmarks:
            #landmarks = alp_landmarks
            #else:
            landmarks = random_ls(g, choice(range(4,20)))
            print "ALT Landmarks chosen:", str(landmarks)
            #print "Landmarks:", landmarks
            for l in landmarks:
                nx.set_node_attributes(g, "ALT_"+str(l), single_source_dijkstra_path_length(g, l))
            pend = time.time()
            prep_time = pend - pstart
            #print "ALT Dict:", ALT_dict
            prep_sql = "INSERT INTO preprocessing VALUES (NULL," + str(tid) + "," + str(heuristics.ALT) + "," + str(gid) + "," + str(prep_time) + ")"
            dba.query_insert(prep_sql)

        #print "Landmarks:", landmarks
        #triangle inequality heuristic
        num_alt_calculations += len(landmarks)
        num_alt_estimates += 1
        return max([abs(g.node[s]['ALT_'+str(l)] - g.node[t]['ALT_'+str(l)]) for l in landmarks])

    def ALP(v, end):
        global alp_landmarks
        global num_alp_calculations
        global num_alp_estimates
        if not alp_landmarks:
            print "Preprocessing"

            pstart = time.time()
            alp_landmarks = alp_random_ls(cluster_nodes)
            print "Landmarks chosen:", alp_landmarks
            #print "Landmarks:", landmarks
            for l in alp_landmarks:
                #fetch the distances between the landmark and all other landmarks and the distances between the landmark and nodes of its cluster
                #print "Target cutoffs for", str(l), ":", cluster_nodes[node_clusters[l]]+alp_landmarks
                nx.set_node_attributes(g, "ALP_"+str(l), single_source_dijkstra_path_length(g, l, target_cutoffs=cluster_nodes[node_clusters[l]]+alp_landmarks))
                #print "Known lengths from", str(l), ":", ALP_dict[l]

                #label each vertex with the id of its landmark
                for a in cluster_nodes[node_clusters[l]]:
                    g.node[a]['landmark'] = l

            pend = time.time()
            prep_time = pend - pstart
            #print "ALT Dict:", ALT_dict
            prep_sql = "INSERT INTO preprocessing VALUES (NULL," + str(tid) + "," + str(heuristics.ALP) + "," + str(gid) + "," + str(prep_time) + ")"
            dba.query_insert(prep_sql)

        #print "Landmarks:", landmarks
        #triangle inequality heuristic

        '''ALP Heuristics for dual landmark'''
        #get the landmark for s
        l_1 = g.node[v]['landmark']
        #get the landmark for t
        l_2 = g.node[end]['landmark']
        '''print "ALP_dict[" + str(l_1) + "]:", ALP_dict[l_1]
        print "ALP_dict[" + str(l_2) + "]:", ALP_dict[l_2]
        print "ALP_dict[l_1][v]:", str(ALP_dict[l_1][v])
        print "ALP_dict[l_1][l_2]:", str(ALP_dict[l_1][l_2])
        print "ALP_dict[l_2][end]:", str(ALP_dict[l_2][end])'''

        estimates = [0, abs(g.node[v]['ALP_'+str(l_1)]-g.node[l_1]['ALP_'+str(l_2)])-g.node[end]['ALP_'+str(l_2)], abs(g.node[v]['ALP_'+str(l_1)]-g.node[end]['ALP_'+str(l_2)])-g.node[l_2]['ALP_'+str(l_1)], abs(g.node[l_2]['ALP_'+str(l_1)]-g.node[end]['ALP_'+str(l_2)])-g.node[v]['ALP_'+str(l_1)]]

        num_alp_calculations += 6
        if l_1 == l_2:
            estimates.extend([abs(g.node[v]['ALP_'+str(l_1)]-g.node[end]['ALP_'+str(l_1)]), abs(g.node[v]['ALP_'+str(l_2)]-g.node[end]['ALP_'+str(l_2)])])
            num_alp_calculations += 2
        else:
            estimates.extend([(abs(g.node[v]['ALP_'+str(l_1)] - g.node[l_2]['ALP_'+str(l_1)]) * abs(g.node[l_2]['ALP_'+str(l_1)] - g.node[end]['ALP_'+str(l_2)]) - g.node[v]['ALP_'+str(l_1)] * g.node[end]['ALP_'+str(l_2)]) / g.node[l_2]['ALP_'+str(l_1)]])
            num_alp_calculations += 6

        num_alp_estimates += 1
        return max(estimates)

        ##########for debugging#########################################################
        '''ALP_1 = abs(ALP_dict[l_1][v]-ALP_dict[l_1][l_2])-ALP_dict[l_2][end]
        ALP_2 = abs(ALP_dict[l_1][v]-ALP_dict[l_2][end])-ALP_dict[l_1][l_2]
        ALP_3 = abs(ALP_dict[l_1][l_2]-ALP_dict[l_2][end])-ALP_dict[l_1][v]
        ALP_4 = abs(ALP_dict[l_1][v]-ALP_dict[l_1][end]) if end in ALP_dict[l_1] else 0
        ALP_5 = abs(ALP_dict[l_2][v]-ALP_dict[l_2][end]) if v in ALP_dict[l_2] else 0
        ptolemys = (abs(ALP_dict[l_1][v]-ALP_dict[l_1][l_2])*abs(ALP_dict[l_1][l_2]-ALP_dict[l_2][end])-ALP_dict[l_1][v]*ALP_dict[l_2][end])/ALP_dict[l_1][l_2] if l_1 != l_2 else 0
        vals = [ALP_1, ALP_2, ALP_3, ALP_4, ALP_5, ptolemys]
        chosen = max(vals)
        #data for mysql
        #print "ALP_1:", str(ALP_1), ", ALP_2:", str(ALP_2), ", ALP_3:", str(ALP_3), ", ALP_4:", str(ALP_4), ", ALP_5:", str(ALP_5), ", Ptolemy's:", str(ptolemys)
        return chosen'''

    first = True
    import calendar
    #time since epoch
    print "Running queries..."
    cur_time = calendar.timegm(time.gmtime())
    sys.stdout.flush()
    for i in range(n):
        num_alp_calculations = 0
        num_alp_estimates = 0
        s = choice(g.nodes())
        t = choice(g.nodes())
        if debug:
            print "S:", str(s)
            print "T:", str(t)
        '''start = time.time()
        print "Dijkstra:", dijkstra_path(g, s, t)
        end = time.time()
        dijkstra_time = end - start
        print str(end - start), "seconds."'''
        if debug:
            print "Running A* with no heuristic (Dijkstra's)"
        start = time.time()
        path, astar_size = astar_path(g,s,t, search_space_size=True)
        if debug:
            print "A* (no heuristic): ", str(path), "Size:", astar_size
        end = time.time()
        empty_a_star_time = end - start
        #print str(end - start), "seconds."
        num_alp_calculations = 0
        num_alp_estimates = 0
        start = time.time()
        path, alp_size = astar_path(g,s,t,ALP, search_space_size=True)
        if debug:
            print "A* (ALP)): ",  str(path), "Size:", alp_size, "Number of ALP calculations:", str(num_alp_calculations), "Numbber of ALP Estimates:", str(num_alp_estimates)
        end = time.time()
        alp_time = end - start
        #print str(end - start), "seconds."
        num_alt_calculations = 0
        num_alt_estimates = 0
        start = time.time()
        path, alt_size = astar_path(g,s,t,ALT, search_space_size=True)
        if debug:
            print "A* (ALT)): ", str(path), "Size:", alt_size, "Number of ALT calculations:", str(num_alt_calculations), "Numbber of ALT Estimates:", str(num_alt_estimates)
        end = time.time()
        alt_time = end - start
        #print str(end - start), "seconds."
        if debug:
            print(str(empty_a_star_time) + ',' + str(astar_size) + ',' + str(alt_time) + ',' + str(alt_size) + ',' + str(num_alt_calculations) + ',' + str(num_alt_estimates) + ',' + str(alp_time) + ',' + str(alp_size) + ',' + str(num_alp_calculations) + ',' + str(num_alp_estimates) + ',' + str(len(path)) + ',' + str(num_landmarks) + '\n')
        if first:
            first = False
        else:
            a_star_sql = "(NULL," + str(tid) + "," + str(heuristics.DIJKSTRA) + ", NULL," + str(s) + "," + str(t) + "," + str(len(path)) + ", NULL," + str(empty_a_star_time) + "," + str(astar_size) + ",0,0)"
            alt_sql = "(NULL," + str(tid) + "," + str(heuristics.ALT) + "," + str(embedding_methods.RANDOM) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(len(landmarks)) + "," + str(alt_time) + "," + str(alt_size) + "," + str(num_alt_calculations) + "," + str(num_alt_estimates) + ")"
            alp_sql = "(NULL," + str(tid) + "," + str(heuristics.ALP) + "," + str(embedding_methods.RANDOM) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(num_landmarks) + "," + str(alp_time) + "," + str(alp_size) + "," + str(num_alp_calculations) + "," + str(num_alp_estimates) + ")"
            query_sql = "INSERT INTO query VALUES " + a_star_sql + ", " + alt_sql + ", " + alp_sql
            if debug:
                print "QUERY SQL:", query_sql
            dba.query_insert(query_sql)
    #insert trial
    tend_time = time.time()
    #trial has ended. wipe the graph
    #start with ALT
    for v in g.nodes():
        for l in landmarks:
            if 'ALT_'+str(l) in g.node[v]:
                del g.node[v]['ALT_'+str(l)]
    #now ALP
    for l in alp_landmarks:
        for v in [n for n in g.nodes() if 'ALP_'+str(l) in g.node[n]]:
            del g.node[v]['ALP_' + str(l)]

dba.query_insert("INSERT INTO experiments (description, start_time) values ('US Road graphs different landmarks', '" + time.strftime('%Y-%m-%d %H:%M:%S') + "')")
result = dba.query("SELECT MAX(experiment_id) FROM experiments")
if result:
    experiment_id = result[0][0]
    num_n = 2000
    for i in range(10):
        #choose different graph size (x2)
        for j in range(10):
            #try:
            #choose different graph
            print "Running trial", str(j)

            g_name = "NULL"
            data_path = "/home/newtonh20/data/roadNet"
            #data_path = "/home/newtonh20/trial2"
            files = os.listdir(data_path)
            files.sort()
            for g_name in files:
                print "Analyzing graph", g_name
                sys.stdout.flush()

                g = nx.read_edgelist(data_path + "/" + g_name)
                if not nx.is_connected(g):
                    g = list(nx.connected_component_subgraphs(g, False))[0]
                    #g = max(nx.connected_component_subgraphs(g, False), key=len)

                result = dba.query("SELECT graph_id FROM graphs WHERE source='" + g_name + "'")
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
                    dba.query_insert(graph_insert)

                sys.stdout.flush()
                #print graph_insert
                print "Running trial on", g_name, "graph."
                sys.stdout.flush()
                result = dba.query("SELECT MAX(graph_id) as id FROM graphs")
                if result:
                    graph_id = result[0][0] if not already else graph_id
                    print "Graph ID:", str(graph_id)
                    for k in range(10):
                        trial_sql = "INSERT INTO alt_alp_comparison_trials VALUES (NULL," + str(experiment_id) + "," + str(graph_id) + ")"
                        #print "Trial SQL:", trial_sql
                        dba.query_insert(trial_sql)
                        result = dba.query("SELECT MAX(trial_id) as id FROM alt_alp_comparison_trials")
                        if result:
                            trial_id = result[0][0]
                            print "Trial ID:", str(trial_id)
                            #choose different landmarks
                            sys.stdout.flush()
                            test(num_n, nx.number_of_nodes(g), g, debug=False, eid=experiment_id, tid=trial_id, gid=graph_id)
                            ALP_landmark_dict = dict()
                            landmarks = []
                            alp_landmarks = []
                            num_alt_calculations = 0
                            num_alp_calculations = 0
                            num_alt_estimates = 0
                            num_alp_estimates = 0
                        else:
                            print "SQL Error: Trouble retrieving last trial query"
                else:
                    print "SQL Error: Trouble retrieving last graph query"
            '''except TypeError as e:
                print str(e)
            except:
                print "Error during or after graph generation", sys.exc_info()[0]'''
        num_n = num_n * 10
    #help(r.randint)
