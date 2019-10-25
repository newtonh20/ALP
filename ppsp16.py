__author__ = 'Newton'
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
from landmarks import *
import path_planning as pp

test = False
DEBUG = False
debug = False
g = None
ALP_landmark_dict = dict()
landmarks = []
alp_landmarks = []
num_alt_calculations = 0
num_alp_calculations = 0
num_alt_estimates = 0
cluster_nodes = None
num_alp_estimates = 0

gid = 0
tid = 0
eid = 0

'''Dijkstra(G =(V,E))'''
def enum(**enums):
    return type('Enum', (object,), enums)

heuristics = enum(DIJKSTRA=1, ALT=2, ALP=3, PCD=4, ALT_OPT=5, ALP_OPT=6)
embedding_methods = enum(RANDOM=1)
#Warning: Continues forever if no path exists

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

# Executes P2P shortest path algorithm between source and target using ALT
def ALT(s, t):
    global g, landmarks, num_alt_calculations, num_alt_estimates, tid,gid

    if not landmarks:
        print "Preprocessing ALT"

        pstart = time.time()
        #run preprocessing algorithm
        #if alp_landmarks:
        #landmarks = alp_landmarks
        #else:
        landmarks = random_ls(g, choice(range(4,20)))
        print "ALT Landmarks chosen:", str(landmarks)
        #print "Landmarks:", landmarks
        for l in landmarks:
            nx.set_node_attributes(g, "ALT_"+str(l), pp.single_source_dijkstra_path_length(g, l))
        pend = time.time()
        prep_time = pend - pstart
        #print "ALT Dict:", ALT_dict
        if not test:
            prep_sql = "INSERT INTO preprocessing VALUES (NULL," + str(tid) + "," + str(heuristics.ALT_OPT) + "," + str(gid) + "," + str(prep_time) + ")"
            dba.query_insert(prep_sql)

        if DEBUG:
            print "Landmarks:", landmarks

    num_alt_calculations += len(landmarks)
    num_alt_estimates += 1

    #triangle inequality heuristic
    return max([abs(g.node[s]['ALT_'+str(l)] - g.node[t]['ALT_'+str(l)]) for l in landmarks])

# Executes P2P shortest path algorithm between source and target using ALP
def ALP(v, end):
    global g, alp_landmarks, cluster_nodes, node_clusters, num_alp_calculations, num_alp_estimates

    if not alp_landmarks:
        print "Preprocessing ALP"

        pstart = time.time()

        alp_landmarks = alp_random_ls_opt(g, cluster_nodes)
        print "Landmarks chosen:", alp_landmarks
        #print "Landmarks:", landmarks
        for l in alp_landmarks:
            #fetch the distances between the landmark and all other landmarks and the distances between the landmark and nodes of its cluster
            nx.set_node_attributes(g, "ALP_"+str(l), pp.single_source_dijkstra_path_length(g, l, target_cutoffs=cluster_nodes[node_clusters[l]]+alp_landmarks))
            #print "Known lengths from", str(l), ":", ALP_dict[l]

            #label each vertex with the id of its landmark
            for a in cluster_nodes[node_clusters[l]]:
                g.node[a]['landmark'] = l

        pend = time.time()
        prep_time = pend - pstart
        #print "ALT Dict:", ALT_dict
        if not test:
            prep_sql = "INSERT INTO preprocessing VALUES (NULL," + str(tid) + "," + str(heuristics.ALP_OPT) + "," + str(gid) + "," + str(prep_time) + ")"
            dba.query_insert(prep_sql)

    #print "Landmarks:", landmarks
    #triangle inequality heuristic

    '''ALP Heuristics for dual landmark'''
    #get the landmark for s
    l_1 = g.node[v]['landmark']
    #get the landmark for t
    l_2 = g.node[end]['landmark']

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


def test(n=1, num_nodes=500000, debug=True, eid=0, tid=0, gid=0):
    global cluster_nodes, node_clusters, g
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

        if debug:
            print "Running A* with no heuristic (Dijkstra's)"
        start = time.time()
        path, astar_size = pp.astar_path(g,s,t, search_space_size=True)
        if debug:
            print "A* (no heuristic): ", str(path), "Size:", astar_size
        end = time.time()
        empty_a_star_time = end - start
        #print str(end - start), "seconds."
        num_alp_calculations = 0
        num_alp_estimates = 0
        start = time.time()
        path, alp_size = pp.astar_path(g,s,t,ALP, search_space_size=True)
        if debug:
            print "A* (ALP)): ",  str(path), "Size:", alp_size, "Number of ALP calculations:", str(num_alp_calculations), "Numbber of ALP Estimates:", str(num_alp_estimates)
        end = time.time()
        alp_time = end - start
        #print str(end - start), "seconds."
        num_alt_calculations = 0
        num_alt_estimates = 0
        start = time.time()
        path, alt_size = pp.astar_path(g,s,t,ALT, search_space_size=True)
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
            alt_sql = "(NULL," + str(tid) + "," + str(heuristics.ALT_OPT) + "," + str(embedding_methods.RANDOM) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(len(landmarks)) + "," + str(alt_time) + "," + str(alt_size) + "," + str(num_alt_calculations) + "," + str(num_alt_estimates) + ")"
            alp_sql = "(NULL," + str(tid) + "," + str(heuristics.ALP_OPT) + "," + str(embedding_methods.RANDOM) + "," + str(s) + "," + str(t) + "," + str(len(path)) + "," + str(num_landmarks) + "," + str(alp_time) + "," + str(alp_size) + "," + str(num_alp_calculations) + "," + str(num_alp_estimates) + ")"
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

dba.query_insert("INSERT INTO experiments (description, start_time) values ('Foreign Road graphs different landmarks', '" + time.strftime('%Y-%m-%d %H:%M:%S') + "')")
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
            data_path = "/home/newtonh20/gdrive/Dissertation/conference/STOC_2015/simulation_environment/stoc1/data/roadNet"
            #data_path = "/home/newtonh20/weighted"
            files = os.listdir(data_path)
            files.sort()
            for g_name in files:
                print "Analyzing graph", g_name
                sys.stdout.flush()

                g = nx.read_weighted_edgelist(data_path + "/" + g_name)
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
        num_n = num_n * 2
    #help(r.randint)
