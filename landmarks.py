__author__ = 'Newton'
#Script containing various landmark selection algorithms to speed up shortest path search
from random import choice, sample
import networkx as nx
import path_planning as pp
import itertools

debug = False

# Random landmark selection for ALP
def alp_random_ls_opt(g, subgraphs, trials_to_nodes_in_graph=.1, queries_per_trial=25):
    #perform the landmark selection within each subgraph
    landmarks = [random_ls_opt(g.subgraph(sub).copy(), 1, trials_to_nodes_in_graph, queries_per_trial) for sub in subgraphs.values()]
    return [l[0] for l in landmarks if l]

#random landmark selection for ALP
#g - graph
#c - clusters
def alp_random_ls_opt2(g, subgraphs, node_clusters, trials_to_nodes_in_graph=.1, queries_per_trial=25):
    temp_ls = []
    n_len = nx.number_of_nodes(g)
    def temp_ALP(v, end):
        '''ALP Heuristics for dual landmark'''
        #get the landmark for s
        l_1 = g.node[v]['landmark']
        #get the landmark for t
        l_2 = g.node[end]['landmark']

        estimates = [0, abs(g.node[v]['ALP_'+str(l_1)]-g.node[l_1]['ALP_'+str(l_2)])-g.node[end]['ALP_'+str(l_2)], abs(g.node[v]['ALP_'+str(l_1)]-g.node[end]['ALP_'+str(l_2)])-g.node[l_2]['ALP_'+str(l_1)], abs(g.node[l_2]['ALP_'+str(l_1)]-g.node[end]['ALP_'+str(l_2)])-g.node[v]['ALP_'+str(l_1)]]

        if l_1 == l_2:
            estimates.extend([abs(g.node[v]['ALP_'+str(l_1)]-g.node[end]['ALP_'+str(l_1)]), abs(g.node[v]['ALP_'+str(l_2)]-g.node[end]['ALP_'+str(l_2)])])
        else:
            estimates.extend([(abs(g.node[v]['ALP_'+str(l_1)] - g.node[l_2]['ALP_'+str(l_1)]) * abs(g.node[l_2]['ALP_'+str(l_1)] - g.node[end]['ALP_'+str(l_2)]) - g.node[v]['ALP_'+str(l_1)] * g.node[end]['ALP_'+str(l_2)]) / g.node[l_2]['ALP_'+str(l_1)]])

        return max(estimates)

    min_ls = []
    num_trials = int(n_len * trials_to_nodes_in_graph)

    #grab the one that has the smallest avg search space
    min_search_space = n_len

    #establish a set of test nodes for trial queries such that we can vet each landmark set
    test_nodes = []
    for j in range(queries_per_trial):
        samps = sample(g.nodes(), 2)
        test_nodes.append((samps[0], samps[1]))

    for i in range(num_trials):
        if debug:
            print "ALP Trial ", str(i)
        #choose a set of landmarks
        temp_ls = [choice(subgraphs[c]) for c in subgraphs.keys()]
        #[choice(g.nodes()) for i in range(k)]

        for l in temp_ls:
            #fetch the distances between the landmark and all other landmarks and the distances between the landmark and nodes of its cluster
            nx.set_node_attributes(g, "ALP_"+str(l), pp.single_source_dijkstra_path_length(g, l, target_cutoffs=subgraphs[node_clusters[l]]+temp_ls))

            #label each vertex with the id of its landmark
            for a in subgraphs[node_clusters[l]]:
                g.node[a]['landmark'] = l

        #issue a series of arbitrary shortest path queries using the temp_ALT algorithm as the heuristic
        #compare
        total_search_space = 0
        for j in range(queries_per_trial):
            for (s,t) in test_nodes:
                path, astar_size = pp.astar_path(g,s,t,temp_ALP, search_space_size=True)
                total_search_space += astar_size
            avg_search_space = total_search_space/queries_per_trial
            if avg_search_space < min_search_space:
                min_search_space = avg_search_space
                min_ls = temp_ls

        #clean up graph labeling
        for l in temp_ls:
            for a in subgraphs[node_clusters[l]]:
                if str("ALP_"+str(l)) in g.node[a]:
                    del g.node[a]["ALP_"+str(l)]
                del g.node[a]['landmark']

    return min_ls

def alp_random_ls(subgraphs):
    return [choice(subgraphs[c]) for c in subgraphs.keys()]

# Chooses landmarks randomly for ALT
def random_ls(g, k):
    #[choice(g.nodes()) for i in range(k)]
    return sample(g.nodes(), k)

#Chooses landmarks randomly over a set of trials; chooses the configuration with the highest speedup for ALT
def random_ls_opt(g, k, trials_to_nodes_in_graph=.1, queries_per_trial=25):
    temp_ls = []
    def temp_ALT(v, end):
        return max([abs(g.node[v]['ALT_'+str(l)] - g.node[end]['ALT_'+str(l)]) for l in temp_ls])

    min_ls = []
    num_trials = int(len(g.nodes()) * trials_to_nodes_in_graph)

    #grab the one that has the smallest avg search space
    min_search_space = len(g.nodes())

    #establish a set of test nodes such that we can vet each landmark set
    test_nodes = []
    for j in range(queries_per_trial):
        samps = sample(g.nodes(), 2)
        test_nodes.append((samps[0], samps[1]))

    for i in range(num_trials):
        #choose a set of landmarks
        temp_ls = sample(g.nodes(), k)
        #[choice(g.nodes()) for i in range(k)]

        for l in temp_ls:
            nx.set_node_attributes(g, "ALT_"+str(l), pp.single_source_dijkstra_path_length(g, l))

        #issue a series of arbitrary shortest path queries using the temp_ALT algorithm as the heuristic
        #compare
        total_search_space = 0
        for j in range(queries_per_trial):
            for (s,t) in test_nodes:
                path, astar_size = pp.astar_path(g,s,t,temp_ALT, search_space_size=True)
                total_search_space += astar_size
            avg_search_space = total_search_space/queries_per_trial
            if avg_search_space < min_search_space:
                min_search_space = avg_search_space
                min_ls = temp_ls

        #clean up graph labeling
        for l in temp_ls:
            for n in g.nodes():
                del g.node[n]["ALT_"+str(l)]

    return min_ls

