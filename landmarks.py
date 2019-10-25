__author__ = 'Acer'

from random import choice, sample
import networkx as nx
import path_planning as pp
import operator
import itertools
import collections
#import matplotlib.pyplot as plt

#TODO: Implement this in C using http://www.cs.yale.edu/homes/aspnes/pinewiki/C(2f)Graphs.html
debug = False
DEBUG = False
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
        landmarks = [random_ls_opt(list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0], 1, trials_to_nodes_in_graph, queries_per_trial) for sub in subgraphs.values()]
    else:
        landmarks = [random_ls_opt(list(nx.connected_component_subgraphs(g.subgraph(sub), False))[0], 1, trials_to_nodes_in_graph, queries_per_trial) for sub in subgraphs.values()]
    return [l[0] for l in landmarks if l]

#performs random landmark selection for ALP without trials
def alp_random_ls_non_opt(g,subgraphs):
    #perform the landmark selection within each subgraph
    if g.is_directed():
        landmarks = [random_ls_non_opt(list(nx.strongly_connected_component_subgraphs(g.subgraph(sub)))[0], 1) for sub in subgraphs.values()]
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
def random_ls_opt(g, k, trials_to_nodes_in_graph=.1, queries_per_trial=25):
    if k < 1:
        raise ValueError("Cannot specify " + str(k) + " as number of landmarks.")

    temp_ls = []
    n_len = nx.number_of_nodes(g)

    if n_len < 3 and k == 1:
        return [g.nodes()[0]]

    def temp_ALT(v, end):
        return max([abs(g.node[v][str(l)] - g.node[end][str(l)]) for l in temp_ls])

    min_ls = []
    num_trials = int(n_len * trials_to_nodes_in_graph)
    #try at least two trials
    num_trials = 2 if num_trials <= 2 else num_trials

    if debug:
        print "Running", str(num_trials), "trials for random landmark selection."
    #grab the one that has the smallest avg search space
    min_search_space = n_len

    #establish a set of test nodes such that we can vet each landmark set
    test_nodes = []
    for j in range(queries_per_trial):
        samps = sample(g.nodes(), 2)
        test_nodes.append((samps[0], samps[1]))

    #for each landmark configuration
    for i in range(num_trials):
        if debug:
            print "Running random landmark selection, trial", str(i)
        #choose a set of landmarks
        temp_ls = sample(g.nodes(), k)
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
                    path, astar_size = pp.astar_path(g,s,t,temp_ALT, search_space_size=True)
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
def alp_farthest_ecc(g, subgraphs, opt = False):
    if not g:
        raise ValueError("Cannot specify graph as null.")
    elif not subgraphs:
        raise ValueError("Cannot specify subgraphs as null.")

    ls = []

    e = nx.eccentricity(g)
    i = 0
    for g_p in subgraphs.values():
        max_nodes = []
        max_ecc = 0
        #first, get the nodes with max eccentricity in the cluster
        for v in g_p:
            if e[v] > max_ecc:
                max_nodes = [v]
                max_ecc = e[v]
            elif e[v] == max_ecc:
                max_nodes.append(v)
        ls.append(choice(max_nodes))
        '''if len(max_nodes) == 1:
            #if there is only one just add it
            ls.append(max_nodes[0])
        elif not ls:
            ls.append(max_nodes[0])
        else:
            #determine which of these nodes is farthest from all landmarks
            ls_dists = []
            for l in ls:
                ls_dists.append(remove_all_but(max_nodes, pp.single_source_shortest_path_length(g, l, target_cutoffs=max_nodes)))

            #find the node furthest from current landmarks
            max_d = 0
            max_node = max_nodes[0]
            #identify node with max eccentricity
            for node in max_nodes:
                total = 0
                for dists in ls_dists:
                    total += dists[node]

                avg_dist = total/len(ls_dists)
                if avg_dist > max_d:
                    max_d = avg_dist
                    max_node = node

            #add to set
            ls.append(max_node)'''
    return ls

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
    for g_p in subgraphs.values():
        i += 1
        #print str(i)
        #print g_p
        if not ls:
            ls.append(choice(g_p))
        else:
            ls_dists = []
            for l in ls:
                if nx.get_edge_attributes(g,"weight") and not opt:
                    ls_dists.append(remove_all_but(g_p, pp.single_source_dijkstra_path_length(g, l, target_cutoffs=g_p[:])))
                else:
                    ls_dists.append(remove_all_but(g_p, pp.single_source_shortest_path_length(g, l, target_cutoffs=g_p[:])))


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
            #print ls
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
        depths = pp.single_source_dijkstra_path_length(g, v_p, target_cutoffs=cf)

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
        depths = pp.single_source_shortest_path_length(g, v_p, cutoff=cf, target_cutoffs=g.nodes())

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
        '''e = nx.eccentricity(subg)
        #get node with highest eccentricity
        if e:
            ls.append(max(e.iteritems(), key=operator.itemgetter(1))[0])
        else:
            print "Could not find eccentricity"'''
        '''max_e = 0
        max_node = 0
        #identify node with max eccentricity
        for node in g_p:
             if e[node] > max_e:
                max_e = e[node]
                max_node = node

        #add to set
        ls.append(max_node)'''
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
	
def increase_landmarks(g, old_landmark_dict, new_landmark_set):
    new_g = g.copy()
	new_landmark_dict = {}
	#print "Landmarks:", landmarks
	#One final shortest path computation amongst the landmark set
	#landmarks_dists = nx.floyd_warshall_numpy(g, alp_landmarks[:])
	#reverse_graph = g.reverse(True)
	def get_landmark_paths(l):
		commun = node_clusters[l]
		sub = cluster_nodes[commun][:]
		for a in sub:
			if DEBUG:
				print "Assigned node", str(a), "to", str(l)
			g.node[a]['landmark'] = l
		alp_label = 'ALP_'+str(l)
		nx.set_node_attributes(g, alp_label, single_source_dijkstra_path_length(g, l, target_cutoffs=sub))
		#nx.set_node_attributes(g, alp_label, single_source_dijkstra_path_length(reverse_graph, l, target_cutoffs=alp_landmarks[:]+sub))
	
	#landmark_dists is now a matrix for the landmark set
	from multiprocessing.pool import ThreadPool
	#pool = ThreadPool(processes=min(int(len(alp_landmarks)/8)+8,15))
	pool = ThreadPool(processes=25)

	#do 8 landmarks at a time
	pool.map(get_landmark_paths, alp_landmarks)
	print "Landmark distances calculated. Storing..."
	#for each new landmark in the set
	#