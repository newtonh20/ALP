import networkx as nx
from alp import alp_preprocess, alp_shortest_path

def main():
    G = nx.cycle_graph(8)
    alp_graph = alp_preprocess(G, num_landmarks=2)
    path = alp_shortest_path(alp_graph, 0, 4)
    print("Shortest path from 0 to 4:", path)

if __name__ == "__main__":
    main()
