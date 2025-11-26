import math
import random

import pytest

nx = pytest.importorskip("networkx", reason="networkx dependency required for ALP heuristic validation")

from alp import alp_preprocess, alp_shortest_path_length
from alp.heuristic import make_alp_heuristic


def _build_weighted_connected_graph(n: int, p: float) -> nx.Graph:
    rng = random.Random(42)
    G = nx.gnp_random_graph(n, p, seed=rng)
    # ensure connectivity by joining components with weight 1 edges
    comps = list(nx.connected_components(G))
    if len(comps) > 1:
        reps = [next(iter(c)) for c in comps]
        for u, v in zip(reps, reps[1:]):
            G.add_edge(u, v)
    for u, v in G.edges():
        G[u][v]["w"] = 1 + rng.random() * 4
    return G


def test_heuristic_is_admissible_undirected():
    G = _build_weighted_connected_graph(12, 0.2)
    alpG = alp_preprocess(G, num_landmarks=4, weight="w")

    # Pick a target and compare heuristic bounds against true distances.
    target = 0
    heuristic = make_alp_heuristic(target=target, landmarks=alpG.landmarks, local_dists=alpG.local_dists)
    true_dist = nx.single_source_dijkstra_path_length(G, target, weight="w")

    for node in G.nodes():
        h_val = heuristic(node)
        # Heuristic must never overestimate the real directed distance.
        assert h_val <= true_dist[node] + 1e-9


def test_heuristic_is_admissible_directed_symmetrized():
    rng = random.Random(7)
    G = nx.DiGraph()
    # Build a bidirectional ring with asymmetric weights to ensure strong connectivity.
    size = 10
    for i in range(size):
        w_forward = 1 + rng.random() * 3
        w_backward = 1 + rng.random() * 3
        G.add_edge(i, (i + 1) % size, w=w_forward)
        G.add_edge((i + 1) % size, i, w=w_backward)

    alpG = alp_preprocess(G, num_landmarks=3, directed=True, weight="w")
    target = 4
    heuristic = make_alp_heuristic(target=target, landmarks=alpG.landmarks, local_dists=alpG.local_dists)

    # Directed distances to target via reversed graph traversal.
    true_dist = nx.single_source_dijkstra_path_length(G.reverse(copy=False), target, weight="w")

    for node in G.nodes():
        h_val = heuristic(node)
        assert h_val <= true_dist[node] + 1e-9


def test_preprocessing_reusable_for_multiple_queries():
    G = nx.grid_2d_graph(3, 4)
    alpG = alp_preprocess(G, num_landmarks=3)

    pairs = [((0, 0), (2, 3)), ((1, 0), (2, 2)), ((0, 2), (2, 0))]
    for source, target in pairs:
        d_alp = alp_shortest_path_length(alpG, source, target)
        d_true = nx.shortest_path_length(G, source, target)
        assert d_alp == d_true


def test_owner_local_distance_matches_dijkstra():
    G = _build_weighted_connected_graph(8, 0.3)
    alpG = alp_preprocess(G, num_landmarks=3, weight="w")

    # Verify stored local distances align with shortest-path lengths from each owner landmark.
    for node, owner in alpG.owner_landmark.items():
        if (owner, node) not in alpG.local_dists:
            raise AssertionError(f"Missing local distance for node {node} owned by {owner}")
        stored = alpG.local_dists[(owner, node)]
        true_length = nx.shortest_path_length(G, owner, node, weight="w")
        assert math.isclose(stored, true_length, rel_tol=1e-9, abs_tol=1e-9)
