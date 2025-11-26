from __future__ import annotations

import random

import pytest

nx = pytest.importorskip("networkx", reason="networkx dependency required for ALP integration tests")

from alp import alp_preprocess, alp_shortest_path, alp_shortest_path_length


def _connect_components(G: nx.Graph) -> None:
    """Ensure connectivity by linking components with unit-weight edges."""
    comps = list(nx.connected_components(G))
    if len(comps) <= 1:
        return
    comp_reps = [next(iter(c)) for c in comps]
    for u, v in zip(comp_reps, comp_reps[1:]):
        G.add_edge(u, v, weight=1)


def test_alp_matches_path_graph_unweighted():
    G = nx.path_graph(8)
    alpG = alp_preprocess(G, num_landmarks=3)
    for s in range(3):
        for t in range(4, 8):
            assert alp_shortest_path_length(alpG, s, t) == nx.shortest_path_length(G, s, t)


def test_alp_matches_weighted_path_graph():
    G = nx.path_graph(6)
    for u, v in G.edges():
        G[u][v]["w"] = 2 + u + v
    alpG = alp_preprocess(G, num_landmarks=2, weight="w")
    for s, t in [(0, 5), (1, 4), (2, 3)]:
        d_alp = alp_shortest_path_length(alpG, s, t, weight="w")
        d_true = nx.shortest_path_length(G, s, t, weight="w")
        assert d_alp == d_true


def test_alp_matches_grid_graph():
    G = nx.grid_2d_graph(3, 3)
    alpG = alp_preprocess(G, num_landmarks=3)
    pairs = [((0, 0), (2, 2)), ((1, 0), (2, 1)), ((0, 2), (2, 0))]
    for s, t in pairs:
        assert alp_shortest_path_length(alpG, s, t) == nx.shortest_path_length(G, s, t)


def test_alp_matches_complete_graph_weighted():
    G = nx.complete_graph(6)
    for u, v in G.edges():
        G[u][v]["w"] = (u + 1) * (v + 2)
    alpG = alp_preprocess(G, num_landmarks=4, weight="w")
    for s in range(3):
        for t in range(3, 6):
            d_alp = alp_shortest_path_length(alpG, s, t, weight="w")
            d_true = nx.shortest_path_length(G, s, t, weight="w")
            assert d_alp == d_true


def test_alp_matches_star_graph():
    G = nx.star_graph(6)
    alpG = alp_preprocess(G, num_landmarks=3)
    for leaf in range(1, 7):
        assert alp_shortest_path_length(alpG, 0, leaf) == nx.shortest_path_length(G, 0, leaf)


def test_alp_matches_directed_cycle():
    G = nx.DiGraph()
    G.add_weighted_edges_from([(i, (i + 1) % 5, 1) for i in range(5)], weight="w")
    alpG = alp_preprocess(G, num_landmarks=2, directed=True, weight="w")
    for s in range(5):
        t = (s + 3) % 5
        d_alp = alp_shortest_path_length(alpG, s, t, weight="w")
        d_true = nx.shortest_path_length(G, s, t, weight="w")
        assert d_alp == d_true


def test_alp_matches_random_connected_small():
    random.seed(0)
    G = nx.gnp_random_graph(12, 0.2, seed=0)
    _connect_components(G)
    alpG = alp_preprocess(G, num_landmarks=4)
    pairs = [(0, 5), (3, 9), (1, 11)]
    for s, t in pairs:
        assert alp_shortest_path_length(alpG, s, t) == nx.shortest_path_length(G, s, t)


def test_alp_path_sequence_matches_unique_path():
    G = nx.path_graph(7)
    alpG = alp_preprocess(G, num_landmarks=2)
    path_alp = alp_shortest_path(alpG, 1, 6)
    assert path_alp == nx.shortest_path(G, 1, 6)


def test_alp_matches_directed_weighted_graph():
    G = nx.DiGraph()
    G.add_weighted_edges_from(
        [
            ("s", "a", 2),
            ("s", "b", 5),
            ("a", "c", 2),
            ("b", "c", 1),
            ("c", "t", 3),
            ("b", "t", 10),
        ],
        weight="w",
    )
    alpG = alp_preprocess(G, num_landmarks=3, directed=True, weight="w")
    d_alp = alp_shortest_path_length(alpG, "s", "t", weight="w")
    d_true = nx.shortest_path_length(G, "s", "t", weight="w")
    assert d_alp == d_true


def test_alp_matches_complete_with_many_landmarks():
    G = nx.complete_graph(5)
    for u, v in G.edges():
        G[u][v]["w"] = u + v + 1
    alpG = alp_preprocess(G, num_landmarks=5, weight="w")
    for s in range(5):
        for t in range(5):
            if s == t:
                continue
            d_alp = alp_shortest_path_length(alpG, s, t, weight="w")
            d_true = nx.shortest_path_length(G, s, t, weight="w")
            assert d_alp == d_true
