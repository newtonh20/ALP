"""Landmark selection strategies for ALP/ALT.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>
"""

from __future__ import annotations

from typing import Any, Iterable, List

import networkx as nx

DEFAULT_WEIGHT_KEY = "weight"


def _best_seed(G: nx.Graph, nodes: Iterable[Any]) -> Any:
    """
    Pick a good initial landmark from a component.

    Preference goes to a high-degree node to avoid choosing isolates when
    multiple candidates exist.

    Args:
        G: Graph whose degrees are inspected.
        nodes: Node ids within a component.

    Returns:
        A chosen node id.

    Raises:
        ValueError: If no nodes are supplied.
    """
    node_list = list(nodes)  # Materialize for multiple scans.
    if not node_list:
        raise ValueError("Cannot choose a landmark from an empty node set.")
    if len(node_list) == 1:
        return node_list[0]
    # Highest-degree node offers better connectivity as an anchor.
    return max(G.degree(node_list), key=lambda item: item[1])[0]


def _farthest_node(G: nx.Graph, seeds: List[Any], weight: str) -> Any:
    """Return the node farthest from the current seed set (multi-source Dijkstra).

    Args:
        G: Input graph.
        seeds: Current landmark set.
        weight: Edge weight attribute.

    Returns:
        Node id farthest from any seed, or None when seeds is empty.
    """
    if not seeds:
        return None
    dist = nx.multi_source_dijkstra_path_length(G, seeds, weight=weight)
    farthest_node = None
    farthest_dist = -1.0
    for node, d in dist.items():
        if node in seeds:
            continue
        if d > farthest_dist:
            farthest_dist = d
            farthest_node = node
    return farthest_node


def choose_landmarks_greedy(
    G: nx.Graph,
    k: int,
    directed: bool = False,
    weight: str = DEFAULT_WEIGHT_KEY,
) -> List[Any]:
    """
    Greedy landmark selection: pick one per component, then iteratively add the farthest node.

    This mirrors the "farthest point" heuristic used in many ALT/ALP implementations,
    giving good coverage while keeping the selection inexpensive.

    Args:
        G: Input graph.
        k: Desired number of landmarks.
        directed: Respect directed components when True and graph is directed.
        weight: Edge weight attribute name.

    Returns:
        List of chosen landmarks (up to k or number of nodes).
    """
    if k <= 0 or G.number_of_nodes() == 0:
        return []

    target_count = min(k, G.number_of_nodes())
    use_directed = directed and G.is_directed()
    components = (
        nx.weakly_connected_components(G)
        if use_directed
        else nx.connected_components(G)
    )

    landmarks: List[Any] = []

    for comp in components:
        if len(landmarks) >= target_count:
            break
        comp_nodes = list(comp)
        if comp_nodes:
            landmarks.append(_best_seed(G, comp_nodes))

    if not landmarks:
        # Fallback in the (unlikely) case no component yielded a seed.
        landmarks.append(next(iter(G.nodes())))

    while len(landmarks) < target_count:
        farthest = _farthest_node(G, landmarks, weight=weight)
        if farthest is None or farthest in landmarks:
            break
        landmarks.append(farthest)

    return landmarks


__all__ = ["choose_landmarks_greedy", "DEFAULT_WEIGHT_KEY"]
