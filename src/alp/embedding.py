"""Distance embedding utilities for ALP preprocessing.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Tuple

import networkx as nx

DEFAULT_WEIGHT_KEY = "weight"


def compute_local_dists(
    G: nx.Graph,
    landmarks: Iterable[Any],
    owner_landmark: Dict[Any, Any],
    directed: bool = False,
    weight: str = DEFAULT_WEIGHT_KEY,
) -> Dict[Tuple[Any, Any], float]:
    """
    Compute landmark -> node distances.

    Distances are recorded for nodes owned by the landmark (per the partition).
    If ownership is empty, distances for every node are stored, which is useful
    for undivided graphs.

    Args:
        G: Input graph.
        landmarks: Landmark node ids.
        owner_landmark: Mapping node -> owning landmark.
        directed: Preserved for API symmetry; graph orientation is inferred from G.
        weight: Edge weight attribute name.

    Returns:
        Mapping (landmark, node) -> distance.
    """
    local_dists: Dict[Tuple[Any, Any], float] = {}
    if not landmarks:
        return local_dists

    for landmark in landmarks:
        # Dijkstra from each landmark to obtain distances within its region.
        lengths = nx.single_source_dijkstra_path_length(G, landmark, weight=weight)
        for node, dist in lengths.items():
            if owner_landmark and owner_landmark.get(node) not in (None, landmark):
                continue
            # Persist only distances that belong to this landmark's Voronoi cell.
            local_dists[(landmark, node)] = dist

    return local_dists


def compute_landmark_distances(
    G: nx.Graph,
    landmarks: List[Any],
    directed: bool = False,
    weight: str = DEFAULT_WEIGHT_KEY,
) -> Dict[Tuple[Any, Any], float]:
    """
    Compute pairwise distances between landmarks.

    For undirected graphs, the result is mirrored; for directed graphs, the
    mapping only includes reachable ordered pairs.

    Args:
        G: Input graph.
        landmarks: Landmark node ids.
        directed: Whether to respect directionality; only used when G is directed.
        weight: Edge weight attribute.

    Returns:
        Mapping (landmark_i, landmark_j) -> shortest path distance.
    """
    pair_distances: Dict[Tuple[Any, Any], float] = {}
    if not landmarks:
        return pair_distances

    use_directed = directed and G.is_directed()

    if use_directed:
        for src in landmarks:
            # Directed graph: compute ordered distances from each landmark.
            lengths = nx.single_source_dijkstra_path_length(G, src, weight=weight)
            for tgt in landmarks:
                if src == tgt:
                    continue
                if tgt in lengths:
                    pair_distances[(src, tgt)] = lengths[tgt]
    else:
        for idx, src in enumerate(landmarks):
            # Undirected graph: compute once per pair and mirror the result.
            lengths = nx.single_source_dijkstra_path_length(G, src, weight=weight)
            for tgt in landmarks[idx + 1 :]:
                if tgt in lengths:
                    dist = lengths[tgt]
                    pair_distances[(src, tgt)] = dist
                    pair_distances[(tgt, src)] = dist

    return pair_distances


__all__ = ["compute_local_dists", "compute_landmark_distances", "DEFAULT_WEIGHT_KEY"]
