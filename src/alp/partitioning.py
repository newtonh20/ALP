"""Graph Voronoi partitioning for ALP preprocessing.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>
"""

from __future__ import annotations

from heapq import heappop, heappush
from typing import Any, Dict, List

import networkx as nx

DEFAULT_WEIGHT_KEY = "weight"


def voronoi_partition(
    G: nx.Graph,
    landmarks: List[Any],
    directed: bool = False,
    weight: str = DEFAULT_WEIGHT_KEY,
) -> Dict[Any, Any]:
    """
    Assign each node to its nearest landmark via a multi-source Dijkstra flood.

    Returns a mapping node -> owning landmark. Distances are implicit; this
    keeps the structure lightweight while giving enough information for
    downstream embedding.

    Args:
        G: Input graph.
        landmarks: Landmark node ids.
        directed: Preserved for API symmetry; graph orientation is inferred from G.
        weight: Edge weight attribute.

    Returns:
        Mapping from node id to its closest landmark.
    """
    if not landmarks:
        return {}

    owner: Dict[Any, Any] = {}
    best_dist: Dict[Any, float] = {}
    frontier: List[tuple[float, Any, Any]] = []  # (distance, node, source landmark)

    for src in landmarks:
        # Seed the heap with every landmark at distance 0.
        heappush(frontier, (0.0, src, src))

    while frontier:
        cur_dist, node, src = heappop(frontier)  # Always expand the closest candidate.
        if node in owner:
            continue
        owner[node] = src  # Claim ownership for this landmark.
        best_dist[node] = cur_dist  # Record best known distance.

        for nbr, edgedata in G[node].items():
            nbr_dist = cur_dist + edgedata.get(weight, 1)  # Accumulate edge weight.
            if nbr in owner:
                continue
            prev_best = best_dist.get(nbr)
            if prev_best is None or nbr_dist < prev_best:
                best_dist[nbr] = nbr_dist  # Tighten bound for neighbor.
                heappush(frontier, (nbr_dist, nbr, src))  # Reconsider with improved key.

    return owner


__all__ = ["voronoi_partition", "DEFAULT_WEIGHT_KEY"]
