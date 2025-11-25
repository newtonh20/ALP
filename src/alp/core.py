"""ALP core orchestration and public API.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import networkx as nx

from .landmarks import choose_landmarks_greedy
from .partitioning import voronoi_partition
from .embedding import compute_landmark_distances, compute_local_dists
from .astar import astar_search_alp


@dataclass
class ALPPreprocessedGraph:
    """Container for ALP preprocessing artifacts."""

    G: nx.Graph
    landmarks: List[Any]
    owner_landmark: Dict[Any, Any]
    landmark_dist: Dict[Tuple[Any, Any], float]
    local_dists: Dict[Tuple[Any, Any], float]


def alp_preprocess(
    G: nx.Graph,
    num_landmarks: int = 8,
    directed: bool = False,
    weight: str = "weight",
) -> ALPPreprocessedGraph:
    """Canonical ALP preprocessing.

    Args:
        G: Input graph.
        num_landmarks: Number of landmarks to select.
        directed: Whether to treat the graph as directed for preprocessing.
        weight: Edge attribute name for weights.

    Returns:
        ALPPreprocessedGraph containing landmarks, ownership, and distance tables.
    """

    landmarks = choose_landmarks_greedy(G, num_landmarks, directed=directed, weight=weight)
    # Assign every node to its nearest landmark (Voronoi over the graph metric).
    owner_landmark = voronoi_partition(G, landmarks, directed=directed, weight=weight)
    # Measure distances from landmarks to their owned vertices (local trees).
    local_dists = compute_local_dists(
        G,
        landmarks,
        owner_landmark,
        directed=directed,
        weight=weight,
    )
    # Compute pairwise landmark distances for triangle-inequality lower bounds.
    landmark_dist = compute_landmark_distances(
        G,
        landmarks,
        directed=directed,
        weight=weight,
    )

    return ALPPreprocessedGraph(
        G=G,
        landmarks=landmarks,
        owner_landmark=owner_landmark,
        landmark_dist=landmark_dist,
        local_dists=local_dists,
    )


def alp_shortest_path_length(
    alp_graph: ALPPreprocessedGraph,
    source: Any,
    target: Any,
    weight: Optional[str] = None,
) -> float:
    """Shortest path length using ALP A*.

    Args:
        alp_graph: Preprocessed ALP graph bundle.
        source: Source node id.
        target: Target node id.
        weight: Edge weight attribute; defaults to NetworkX behavior when None.

    Returns:
        Shortest path length from source to target.
    """

    # Delegate to the shared A* driver with landmark-informed heuristic.
    return astar_search_alp(
        alp_graph.G,
        source,
        target,
        alp_graph,
        weight=weight,
        return_path=False,
    )


def alp_shortest_path(
    alp_graph: ALPPreprocessedGraph,
    source: Any,
    target: Any,
    weight: Optional[str] = None,
) -> List[Any]:
    """Full path using ALP A*.

    Args:
        alp_graph: Preprocessed ALP graph bundle.
        source: Source node id.
        target: Target node id.
        weight: Edge weight attribute; defaults to NetworkX behavior when None.

    Returns:
        List of nodes along the shortest path.
    """

    # Return the explicit node sequence rather than only the distance.
    return astar_search_alp(
        alp_graph.G,
        source,
        target,
        alp_graph,
        weight=weight,
        return_path=True,
    )
