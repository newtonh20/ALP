"""A* search driver wired to the ALP heuristic.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>
"""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING

import networkx as nx

from .heuristic import make_alp_heuristic

if TYPE_CHECKING:
    from .core import ALPPreprocessedGraph


def astar_search_alp(
    G: nx.Graph,
    source: Any,
    target: Any,
    alp_graph: "ALPPreprocessedGraph",
    weight: Optional[str] = None,
    return_path: bool = True,
):
    """
    Run A* using the ALP landmark heuristic.

    Args:
        G: Graph to search.
        source: Source node.
        target: Target node.
        alp_graph: Preprocessed ALP data (landmarks, local distances, etc.).
        weight: Edge attribute for weights (defaults to NetworkX default when None).
        return_path: When True, return the node sequence; otherwise, only path length.

    Returns:
        Path (list of nodes) or path length depending on return_path.
    """
    heuristic = make_alp_heuristic(
        target=target,
        landmarks=alp_graph.landmarks,
        local_dists=alp_graph.local_dists,
    )

    if return_path:
        # Use NetworkX's A* to recover the node sequence.
        return nx.astar_path(
            G,
            source,
            target,
            heuristic=heuristic,
            weight=weight,
        )

    # Otherwise only the path length is required.
    return nx.astar_path_length(
        G,
        source,
        target,
        heuristic=heuristic,
        weight=weight,
    )


__all__ = ["astar_search_alp"]
