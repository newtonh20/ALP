"""Heuristic construction for ALP A* search.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, Optional, Tuple

DEFAULT_WEIGHT_KEY = "weight"


def _lookup_local_distance(
    local_dists: Dict[Tuple[Any, Any], float],
    landmark: Any,
    node: Any,
) -> Optional[float]:
    """Return the precomputed distance from landmark to node (or None).

    Args:
        local_dists: Mapping of (landmark, node) -> distance.
        landmark: Landmark node id.
        node: Query node id.

    Returns:
        Distance if present, else None.
    """
    if landmark is None or node is None:
        return None
    if landmark == node:
        return 0.0
    return local_dists.get((landmark, node))


def make_alp_heuristic(
    target: Any,
    landmarks: Iterable[Any],
    local_dists: Dict[Tuple[Any, Any], float],
) -> Callable[[Any], float]:
    """
    Build an ALP/ALT-style heuristic function for A*.

    The heuristic uses the polygon (triangle) inequality: for any landmark L
    with known distances to both the current node and the target, a lower bound
    is |d(L, target) - d(L, node)|. Taking the max over landmarks preserves
    admissibility and consistency.

    Args:
        target: Target node id.
        landmarks: Landmark node ids.
        local_dists: Mapping (landmark, node) -> distance.

    Returns:
        Callable that computes an admissible/consistent heuristic h(node).
    """
    target_landmark_dists: Dict[Any, float] = {}
    for landmark in landmarks:
        # Capture only landmarks that can reach the target.
        dist = _lookup_local_distance(local_dists, landmark, target)
        if dist is not None:
            target_landmark_dists[landmark] = dist

    if not target_landmark_dists:
        return lambda _node, _target=None: 0.0

    node_dist_cache: Dict[Tuple[Any, Any], Optional[float]] = {}

    def _distance_to_landmark(node: Any, landmark: Any) -> Optional[float]:
        key = (node, landmark)
        if key in node_dist_cache:
            return node_dist_cache[key]
        # Lazy lookup caches distances to avoid repeated dictionary hits.
        dist = _lookup_local_distance(local_dists, landmark, node)
        node_dist_cache[key] = dist
        return dist

    def heuristic(node: Any, _target: Any = None) -> float:
        if node == target:
            return 0.0

        best = 0.0
        for landmark, target_dist in target_landmark_dists.items():
            node_dist = _distance_to_landmark(node, landmark)
            if node_dist is None:
                continue
            # Triangle inequality lower bound for this landmark.
            diff = abs(target_dist - node_dist)
            if diff > best:
                best = diff
        return best

    return heuristic


__all__ = ["make_alp_heuristic", "DEFAULT_WEIGHT_KEY"]
