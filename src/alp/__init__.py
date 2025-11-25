"""ALP: A* with Landmarks and Polygon inequalities.

Author: Newton Campbell <ncampbell@atlanticcouncil.org>

This package implements Newton Campbell's ALP shortest path method
on top of general graphs.
"""

from .core import (
    alp_preprocess,
    alp_shortest_path,
    alp_shortest_path_length,
    ALPPreprocessedGraph,
)

__all__ = [
    "alp_preprocess",
    "alp_shortest_path",
    "alp_shortest_path_length",
    "ALPPreprocessedGraph",
]

__version__ = "0.1.0"
