#ifndef IGRAPH_ALP_H
#define IGRAPH_ALP_H

/*
 * igraph_alp.h
 *
 * Description: Public API for the ALP (A*, Landmarks, Polygon inequalities)
 * heuristic used to accelerate shortest-path queries in igraph while keeping
 * memory usage at one distance per vertex plus an L^2 landmark matrix.
 *
 * Author: Newton Campbell <ncampbell@atlanticcouncil.org>
 * Date: 2025-02-07
 */

#include <igraph/igraph.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Stores the per-landmark partition information for ALP.
 *
 * Each partition owns a disjoint subset of vertices. For every owned vertex we
 * keep exactly one distance value from the owning landmark, preserving ALP's
 * low-memory objective compared to ALT's L x |V| footprint.
 */
typedef struct igraph_alp_partition_t {
    igraph_vector_int_t vertices; /**< Vertices owned by the landmark. */
    igraph_vector_t dist;         /**< Distances from the landmark to the owned vertices. */
} igraph_alp_partition_t;

/**
 * @brief Opaque handle for ALP preprocessing results.
 *
 * The handle captures the chosen landmarks, the Voronoi-style partition each
 * vertex belongs to, a per-vertex index into the owning partition, and an
 * L x L matrix of landmark-to-landmark distances used by the polygon
 * inequalities in the heuristic.
 */
typedef struct igraph_alp_t {
    igraph_integer_t nof_landmarks; /**< Number of landmarks. */
    igraph_bool_t directed;         /**< Whether the input graph was directed. */

    igraph_vector_int_t landmarks;  /**< Landmark vertex ids. */
    igraph_vector_int_t owner;      /**< owner[v] gives the owning landmark index. */
    igraph_vector_int_t vertex_pos; /**< Position of v inside its partition vector. */

    igraph_matrix_t landmark_dist;  /**< Landmark-to-landmark distances in auxiliary metric. */
    igraph_vector_ptr_t partitions; /**< Vector of igraph_alp_partition_t*. */
} igraph_alp_t;

/**
 * @brief Initializes an ALP handle.
 *
 * Initializes all vectors and matrices to empty so the handle can be used with
 * preprocessing.
 *
 * @param alp Pointer to the handle to initialize.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 *
 * @note Must be paired with igraph_alp_destroy() to release allocated memory.
 */
IGRAPH_EXPORT int igraph_alp_init(igraph_alp_t *alp);

/**
 * @brief Releases all memory held by an ALP handle.
 *
 * @param alp Pointer to the handle to destroy. NULL is allowed and ignored.
 */
IGRAPH_EXPORT void igraph_alp_destroy(igraph_alp_t *alp);

/**
 * @brief Runs ALP preprocessing with automatically selected landmarks.
 *
 * @param graph          Input graph.
 * @param nof_landmarks  Number of landmarks to select randomly.
 * @param directed       Whether the input graph is directed. Preprocessing
 *                       always uses an undirected metric for admissibility.
 * @param alp            Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code on failure.
 */
IGRAPH_EXPORT int igraph_alp_preprocess(
    const igraph_t *graph,
    igraph_integer_t nof_landmarks,
    igraph_bool_t directed,
    igraph_alp_t *alp);

/**
 * @brief Runs ALP preprocessing with user supplied landmarks and partitions.
 *
 * @param graph             Input graph.
 * @param landmarks         Optional landmark vertex ids (NULL for auto).
 * @param partition_labels  Optional owner labels per vertex sized |V| (NULL for auto).
 * @param directed          Whether the input graph is directed.
 * @param alp               Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code on failure.
 */
IGRAPH_EXPORT int igraph_alp_preprocess_with_partitions(
    const igraph_t *graph,
    const igraph_vector_int_t *landmarks,
    const igraph_vector_int_t *partition_labels,
    igraph_bool_t directed,
    igraph_alp_t *alp);

/**
 * @brief Computes shortest path distances using ALP-guided A*.
 *
 * The heuristic is derived from the symmetric auxiliary metric built during
 * preprocessing, keeping admissibility even when the input graph is directed.
 *
 * @param graph    Graph to search.
 * @param alp      Preprocessed ALP handle.
 * @param res      Output vector of distances (length equals size of @p to).
 * @param from     Source vertex id.
 * @param to       Vector of target vertex ids.
 * @param weights  Optional edge weights (NULL for unweighted).
 * @param mode     Neighbor mode to use during search on directed graphs.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_get_shortest_paths_alp_d(
    const igraph_t *graph,
    const igraph_alp_t *alp,
    igraph_vector_t *res,
    igraph_integer_t from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_neimode_t mode);

/**
 * @brief Computes shortest paths using ALP-guided A* and returns the routes.
 *
 * @param graph      Graph to search.
 * @param alp        Preprocessed ALP handle.
 * @param res_paths  Output vector of path vectors (length equals size of @p to).
 * @param from       Source vertex id.
 * @param to         Vector of target vertex ids.
 * @param weights    Optional edge weights (NULL for unweighted).
 * @param mode       Neighbor mode to use during search on directed graphs.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_get_shortest_paths_alp(
    const igraph_t *graph,
    const igraph_alp_t *alp,
    igraph_vector_ptr_t *res_paths,
    igraph_integer_t from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_neimode_t mode);

#ifdef __cplusplus
}
#endif

#endif /* IGRAPH_ALP_H */
