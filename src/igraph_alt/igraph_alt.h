#ifndef IGRAPH_ALT_H
#define IGRAPH_ALT_H

/*
 * igraph_alt.h
 *
 * Description: Public API for the ALT (A*, Landmarks, Triangle inequality)
 * heuristic integration in igraph. Unlike ALP, ALT stores a full L x |V|
 * distance table to deliver tighter bounds at the cost of additional memory.
 *
 * Author: OpenAI Assistant
 * Date: 2025-02-07
 */

#include <igraph/igraph.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Opaque handle for ALT preprocessing results.
 *
 * Stores the chosen landmarks, a flag indicating whether the input graph was
 * directed, and a landmark-to-vertex distance matrix sized L x |V| that is
 * consumed by the ALT heuristic.
 */
typedef struct igraph_alt_t {
    igraph_integer_t nof_landmarks; /**< Number of landmarks. */
    igraph_bool_t directed;         /**< Whether the input graph was directed. */

    igraph_vector_int_t landmarks;  /**< Landmark vertex ids. */
    igraph_matrix_t landmark_vertex_dist; /**< L x |V| matrix of distances. */
} igraph_alt_t;

/**
 * @brief Initializes an ALT handle.
 *
 * @param alt Pointer to the handle to initialize.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_alt_init(igraph_alt_t *alt);

/**
 * @brief Releases all memory held by an ALT handle.
 *
 * @param alt Pointer to the handle to destroy (NULL is ignored).
 */
IGRAPH_EXPORT void igraph_alt_destroy(igraph_alt_t *alt);

/**
 * @brief Runs ALT preprocessing with automatically selected landmarks.
 *
 * @param graph          Input graph.
 * @param nof_landmarks  Number of landmarks to select randomly.
 * @param directed       Whether the caller plans to run directed searches.
 * @param alt            Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_alt_preprocess(
    const igraph_t *graph,
    igraph_integer_t nof_landmarks,
    igraph_bool_t directed,
    igraph_alt_t *alt);

/**
 * @brief Runs ALT preprocessing with user supplied landmarks.
 *
 * @param graph      Input graph.
 * @param landmarks  Optional vector of landmark vertex ids (NULL for random).
 * @param directed   Whether the input graph is directed.
 * @param alt        Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_alt_preprocess_with_landmarks(
    const igraph_t *graph,
    const igraph_vector_int_t *landmarks,
    igraph_bool_t directed,
    igraph_alt_t *alt);

/**
 * @brief Computes shortest path distances using ALT-guided A*.
 *
 * @param graph    Graph to search.
 * @param alt      Preprocessed ALT handle.
 * @param res      Output vector of distances (length equals size of @p to).
 * @param from     Source vertex id.
 * @param to       Vector of target vertex ids.
 * @param weights  Optional edge weights (NULL for unweighted).
 * @param mode     Neighbor mode to use during search on directed graphs.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_get_shortest_paths_alt_d(
    const igraph_t *graph,
    const igraph_alt_t *alt,
    igraph_vector_t *res,
    igraph_integer_t from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_neimode_t mode);

/**
 * @brief Computes shortest paths using ALT-guided A* and returns the routes.
 *
 * @param graph      Graph to search.
 * @param alt        Preprocessed ALT handle.
 * @param res_paths  Output vector of path vectors (length equals size of @p to).
 * @param from       Source vertex id.
 * @param to         Vector of target vertex ids.
 * @param weights    Optional edge weights (NULL for unweighted).
 * @param mode       Neighbor mode to use during search on directed graphs.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
IGRAPH_EXPORT int igraph_get_shortest_paths_alt(
    const igraph_t *graph,
    const igraph_alt_t *alt,
    igraph_vector_ptr_t *res_paths,
    igraph_integer_t from,
    const igraph_vector_int_t *to,
    const igraph_vector_t *weights,
    igraph_neimode_t mode);

#ifdef __cplusplus
}
#endif

#endif /* IGRAPH_ALT_H */
