/*
 * igraph_alp.c
 *
 * Description: Implements preprocessing and query-time A* heuristics for the
 * ALP (A*, Landmarks, Polygon inequalities) variant within igraph. The code
 * keeps memory bounded to one stored distance per vertex plus the landmark
 * matrix while remaining admissible for directed graphs by using a symmetric
 * auxiliary metric during preprocessing.
 *
 * Author: Newton Campbell <ncampbell@atlanticcouncil.org>
 * Date: 2025-02-07
 */

#include "igraph_alp.h"
#include <math.h>
#include <stdlib.h>

static int igraph_i_alp_choose_landmarks(const igraph_t *graph, igraph_integer_t L, igraph_vector_int_t *landmarks);
static int igraph_i_alp_build_partitions(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                         const igraph_vector_int_t *partition_labels, igraph_bool_t directed,
                                         igraph_alp_t *alp);
static int igraph_i_alp_build_landmark_dist(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                            igraph_bool_t directed, igraph_matrix_t *mat);
static igraph_real_t igraph_i_alp_dist_from_partition(const igraph_alp_t *alp, igraph_integer_t landmark_idx,
                                                      igraph_integer_t v);
static igraph_real_t igraph_i_alp_heuristic(const igraph_alp_t *alp, igraph_integer_t s,
                                            igraph_integer_t t, igraph_integer_t v);
static int igraph_i_get_shortest_path_alp_core(const igraph_t *graph, const igraph_alp_t *alp,
                                               igraph_integer_t from, igraph_integer_t to,
                                               const igraph_vector_t *weights, igraph_neimode_t mode,
                                               igraph_vector_int_t *path_out, igraph_real_t *distance_out);

/**
 * @brief Initializes an igraph_alp_t handle for preprocessing or querying.
 *
 * @param alp Target handle to initialize.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
int igraph_alp_init(igraph_alp_t *alp) {
    if (!alp) {
        return IGRAPH_EINVAL;
    }
    alp->nof_landmarks = 0;
    alp->directed = 0;
    IGRAPH_CHECK(igraph_vector_int_init(&alp->landmarks, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&alp->owner, 0));
    IGRAPH_CHECK(igraph_vector_int_init(&alp->vertex_pos, 0));
    IGRAPH_CHECK(igraph_matrix_init(&alp->landmark_dist, 0, 0));
    IGRAPH_CHECK(igraph_vector_ptr_init(&alp->partitions, 0));
    return IGRAPH_SUCCESS;
}

/**
 * @brief Releases all storage held by an ALP handle.
 *
 * Frees all partition vectors and resets the handle state. A NULL pointer is
 * accepted and ignored to simplify cleanup paths.
 *
 * @param alp Handle to release.
 */
void igraph_alp_destroy(igraph_alp_t *alp) {
    if (!alp) {
        return;
    }
    for (long i = 0; i < igraph_vector_ptr_size(&alp->partitions); ++i) {
        igraph_alp_partition_t *part = VECTOR(alp->partitions)[i];
        if (part) {
            igraph_vector_int_destroy(&part->vertices);
            igraph_vector_destroy(&part->dist);
            igraph_free(part);
        }
    }
    igraph_vector_ptr_destroy(&alp->partitions);
    igraph_matrix_destroy(&alp->landmark_dist);
    igraph_vector_int_destroy(&alp->vertex_pos);
    igraph_vector_int_destroy(&alp->owner);
    igraph_vector_int_destroy(&alp->landmarks);
}

/**
 * @brief Executes ALP preprocessing with randomly chosen landmarks.
 *
 * @param graph          Input graph.
 * @param nof_landmarks  Number of landmarks to choose uniformly at random.
 * @param directed       Whether the caller intends to run directed searches.
 * @param alp            Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
int igraph_alp_preprocess(const igraph_t *graph, igraph_integer_t nof_landmarks,
                          igraph_bool_t directed, igraph_alp_t *alp) {
    if (alp) {
        alp->nof_landmarks = nof_landmarks;
    }
    return igraph_alp_preprocess_with_partitions(graph, NULL, NULL, directed, alp);
}

/**
 * @brief Executes ALP preprocessing with optional user supplied landmarks/partitions.
 *
 * When @p partition_labels is NULL, the function builds a Voronoi-style
 * partition around the landmarks using a multi-source BFS over the auxiliary
 * undirected metric. Otherwise, it respects the user-provided labels and only
 * fills local distance arrays for vertices mapped to each landmark.
 *
 * @param graph             Input graph.
 * @param landmarks         Optional landmark vertex ids (NULL for random).
 * @param partition_labels  Optional owner label per vertex (NULL to auto-compute).
 * @param directed          Whether the input graph is directed.
 * @param alp               Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
int igraph_alp_preprocess_with_partitions(const igraph_t *graph,
                                          const igraph_vector_int_t *landmarks,
                                          const igraph_vector_int_t *partition_labels,
                                          igraph_bool_t directed,
                                          igraph_alp_t *alp) {
    int rc = 0;
    if (!graph || !alp) {
        return IGRAPH_EINVAL;
    }

    IGRAPH_CHECK(igraph_alp_init(alp));
    alp->directed = directed;

    igraph_integer_t L = landmarks ? igraph_vector_int_size(landmarks) : alp->nof_landmarks;
    if (!landmarks) {
        rc = igraph_i_alp_choose_landmarks(graph, L, &alp->landmarks);
        if (rc != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    } else {
        rc = igraph_vector_int_resize(&alp->landmarks, L);
        if (rc != IGRAPH_SUCCESS) {
            goto cleanup;
        }
        igraph_vector_int_update(&alp->landmarks, landmarks);
    }
    alp->nof_landmarks = L;

    rc = igraph_i_alp_build_partitions(graph, &alp->landmarks, partition_labels, directed, alp);
    if (rc != IGRAPH_SUCCESS) {
        goto cleanup;
    }
    rc = igraph_i_alp_build_landmark_dist(graph, &alp->landmarks, directed, &alp->landmark_dist);
    if (rc != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    return IGRAPH_SUCCESS;

cleanup:
    igraph_alp_destroy(alp);
    return rc;
}

/**
 * @brief Selects @p L distinct random landmarks from the graph.
 *
 * Uses igraph's default RNG to maintain reproducibility with other igraph
 * sampling utilities.
 *
 * @param graph      Graph to sample vertices from.
 * @param L          Number of landmarks to select.
 * @param landmarks  Output vector resized to @p L and filled with vertex ids.
 * @return IGRAPH_SUCCESS on success or an igraph error code on failure.
 */
static int igraph_i_alp_choose_landmarks(const igraph_t *graph, igraph_integer_t L, igraph_vector_int_t *landmarks) {
    igraph_integer_t n = igraph_vcount(graph);
    if (L <= 0 || L > n) {
        return IGRAPH_EINVAL;
    }
    IGRAPH_CHECK(igraph_vector_int_resize(landmarks, L));

    igraph_rng_t *rng = igraph_rng_default();
    igraph_bool_t *used = IGRAPH_CALLOC(n, igraph_bool_t);
    if (!used) {
        return IGRAPH_ENOMEM;
    }

    for (igraph_integer_t i = 0; i < L; ++i) {
        igraph_integer_t v = igraph_rng_get_integer(rng, 0, n - 1);
        while (used[v]) {
            v = igraph_rng_get_integer(rng, 0, n - 1);
        }
        VECTOR(*landmarks)[i] = v;
        used[v] = 1;
    }
    IGRAPH_FREE(used);
    return IGRAPH_SUCCESS;
}

/**
 * @brief Builds Voronoi-like partitions around landmarks and populates distances.
 *
 * When @p partition_labels is NULL, this runs a multi-source BFS over the
 * undirected view of the graph, assigning each vertex to the nearest landmark
 * and recording the distance. When labels are provided, it fills distances only
 * for the vertices mapped to each landmark while still using BFS to compute the
 * metric.
 *
 * @param graph             Graph to traverse.
 * @param landmarks         Vector of landmark vertex ids.
 * @param partition_labels  Optional owner labels per vertex (NULL for auto).
 * @param directed          Whether the input graph is directed (preprocessing uses IGRAPH_ALL).
 * @param alp               Handle to populate with partitions, owners, and vertex positions.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
static int igraph_i_alp_build_partitions(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                         const igraph_vector_int_t *partition_labels, igraph_bool_t directed,
                                         igraph_alp_t *alp) {
    (void) directed;
    igraph_integer_t n = igraph_vcount(graph);
    igraph_integer_t L = igraph_vector_int_size(landmarks);
    igraph_integer_t alloc_count = 0;
    igraph_bool_t neis_initialized = 0;
    int rc = IGRAPH_SUCCESS;

    if ((rc = igraph_vector_int_resize(&alp->owner, n))) {
        return rc;
    }
    if ((rc = igraph_vector_int_resize(&alp->vertex_pos, n))) {
        return rc;
    }

    for (igraph_integer_t i = 0; i < n; ++i) {
        VECTOR(alp->owner)[i] = -1;
        VECTOR(alp->vertex_pos)[i] = -1;
    }

    igraph_vector_t dist;
    if ((rc = igraph_vector_init(&dist, n))) {
        return rc;
    }
    igraph_vector_fill(&dist, IGRAPH_INFINITY);

    if ((rc = igraph_vector_ptr_resize(&alp->partitions, L))) {
        igraph_vector_destroy(&dist);
        return rc;
    }
    for (igraph_integer_t k = 0; k < L; ++k) {
        igraph_alp_partition_t *part = IGRAPH_CALLOC(1, igraph_alp_partition_t);
        if (!part) {
            rc = IGRAPH_ENOMEM;
            goto cleanup_error;
        }
        if ((rc = igraph_vector_int_init(&part->vertices, 0))) {
            igraph_free(part);
            goto cleanup_error;
        }
        if ((rc = igraph_vector_init(&part->dist, 0))) {
            igraph_vector_int_destroy(&part->vertices);
            igraph_free(part);
            goto cleanup_error;
        }
        VECTOR(alp->partitions)[k] = part;
        ++alloc_count;
    }

    igraph_vector_int_t neis;
    if ((rc = igraph_vector_int_init(&neis, 0))) {
        goto cleanup_error;
    }
    neis_initialized = 1;

    if (partition_labels == NULL) {
        igraph_dqueue_int_t q;
        if ((rc = igraph_dqueue_int_init(&q, 0))) {
            goto cleanup_error;
        }

        for (igraph_integer_t k = 0; k < L; ++k) {
            igraph_integer_t v = VECTOR(*landmarks)[k];
            VECTOR(dist)[v] = 0.0;
            VECTOR(alp->owner)[v] = k;
            igraph_dqueue_int_push(&q, v);
        }

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t v = igraph_dqueue_int_pop(&q);
            igraph_integer_t owner = VECTOR(alp->owner)[v];
            igraph_real_t dv = VECTOR(dist)[v];

            if ((rc = igraph_neighbors(graph, &neis, v, IGRAPH_ALL))) {
                igraph_dqueue_int_destroy(&q);
                goto cleanup_error;
            }
            long nlen = igraph_vector_int_size(&neis);
            for (long i = 0; i < nlen; ++i) {
                igraph_integer_t u = VECTOR(neis)[i];
                igraph_real_t alt = dv + 1.0;
                if (alt < VECTOR(dist)[u]) {
                    VECTOR(dist)[u] = alt;
                    VECTOR(alp->owner)[u] = owner;
                    igraph_dqueue_int_push(&q, u);
                }
            }
        }

        igraph_dqueue_int_destroy(&q);
    } else {
        if (igraph_vector_int_size(partition_labels) != n) {
            rc = IGRAPH_EINVAL;
            goto cleanup_error;
        }
        for (igraph_integer_t k = 0; k < L; ++k) {
            igraph_integer_t source = VECTOR(*landmarks)[k];
            igraph_dqueue_int_t q;
            if ((rc = igraph_dqueue_int_init(&q, 0))) {
                goto cleanup_error;
            }
            igraph_vector_fill(&dist, IGRAPH_INFINITY);
            VECTOR(dist)[source] = 0.0;
            igraph_dqueue_int_push(&q, source);

            while (!igraph_dqueue_int_empty(&q)) {
                igraph_integer_t v = igraph_dqueue_int_pop(&q);
                igraph_real_t dv = VECTOR(dist)[v];
                if ((rc = igraph_neighbors(graph, &neis, v, IGRAPH_ALL))) {
                    igraph_dqueue_int_destroy(&q);
                    goto cleanup_error;
                }
                long nlen = igraph_vector_int_size(&neis);
                for (long i = 0; i < nlen; ++i) {
                    igraph_integer_t u = VECTOR(neis)[i];
                    igraph_real_t alt = dv + 1.0;
                    if (alt < VECTOR(dist)[u]) {
                        VECTOR(dist)[u] = alt;
                        igraph_dqueue_int_push(&q, u);
                    }
                }
            }

            for (igraph_integer_t v = 0; v < n; ++v) {
                if (VECTOR(*partition_labels)[v] != k) {
                    continue;
                }
                VECTOR(alp->owner)[v] = k;
                igraph_alp_partition_t *part = VECTOR(alp->partitions)[k];
                igraph_integer_t pos = igraph_vector_size(&part->dist);
                if ((rc = igraph_vector_int_push_back(&part->vertices, v))) {
                    igraph_dqueue_int_destroy(&q);
                    goto cleanup_error;
                }
                if ((rc = igraph_vector_push_back(&part->dist, VECTOR(dist)[v]))) {
                    igraph_dqueue_int_destroy(&q);
                    goto cleanup_error;
                }
                VECTOR(alp->vertex_pos)[v] = pos;
            }

            igraph_dqueue_int_destroy(&q);
        }
    }

    if (partition_labels == NULL) {
        for (igraph_integer_t v = 0; v < n; ++v) {
            igraph_integer_t k = VECTOR(alp->owner)[v];
            if (k < 0) {
                continue;
            }
            igraph_alp_partition_t *part = VECTOR(alp->partitions)[k];
            igraph_integer_t pos = igraph_vector_size(&part->dist);
            if ((rc = igraph_vector_int_push_back(&part->vertices, v))) {
                goto cleanup_error;
            }
            if ((rc = igraph_vector_push_back(&part->dist, VECTOR(dist)[v]))) {
                goto cleanup_error;
            }
            VECTOR(alp->vertex_pos)[v] = pos;
        }
    }

    if (neis_initialized) {
        igraph_vector_int_destroy(&neis);
    }
    igraph_vector_destroy(&dist);
    return IGRAPH_SUCCESS;

cleanup_error:
    if (neis_initialized) {
        igraph_vector_int_destroy(&neis);
    }
    for (igraph_integer_t i = 0; i < alloc_count; ++i) {
        igraph_alp_partition_t *part = VECTOR(alp->partitions)[i];
        if (part) {
            igraph_vector_int_destroy(&part->vertices);
            igraph_vector_destroy(&part->dist);
            igraph_free(part);
        }
        VECTOR(alp->partitions)[i] = NULL;
    }
    igraph_vector_ptr_resize(&alp->partitions, 0);
    igraph_vector_destroy(&dist);
    return rc;
}

/**
 * @brief Computes all-pairs distances between landmarks using the aux metric.
 *
 * Runs a BFS from each landmark over the undirected neighbor set and records
 * the distance to every other landmark, storing the result in @p mat.
 *
 * @param graph     Input graph.
 * @param landmarks Landmark vertex ids.
 * @param directed  Whether the input graph is directed (ignored for preprocessing).
 * @param mat       Output L x L matrix filled with distances or IGRAPH_INFINITY.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
static int igraph_i_alp_build_landmark_dist(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                            igraph_bool_t directed, igraph_matrix_t *mat) {
    (void) directed;
    igraph_integer_t L = igraph_vector_int_size(landmarks);
    int rc = IGRAPH_SUCCESS;
    if ((rc = igraph_matrix_resize(mat, L, L))) {
        return rc;
    }
    igraph_matrix_fill(mat, IGRAPH_INFINITY);

    igraph_integer_t n = igraph_vcount(graph);
    igraph_vector_t dist;
    igraph_vector_int_t neis;
    if ((rc = igraph_vector_init(&dist, n))) {
        return rc;
    }
    if ((rc = igraph_vector_int_init(&neis, 0))) {
        igraph_vector_destroy(&dist);
        return rc;
    }

    for (igraph_integer_t i = 0; i < L; ++i) {
        igraph_integer_t source = VECTOR(*landmarks)[i];
        igraph_vector_fill(&dist, IGRAPH_INFINITY);
        igraph_dqueue_int_t q;
        if ((rc = igraph_dqueue_int_init(&q, 0))) {
            goto cleanup_error;
        }
        VECTOR(dist)[source] = 0.0;
        igraph_dqueue_int_push(&q, source);

        while (!igraph_dqueue_int_empty(&q)) {
            igraph_integer_t v = igraph_dqueue_int_pop(&q);
            igraph_real_t dv = VECTOR(dist)[v];
            if ((rc = igraph_neighbors(graph, &neis, v, IGRAPH_ALL))) {
                igraph_dqueue_int_destroy(&q);
                goto cleanup_error;
            }
            long nlen = igraph_vector_int_size(&neis);
            for (long j = 0; j < nlen; ++j) {
                igraph_integer_t u = VECTOR(neis)[j];
                igraph_real_t alt = dv + 1.0;
                if (alt < VECTOR(dist)[u]) {
                    VECTOR(dist)[u] = alt;
                    igraph_dqueue_int_push(&q, u);
                }
            }
        }

        for (igraph_integer_t j = 0; j < L; ++j) {
            igraph_integer_t target = VECTOR(*landmarks)[j];
            igraph_real_t d = VECTOR(dist)[target];
            MATRIX(*mat, i, j) = d;
        }
        igraph_dqueue_int_destroy(&q);
    }

    igraph_vector_int_destroy(&neis);
    igraph_vector_destroy(&dist);
    return IGRAPH_SUCCESS;

cleanup_error:
    igraph_vector_int_destroy(&neis);
    igraph_vector_destroy(&dist);
    return rc;
}

/**
 * @brief Returns the stored distance from a landmark to a vertex inside its partition.
 *
 * @param alp           Preprocessed ALP handle.
 * @param landmark_idx  Index of the landmark in the handle.
 * @param v             Vertex id to query.
 * @return Stored distance or IGRAPH_INFINITY when the vertex is not owned.
 */
static igraph_real_t igraph_i_alp_dist_from_partition(const igraph_alp_t *alp, igraph_integer_t landmark_idx,
                                                      igraph_integer_t v) {
    if (landmark_idx < 0 || landmark_idx >= alp->nof_landmarks) {
        return IGRAPH_INFINITY;
    }
    if (VECTOR(alp->owner)[v] != landmark_idx) {
        return IGRAPH_INFINITY;
    }
    igraph_integer_t pos = VECTOR(alp->vertex_pos)[v];
    if (pos < 0) {
        return IGRAPH_INFINITY;
    }
    igraph_alp_partition_t *part = VECTOR(alp->partitions)[landmark_idx];
    return VECTOR(part->dist)[pos];
}

/**
 * @brief Computes the ALP heuristic for a node @p v on the path from @p s to @p t.
 *
 * Applies triangle and quadrilateral inequalities using stored partition and
 * landmark distances. Any term involving an infinite operand is skipped to
 * preserve admissibility on disconnected components.
 *
 * @param alp Preprocessed ALP handle.
 * @param s   Source vertex id.
 * @param t   Target vertex id.
 * @param v   Current vertex being expanded.
 * @return Non-negative admissible heuristic value.
 */
static igraph_real_t igraph_i_alp_heuristic(const igraph_alp_t *alp, igraph_integer_t s,
                                            igraph_integer_t t, igraph_integer_t v) {
    igraph_integer_t k_s = VECTOR(alp->owner)[s];
    igraph_integer_t k_t = VECTOR(alp->owner)[t];

    igraph_real_t h = 0.0;

    igraph_real_t d_ls_t = igraph_i_alp_dist_from_partition(alp, k_s, t);
    igraph_real_t d_ls_v = igraph_i_alp_dist_from_partition(alp, k_s, v);
    igraph_real_t d_ls_s = igraph_i_alp_dist_from_partition(alp, k_s, s);
    igraph_real_t d_lt_s = igraph_i_alp_dist_from_partition(alp, k_t, s);
    igraph_real_t d_lt_v = igraph_i_alp_dist_from_partition(alp, k_t, v);
    igraph_real_t d_lt_t = igraph_i_alp_dist_from_partition(alp, k_t, t);
    igraph_real_t d_ls_lt = (k_s >= 0 && k_t >= 0) ? MATRIX(alp->landmark_dist, k_s, k_t) : IGRAPH_INFINITY;

    if (isfinite(d_ls_t) && isfinite(d_ls_v)) {
        igraph_real_t h1 = fabs(d_ls_t - d_ls_v);
        if (h1 > h) h = h1;
    }
    if (isfinite(d_lt_s) && isfinite(d_lt_v)) {
        igraph_real_t h2 = fabs(d_lt_s - d_lt_v);
        if (h2 > h) h = h2;
    }
    if (isfinite(d_ls_lt) && isfinite(d_ls_v) && isfinite(d_lt_t)) {
        igraph_real_t h3 = d_ls_lt - d_ls_v - d_lt_t;
        if (h3 < 0) h3 = 0;
        if (h3 > h) h = h3;
    }
    if (isfinite(d_ls_lt) && isfinite(d_ls_s) && isfinite(d_lt_v)) {
        igraph_real_t h4 = d_ls_lt - d_ls_s - d_lt_v;
        if (h4 < 0) h4 = 0;
        if (h4 > h) h = h4;
    }

    return h;
}

/**
 * @brief Runs A* search using the ALP heuristic.
 *
 * Allocates temporary vectors for scores and predecessors, reusing a neighbor
 * buffer for every expansion to keep allocations predictable. The search
 * respects @p mode for directed graphs while the heuristic remains admissible
 * due to its construction on the symmetrized metric.
 *
 * @param graph        Graph to search.
 * @param alp          Preprocessed ALP handle supplying the heuristic.
 * @param from         Source vertex id.
 * @param to           Target vertex id.
 * @param weights      Optional edge weights (NULL for unweighted searches).
 * @param mode         Neighbor direction (ignored on undirected graphs).
 * @param path_out     Optional vector to fill with the resulting path.
 * @param distance_out Optional distance output pointer.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
static int igraph_i_get_shortest_path_alp_core(const igraph_t *graph, const igraph_alp_t *alp,
                                               igraph_integer_t from, igraph_integer_t to,
                                               const igraph_vector_t *weights, igraph_neimode_t mode,
                                               igraph_vector_int_t *path_out, igraph_real_t *distance_out) {
    igraph_integer_t n = igraph_vcount(graph);
    igraph_vector_t g;
    igraph_vector_int_t parent;
    igraph_vector_char_t closed;
    igraph_vector_int_t neis;
    igraph_heap_t heap;
    int rc = IGRAPH_SUCCESS;

    if ((rc = igraph_vector_init(&g, n))) {
        return rc;
    }
    if ((rc = igraph_vector_int_init(&parent, n))) {
        igraph_vector_destroy(&g);
        return rc;
    }
    if ((rc = igraph_vector_char_init(&closed, n))) {
        igraph_vector_int_destroy(&parent);
        igraph_vector_destroy(&g);
        return rc;
    }
    if ((rc = igraph_vector_int_init(&neis, 0))) {
        igraph_vector_char_destroy(&closed);
        igraph_vector_int_destroy(&parent);
        igraph_vector_destroy(&g);
        return rc;
    }

    igraph_vector_fill(&g, IGRAPH_INFINITY);
    for (igraph_integer_t i = 0; i < n; ++i) {
        VECTOR(parent)[i] = -1;
        VECTOR(closed)[i] = 0;
    }

    VECTOR(g)[from] = 0.0;
    igraph_real_t h_from = igraph_i_alp_heuristic(alp, from, to, from);
    if ((rc = igraph_heap_init(&heap, n))) {
        igraph_vector_int_destroy(&neis);
        igraph_vector_char_destroy(&closed);
        igraph_vector_int_destroy(&parent);
        igraph_vector_destroy(&g);
        return rc;
    }
    igraph_heap_push_with_index(&heap, from, -(VECTOR(g)[from] + h_from));

    while (!igraph_heap_empty(&heap)) {
        igraph_integer_t v = igraph_heap_max_index(&heap);
        igraph_heap_pop(&heap);
        if (VECTOR(closed)[v]) {
            continue;
        }
        VECTOR(closed)[v] = 1;
        if (v == to) {
            break;
        }

        if ((rc = igraph_neighbors(graph, &neis, v, alp->directed ? mode : IGRAPH_ALL))) {
            goto cleanup;
        }
        long nlen = igraph_vector_int_size(&neis);
        for (long i = 0; i < nlen; ++i) {
            igraph_integer_t u = VECTOR(neis)[i];
            if (VECTOR(closed)[u]) continue;

            igraph_real_t w = 1.0;
            if (weights) {
                igraph_integer_t eid;
                int eid_rc = igraph_get_eid(graph, &eid, v, u, alp->directed ? (mode == IGRAPH_OUT ? 1 : 0) : 0, 0);
                if (eid_rc != IGRAPH_SUCCESS) {
                    continue;
                }
                w = VECTOR(*weights)[eid];
            }

            igraph_real_t alt = VECTOR(g)[v] + w;
            if (alt < VECTOR(g)[u]) {
                VECTOR(g)[u] = alt;
                VECTOR(parent)[u] = v;
                igraph_real_t h = igraph_i_alp_heuristic(alp, from, to, u);
                igraph_heap_push_with_index(&heap, u, -(alt + h));
            }
        }
    }

    if (distance_out) {
        *distance_out = VECTOR(g)[to];
    }

    if (path_out) {
        igraph_vector_int_clear(path_out);
        if (isfinite(VECTOR(g)[to])) {
            igraph_integer_t cur = to;
            while (cur != -1) {
                igraph_vector_int_push_back(path_out, cur);
                cur = VECTOR(parent)[cur];
            }
            igraph_vector_int_reverse(path_out);
        }
    }

cleanup:
    igraph_heap_destroy(&heap);
    igraph_vector_int_destroy(&neis);
    igraph_vector_char_destroy(&closed);
    igraph_vector_int_destroy(&parent);
    igraph_vector_destroy(&g);
    return rc;
}

int igraph_get_shortest_paths_alp_d(const igraph_t *graph, const igraph_alp_t *alp,
                                    igraph_vector_t *res, igraph_integer_t from,
                                    const igraph_vector_int_t *to, const igraph_vector_t *weights,
                                    igraph_neimode_t mode) {
    if (!graph || !alp || !res || !to) {
        return IGRAPH_EINVAL;
    }
    igraph_integer_t m = igraph_vector_int_size(to);
    IGRAPH_CHECK(igraph_vector_resize(res, m));
    for (igraph_integer_t i = 0; i < m; ++i) {
        igraph_integer_t target = VECTOR(*to)[i];
        igraph_real_t d = 0.0;
        IGRAPH_CHECK(igraph_i_get_shortest_path_alp_core(graph, alp, from, target,
                                                         weights, mode, NULL, &d));
        VECTOR(*res)[i] = d;
    }
    return IGRAPH_SUCCESS;
}

int igraph_get_shortest_paths_alp(const igraph_t *graph, const igraph_alp_t *alp,
                                  igraph_vector_ptr_t *res_paths, igraph_integer_t from,
                                  const igraph_vector_int_t *to, const igraph_vector_t *weights,
                                  igraph_neimode_t mode) {
    if (!graph || !alp || !res_paths || !to) {
        return IGRAPH_EINVAL;
    }
    igraph_integer_t m = igraph_vector_int_size(to);
    int rc = igraph_vector_ptr_resize(res_paths, m);
    if (rc != IGRAPH_SUCCESS) {
        return rc;
    }
    for (igraph_integer_t i = 0; i < m; ++i) {
        VECTOR(*res_paths)[i] = NULL;
    }

    for (igraph_integer_t i = 0; i < m; ++i) {
        igraph_integer_t target = VECTOR(*to)[i];
        igraph_vector_int_t *path = IGRAPH_CALLOC(1, igraph_vector_int_t);
        if (!path) {
            rc = IGRAPH_ENOMEM;
            goto cleanup_paths;
        }
        igraph_vector_int_init(path, 0);
        igraph_real_t d = 0.0;
        int core_rc = igraph_i_get_shortest_path_alp_core(graph, alp, from, target, weights, mode, path, &d);
        if (core_rc != IGRAPH_SUCCESS) {
            igraph_vector_int_destroy(path);
            IGRAPH_FREE(path);
            rc = core_rc;
            goto cleanup_paths;
        }
        VECTOR(*res_paths)[i] = path;
    }
    return IGRAPH_SUCCESS;

cleanup_paths:
    for (igraph_integer_t j = 0; j < m; ++j) {
        igraph_vector_int_t *stored = VECTOR(*res_paths)[j];
        if (stored) {
            igraph_vector_int_destroy(stored);
            IGRAPH_FREE(stored);
            VECTOR(*res_paths)[j] = NULL;
        }
    }
    return rc;
}

