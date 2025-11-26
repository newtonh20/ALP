/*
 * igraph_alt.c
 *
 * Description: Implements the ALT (A*, Landmarks, Triangle inequality)
 * heuristic for igraph. This version keeps a full landmark-to-vertex distance
 * table to favor tighter heuristics while making directed searches admissible
 * by building the heuristic on the symmetrized graph metric.
 *
 * Author: OpenAI Assistant
 * Date: 2025-02-07
 */

#include "igraph_alt.h"
#include <math.h>
#include <stdlib.h>

static int igraph_i_alt_choose_landmarks(const igraph_t *graph, igraph_integer_t L, igraph_vector_int_t *landmarks);
static int igraph_i_alt_build_tables(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                     igraph_bool_t directed, igraph_matrix_t *mat);
static igraph_real_t igraph_i_alt_heuristic(const igraph_alt_t *alt, igraph_integer_t t, igraph_integer_t v);
static int igraph_i_get_shortest_path_alt_core(const igraph_t *graph, const igraph_alt_t *alt,
                                               igraph_integer_t from, igraph_integer_t to,
                                               const igraph_vector_t *weights, igraph_neimode_t mode,
                                               igraph_vector_int_t *path_out, igraph_real_t *distance_out);

/**
 * @brief Initializes an ALT handle.
 *
 * @param alt Target handle to initialize.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
int igraph_alt_init(igraph_alt_t *alt) {
    if (!alt) {
        return IGRAPH_EINVAL;
    }
    alt->nof_landmarks = 0;
    alt->directed = 0;
    IGRAPH_CHECK(igraph_vector_int_init(&alt->landmarks, 0));
    IGRAPH_CHECK(igraph_matrix_init(&alt->landmark_vertex_dist, 0, 0));
    return IGRAPH_SUCCESS;
}

/**
 * @brief Releases all memory held by an ALT handle.
 *
 * @param alt Handle to release (NULL is allowed).
 */
void igraph_alt_destroy(igraph_alt_t *alt) {
    if (!alt) {
        return;
    }
    igraph_matrix_destroy(&alt->landmark_vertex_dist);
    igraph_vector_int_destroy(&alt->landmarks);
}

/**
 * @brief Executes ALT preprocessing with randomly selected landmarks.
 *
 * @param graph          Input graph.
 * @param nof_landmarks  Number of landmarks to choose uniformly at random.
 * @param directed       Whether the caller intends to run directed searches.
 * @param alt            Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
int igraph_alt_preprocess(const igraph_t *graph, igraph_integer_t nof_landmarks,
                          igraph_bool_t directed, igraph_alt_t *alt) {
    if (alt) {
        alt->nof_landmarks = nof_landmarks;
    }
    return igraph_alt_preprocess_with_landmarks(graph, NULL, directed, alt);
}

/**
 * @brief Executes ALT preprocessing with optional user supplied landmarks.
 *
 * @param graph      Input graph.
 * @param landmarks  Optional vector of landmark vertex ids (NULL for random).
 * @param directed   Whether the input graph is directed.
 * @param alt        Output handle to populate.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
int igraph_alt_preprocess_with_landmarks(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                         igraph_bool_t directed, igraph_alt_t *alt) {
    int rc = 0;
    if (!graph || !alt) {
        return IGRAPH_EINVAL;
    }

    IGRAPH_CHECK(igraph_alt_init(alt));
    alt->directed = directed;

    igraph_integer_t L = landmarks ? igraph_vector_int_size(landmarks) : alt->nof_landmarks;
    if (!landmarks) {
        rc = igraph_i_alt_choose_landmarks(graph, L, &alt->landmarks);
        if (rc != IGRAPH_SUCCESS) {
            goto cleanup;
        }
    } else {
        rc = igraph_vector_int_resize(&alt->landmarks, L);
        if (rc != IGRAPH_SUCCESS) {
            goto cleanup;
        }
        igraph_vector_int_update(&alt->landmarks, landmarks);
    }
    alt->nof_landmarks = L;

    rc = igraph_i_alt_build_tables(graph, &alt->landmarks, directed, &alt->landmark_vertex_dist);
    if (rc != IGRAPH_SUCCESS) {
        goto cleanup;
    }

    return IGRAPH_SUCCESS;

cleanup:
    igraph_alt_destroy(alt);
    return rc;
}

/**
 * @brief Selects @p L distinct random landmarks from the graph.
 *
 * @param graph      Graph to sample from.
 * @param L          Number of landmarks to select.
 * @param landmarks  Output vector resized to @p L and filled with vertex ids.
 * @return IGRAPH_SUCCESS on success or an igraph error code on failure.
 */
static int igraph_i_alt_choose_landmarks(const igraph_t *graph, igraph_integer_t L, igraph_vector_int_t *landmarks) {
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
 * @brief Builds the ALT landmark-to-vertex distance table.
 *
 * Runs a BFS from each landmark over the undirected neighbor set to maintain an
 * admissible symmetric metric even when the caller later performs directed
 * searches.
 *
 * @param graph      Input graph.
 * @param landmarks  Vector of landmark vertex ids.
 * @param directed   Whether the caller treats the graph as directed (ignored here).
 * @param mat        Output L x |V| matrix of distances.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
static int igraph_i_alt_build_tables(const igraph_t *graph, const igraph_vector_int_t *landmarks,
                                     igraph_bool_t directed, igraph_matrix_t *mat) {
    (void) directed;
    igraph_integer_t L = igraph_vector_int_size(landmarks);
    igraph_integer_t n = igraph_vcount(graph);
    int rc = IGRAPH_SUCCESS;

    if ((rc = igraph_matrix_resize(mat, L, n))) {
        return rc;
    }
    igraph_matrix_fill(mat, IGRAPH_INFINITY);

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

        for (igraph_integer_t v = 0; v < n; ++v) {
            MATRIX(*mat, i, v) = VECTOR(dist)[v];
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
 * @brief Computes the ALT heuristic for a node @p v when targeting @p t.
 *
 * Applies the triangle inequality term |d(L, t) - d(L, v)| for every landmark
 * row in the precomputed distance table, ignoring any infinite entries.
 *
 * @param alt Preprocessed ALT handle.
 * @param t   Target vertex id.
 * @param v   Current vertex being expanded.
 * @return Non-negative admissible heuristic value.
 */
static igraph_real_t igraph_i_alt_heuristic(const igraph_alt_t *alt, igraph_integer_t t, igraph_integer_t v) {
    igraph_real_t h = 0.0;
    for (igraph_integer_t k = 0; k < alt->nof_landmarks; ++k) {
        igraph_real_t d_lt = MATRIX(alt->landmark_vertex_dist, k, t);
        igraph_real_t d_lv = MATRIX(alt->landmark_vertex_dist, k, v);
        if (isfinite(d_lt) && isfinite(d_lv)) {
            igraph_real_t h_k = fabs(d_lt - d_lv);
            if (h_k > h) {
                h = h_k;
            }
        }
    }
    return h;
}

/**
 * @brief Runs A* search using the ALT heuristic.
 *
 * Reuses a neighbor buffer across expansions and keeps the heuristic admissible
 * on directed graphs by basing it on the symmetrized metric computed in
 * preprocessing.
 *
 * @param graph        Graph to search.
 * @param alt          Preprocessed ALT handle.
 * @param from         Source vertex id.
 * @param to           Target vertex id.
 * @param weights      Optional edge weights (NULL for unweighted searches).
 * @param mode         Neighbor direction (ignored on undirected graphs).
 * @param path_out     Optional vector to fill with the resulting path.
 * @param distance_out Optional distance output pointer.
 * @return IGRAPH_SUCCESS on success or an igraph error code otherwise.
 */
static int igraph_i_get_shortest_path_alt_core(const igraph_t *graph, const igraph_alt_t *alt,
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
    igraph_real_t h_from = igraph_i_alt_heuristic(alt, to, from);
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

        if ((rc = igraph_neighbors(graph, &neis, v, alt->directed ? mode : IGRAPH_ALL))) {
            goto cleanup;
        }
        long nlen = igraph_vector_int_size(&neis);
        for (long i = 0; i < nlen; ++i) {
            igraph_integer_t u = VECTOR(neis)[i];
            if (VECTOR(closed)[u]) {
                continue;
            }

            igraph_real_t w = 1.0;
            if (weights) {
                igraph_integer_t eid;
                int eid_rc = igraph_get_eid(graph, &eid, v, u, alt->directed ? (mode == IGRAPH_OUT ? 1 : 0) : 0, 0);
                if (eid_rc != IGRAPH_SUCCESS) {
                    continue;
                }
                w = VECTOR(*weights)[eid];
            }

            igraph_real_t alt_cost = VECTOR(g)[v] + w;
            if (alt_cost < VECTOR(g)[u]) {
                VECTOR(g)[u] = alt_cost;
                VECTOR(parent)[u] = v;
                igraph_real_t h = igraph_i_alt_heuristic(alt, to, u);
                igraph_heap_push_with_index(&heap, u, -(alt_cost + h));
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

int igraph_get_shortest_paths_alt_d(const igraph_t *graph, const igraph_alt_t *alt,
                                    igraph_vector_t *res, igraph_integer_t from,
                                    const igraph_vector_int_t *to, const igraph_vector_t *weights,
                                    igraph_neimode_t mode) {
    if (!graph || !alt || !res || !to) {
        return IGRAPH_EINVAL;
    }
    igraph_integer_t m = igraph_vector_int_size(to);
    IGRAPH_CHECK(igraph_vector_resize(res, m));
    for (igraph_integer_t i = 0; i < m; ++i) {
        igraph_integer_t target = VECTOR(*to)[i];
        igraph_real_t d = 0.0;
        IGRAPH_CHECK(igraph_i_get_shortest_path_alt_core(graph, alt, from, target,
                                                         weights, mode, NULL, &d));
        VECTOR(*res)[i] = d;
    }
    return IGRAPH_SUCCESS;
}

int igraph_get_shortest_paths_alt(const igraph_t *graph, const igraph_alt_t *alt,
                                  igraph_vector_ptr_t *res_paths, igraph_integer_t from,
                                  const igraph_vector_int_t *to, const igraph_vector_t *weights,
                                  igraph_neimode_t mode) {
    if (!graph || !alt || !res_paths || !to) {
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
        int core_rc = igraph_i_get_shortest_path_alt_core(graph, alt, from, target, weights, mode, path, &d);
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

