# graph-helpers

## Functions to help manipulate dagitty and tidygraph objects

library(dagitty)
library(tidygraph)
library(igraph)
library(dplyr)
library(purrr)
library(tidyr)

# ===== Translators: DAGitty; tidygraph ==========

#' tbl_graph (directed, with a `name` node attribute) -> dagitty DAG
tidygraph_to_dagitty <- function(tg, node_var = "name") {
    stopifnot(inherits(tg, "tbl_graph"))
    nm <- as_tibble(activate(tg, nodes))[[node_var]]
    el <- tg |> activate(edges) |> as_tibble()
    lines <- c(nm, if (nrow(el) > 0) paste(nm[el$from], "->", nm[el$to]))
    dagitty(paste0("dag {\n", paste(lines, collapse = "\n"), "\n}"))
}

tg_to_dag <- tidygraph_to_dagitty

#' dagitty graph (DAG or MAG) -> tbl_graph, edge `type` preserved
dagitty_to_tidygraph <- function(g) {
    v <- names(g)
    e <- dagitty::edges(g)
    tbl_graph(
        nodes = tibble(name = v),
        edges = tibble(
            from = match(e$v, v),
            to   = match(e$w, v),
            type = e$e                       # "->", "<->", "--"
        ),
        directed = TRUE
    )
}

dag_to_tg <- dagitty_to_tidygraph

dagitty_from_edge_df <- function(edges) {
    # edges is a df with "from" and "to"
    dagitty(str_c(
        "dag{", 
        str_c(edges |> 
                  mutate(edge_string = str_c(from, "->", to)) |> 
                  pull(edge_string),
              collapse = ";\n")
        ,"}"))
}



# ======= Graph helpers ===========

gcm_nodelist <- function(g) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    # tidygraph to a tibble of nodes
    g %N>%
        as_tibble() |> 
        select(name) 
}

gcm_edgelist <- function(g){ 
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    # tidygraph to tibble of edges
    g %E>%
        as_tibble() %>%
        mutate(
            from = g %N>% pull(name) %>% .[from],
            to   = g %N>% pull(name) %>% .[to]
        ) %>%
        select(from, to)
}

gcm_possible_edges <- function(g) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    # m <- as.character(substitute(g)) |>  
    #     str_remove("^tdag_")
    gcm_nodelist(g) |> 
        rename(from = name) |> 
        cross_join(gcm_nodelist(g) |> 
                       rename(to = name))
}

gcm_all_descendants_edge_list <- function(g) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    # Get descendant edges
    expanded <- igraph::distances(g, mode = "out") < Inf
    closure_edges <- which(expanded, arr.ind = TRUE)
    
    expanded_edges <- data.frame(
        from = V(g)[closure_edges[,1]]$name,
        to   = V(g)[closure_edges[,2]]$name
    ) |> dplyr::filter(from != to)
    return(expanded_edges)
}




gcm_count_paths <- function(g) {
    # if (class(g)[[1]] != "dagitty") g <- tidygraph_to_dagitty(g)
    tibble(
        total     = nrow(as_tibble(dagitty::paths(g, limit = 1e4))),
        backdoor  = nrow(as_tibble(dagitty::paths(dagitty::backDoorGraph(g), 
                                                  limit = 1e4))),
        frontdoor = nrow(as_tibble(dagitty::paths(g, directed = TRUE, 
                                                  limit = 1e4)))
    )
}


adorn_switches_exist <- function(d, X, Y) {
    # takes a data frame and adorns a column "switched" based on
    # two other columns. Switched is TRUE for the variables X, Y if
    # somewhere else in the data frame Y, X is also present.
    d |> 
        mutate(
            switched := str_c({{X}}, "___", {{Y}}) %in% str_c({{Y}}, "___", {{X}})
        )
}

# Merging two graphs --------------

gcm_merge <- function(g1, g2) {
    if (class(g1)[[1]] == "dagitty") g1 <- dagitty_to_tidygraph(g1)
    if (class(g2)[[1]] == "dagitty") g2 <- dagitty_to_tidygraph(g2)
    nodes <- bind_rows(
        g1 %N>% as_tibble(),
        g2 %N>% as_tibble()
    ) |> distinct(name, .keep_all = TRUE)
    
    edges <- bind_rows(
        g1 %E>% as_tibble() |> 
            mutate(across(from:to, \(i) igraph::V(as.igraph(g1))$name[i])),
        g2 %E>% as_tibble() |> 
            mutate(across(from:to, \(i) igraph::V(as.igraph(g2))$name[i]))
    ) |> 
        group_by(from, to) |> 
        summarise( # Any meta-data on nodes needs to be collapsed
            across(where(is.numeric), \(x)sum(x)),
            across(where(is.character), \(c) str_c(c, collapse = ";")),
            across(where(is.factor), \(c) str_c(c, collapse = ";")),
            .groups = "drop")
    
    tbl_graph(nodes = nodes, edges = edges)
}


gcm_to_mag <- function(g, S) {
    # S is the set of nodes to keep
    if (inherits(g, "tbl_graph")) g <- tidygraph_to_dagitty(g)
    an <- setNames(lapply(S, function(v) dagitty::ancestors(g, v)), S)
    edge <- function(p) {
        x <- p[1]; y <- p[2]; rest <- setdiff(S, p)
        for (k in 0:length(rest))
            for (Z in combn(rest, k, simplify = FALSE))
                if (dagitty::dseparated(g, x, y, Z)) return(NA_character_)
        if (x %in% an[[y]])      sprintf("%s -> %s", x, y)
        else if (y %in% an[[x]]) sprintf("%s -> %s", y, x)
        else                     sprintf("%s <-> %s", x, y)
    }
    e <- na.omit(vapply(combn(S, 2, simplify = FALSE), edge, character(1)))
    dagitty::dagitty(paste0("mag { ", paste(c(S, e), collapse = " ; "), " }"))
}



# Causality helpers -------------------------

# gcm_has_directed_path <- function(g, from, to) {
#     if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
#     ig   <- as.igraph(g)
#     from_idx <- which(igraph::V(ig)$name == from)
#     to_idx   <- which(igraph::V(ig)$name == to)
#     igraph::distances(ig, v = from_idx, to = to_idx, mode = "out") |> 
#         is.finite() 
# }
# 
# gcm_get_adjustment_sets_from_node_list <- function(g, node_list) {
#     # Returns as a data frame with exposure, outcome, list of 
#     # adjustment sets (list of lists)
#     if (class(g)[[1]] != "dagitty") g <- tidygraph_to_dagitty(g)
#     df_of_tests <-
#         cross_join(
#             tibble(exposure = node_list),
#             tibble(outcome = node_list)
#         ) |> 
#         filter(exposure != outcome) |> 
#         rowwise() |> 
#         mutate(has_path = map2_lgl(
#             exposure, 
#             outcome, \(e, o) gcm_has_directed_path(g, e, o)))
#     
#     
#     df_of_tests |> 
#         mutate(
#             adjustment_sets = pmap(
#                 list(exposure, outcome),
#                 \(exp, out) g_get_adjustment_sets(g, exp, out)
#             ))
# }


# ===== Clustering DAGs ==========

gcm_cluster_graph <- function(g, .f,
                          quietly = FALSE # true to suppress messages
) {
    # Clusters a graph according to a function, .f, that remaps
    # the old node names to new ones
    # Would be good to add restrictions based on CDAG and PCDAG theory
    # This would need a "proposed" new nodes list, and then a check
    
    # might be best as an 'exception' list, and then adjust .f accordingly
    
    # could also use 'partitioning' functions that generate a partition of
    # the nodes according to some of the other algorithms
    
    # also, it might be good to track error metrics here. Options:
    # 1. adjustment set error: so list of fine-grained adjustment sets, apply .f
    # to the variables in those - are these the same as the adjustment sets in 
    # the cluster DAG?
    # 2. Same as above, but for conditional independence / d-separation
    # 3. Same, but for sets of descendants and ancestors
    # Q: Can these metrics be used for general DAG comparison?
    
    # Would also be nice to message / flag any cycles induced by the clustering
    # or bi-directional edges induced
    
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    
    new_nodes <- g %N>% 
        as_tibble() |> 
        mutate(name = .f(name)) |> 
        group_by(name) |> 
        summarise( # Any meta-data on nodes needs to be collapsed
            across(where(is.numeric), \(x)sum(x)),
            across(where(is.character), \(c) str_c(c, collapse = ";")),
            across(where(is.factor), \(c) str_c(c, collapse = ";")),
            .groups = "drop")
    
    new_edges <- g |> 
        gcm_edgelist() |> 
        mutate(across(from:to, .f)) |> 
        count(from, to) |> 
        filter(from != to) |> 
        rename(edge_w = n)
    
    tbl_graph(nodes = new_nodes, edges = new_edges)
}

#' Cluster membership induced by .f on a graph's nodes
gcm_cluster_memberships <- function(g, .f) {
    nm <- g |> activate(nodes) |> pull(name)
    split(nm, .f(nm))
}






# ======== Causal consistency measures ===========
# These are for examining causal consistency between
# any two graphs. 

# Getting L1 level conditions

gcm_implied_conditional_independencies <- function(g) {
    if ("tbl_graph" %in% class(g)) g <- tidygraph_to_dagitty(g)
    g |> 
        dagitty::impliedConditionalIndependencies() |> 
        map(\(x) tibble(X = x$X, Y = x$Y, Z = list(x$Z))) |> 
        list_rbind() |> 
        arrange(X, Y)
}

gcm_optimal_adjustment_set <- function(g, exposure, outcome) {
    # !! Code generated with Claude for this function (Fable). 
    # Seems ok on testing. 
    # O-set (Henckel, Perković & Maathuis 2022): O = pa(cn) \ forb
    # Returns a LIST mirroring dagitty::adjustmentSets():
    #   list()               -> not identifiable by adjustment
    #   list(character(0))   -> identifiable, empty adjustment set
    #   list(<chr vector>)   -> identifiable, the O-set
    
    if (!(class(g)[[1]] == "dagitty")) g <- tidygraph_to_dagitty(g)
    # cn: nodes on proper causal paths x -> y (includes y, excludes x)
    cn <- setdiff(intersect(descendants(g, exposure), 
                            ancestors(g, outcome)), 
                  exposure)
    
    if (length(cn) == 0) {
        # Null effect: O-set construction degenerates. Fall back to the
        # canonical set for exact dagitty parity (complete for any pair).
        s <- adjustmentSets(g, exposure, outcome, type = "canonical")
        if (length(s) == 0) return(list())
        return(list(sort(as.character(unname(unlist(s))))))
    }
    
    forb <- union(
        unique(unlist(lapply(cn, \(v) descendants(g, v)))),
        exposure)
    O <- setdiff(
        unique(unlist(lapply(cn, \(v) parents(g, v)))),
        forb)
    
    # Completeness (for x -> y with a causal path): O is valid iff ANY
    # valid adjustment set exists, so this doubles as the id indicator
    if (!isAdjustmentSet(g, O, exposure = exposure, outcome = outcome))
        return(list())
    
    list(sort(O))
}






## L0 consistency
# Not really a thing on the causal hierarchy, but instead just a measure
# of how many nodes (relatively) the two models share
# Here the "fingerprint" is simply the list of nodes
gcm_L0_consistency <- function(gX, gY) {
    if (inherits(gX, "tbl_graph")) gX <- tidygraph_to_dagitty(gX)
    if (inherits(gY, "tbl_graph")) gY <- tidygraph_to_dagitty(gY)
    
    n_X <- length(names(gX))
    n_Y <- length(names(gY))
    shared <- intersect(names(gX), names(gY))
    union <- union(names(gX), names(gY))
    
    jaccard <- length(shared) / length(union)
    dXY <- length(shared) / n_X
    dYX <- length(shared) / n_Y
    
    list(
        d = 1 - jaccard,
        dXY = dXY,
        dYX = dYX,
        shared_nodes = shared
    )
}






# ===================================================================
# L1 consistency, aligned with the L2 architecture
#
#   claim        one CI statement (x _||_ y | Z) over a vocabulary
#   fingerprint  the set of claims a graph makes    (gcm_ci_fingerprint)
#   scorer       directional containment            (L1_score_direction)
#   wrappers     pairwise shared-margin comparison  (gcm_L1_consistency)
#                coarsening degradation             (gcm_cluster_L1_degradation)
#
# Correctness notes vs the previous implementations:
#  * full semantic enumeration, NOT basis sets (bases are not canonical
#    across Markov-equivalent graphs)
#  * marginalization = keep claims with Z entirely inside the shared
#    vocabulary (== m-separation of the latent projection); NEVER trim
#    private nodes out of a conditioning set (that fabricates claims)
#  * both-empty fingerprints = AGREEMENT (d = 0): "no independencies"
#    is itself a shared, testable assertion (cf. the Markov-equivalent
#    triangle pair, which must score L1 distance 0)
#  * pooled per-claim aggregation, directional, mirroring L2
#
# Degradation theorems (Anand et al. 2023, C-DAG d-separation
# soundness): every quotient claim lifts to a true fine set-level
# claim, hence d_{coarse->fine} == 0 identically; nonzero values are
# bug detectors (the old VIOLATION status). d_{fine->coarse} is the
# old ci_loss with the corrected denominator (independence claims,
# not all queries).
# ===================================================================

# ---- 1. Fingerprint: CI claims over a vocabulary --------------------

gcm_ci_fingerprint <- function(g, vocab = NULL, max_cond = Inf,
                               z_as_list = FALSE) {
    # Enumerates x _||_ y | Z for unordered pairs {x, y} in vocab and
    # ALL Z subseteq vocab \ {x, y} (up to max_cond). Evaluated by
    # d-separation on the FULL graph, so with vocab = shared this is
    # exactly the marginal (latent-projected) independence model.
    # Exponential in |vocab|; cap max_cond for larger margins.
    if (inherits(g, "tbl_graph")) g <- tidygraph_to_dagitty(g)
    if (is.null(vocab)) vocab <- names(g)
    stopifnot(all(vocab %in% names(g)))
    vocab <- sort(vocab)
    
    rows <- list()
    for (p in combn(vocab, 2, simplify = FALSE)) {
        rest <- setdiff(vocab, p)
        for (k in 0:min(max_cond, length(rest))) {
            Zs <- if (k == 0) list(character(0))
            else combn(rest, k, simplify = FALSE)
            for (Z in Zs) {
                if (dseparated(g, p[1], p[2], Z)) {
                    if (z_as_list) {
                        rows[[length(rows) + 1]] <- tibble(
                            x = p[1], y = p[2],
                            # It is critical that Z is sorted here!!!
                            Z = list(sort(Z)))
                    } else {
                        rows[[length(rows) + 1]] <- tibble(
                            x = p[1], y = p[2],
                            # It is critical that Z is sorted here!!!
                            Z = paste(sort(Z), collapse = ","))
                    }
                }
            }
        }
    }
    if (length(rows) == 0) {
        if (z_as_list) {
            tibble(x = character(), y = character(), Z = list())
        } else {
            tibble(x = character(), y = character(), Z = character())
        }
    } else bind_rows(rows)
}

gcm_claim_keys <- function(fp) {
    if (nrow(fp) == 0) return(character(0))
    a  <- pmin(fp$x, fp$y)          # elementwise lexicographic min/max
    b  <- pmax(fp$x, fp$y)
    Zc <- map_chr(strsplit(fp$Z, ","), 
                  \(z) paste(sort(z, method = "radix"), collapse = ","))
    paste(a, "_||_", b, "|", Zc)
}

# ---- 2. Scorer: directional containment of claim sets ---------------

L1_score_direction <- function(fp_ref, fp_oth) {
    # d = share of the reference's claims ABSENT from the other model.
    # Unit weight per claim (each CI statement = one unit of evidence),
    # pooled -- the L1 analogue of L2_score_direction.
    #
    # Empty-reference convention (mirrors L2):
    #   ref empty, oth empty -> d = 0, one unit ("no independencies",
    #                            agreed: the saturated-model anchor)
    #   ref empty, oth not   -> d = NA, zero units (agnostic in this
    #                            direction; oth's claims are counted as
    #                            lost in the REVERSE direction)
    kr <- gcm_claim_keys(fp_ref)
    ko <- gcm_claim_keys(fp_oth)
    
    if (length(kr) == 0) {
        agree_empty <- length(ko) == 0
        return(list(
            d = if (agree_empty) 0 else NA_real_,
            n_claims = 0L,
            n_lost   = 0L,
            lost     = character(0),
            total_weight = if (agree_empty) 1 else 0
        ))
    }
    lost <- setdiff(kr, ko)
    list(
        d = length(lost) / length(kr),
        n_claims = length(kr),
        n_lost   = length(lost),
        lost     = lost,
        total_weight = length(kr)
    )
}

# ---- 3. Pairwise wrapper: two graphs, shared margin ------------------

gcm_L1_consistency <- function(gX, gY, max_cond = Inf,
                               fpX = NULL, fpY = NULL) {
    # fpX/fpY: optionally precomputed via gcm_ci_fingerprint over the
    # correct shared vocabulary (caller's responsibility). Unlike L2
    # fingerprints, L1 fingerprints are NOT reusable across different
    # comparison partners: the claim set depends on the margin, so a
    # full-vocabulary fingerprint cannot simply be row-filtered here.
    if (inherits(gX, "tbl_graph")) gX <- tidygraph_to_dagitty(gX)
    if (inherits(gY, "tbl_graph")) gY <- tidygraph_to_dagitty(gY)
    
    shared <- intersect(names(gX), names(gY))
    if (length(shared) < 2)
        stop("fewer than 2 shared nodes: no CI claims are comparable")
    
    if (is.null(fpX)) fpX <- gcm_ci_fingerprint(gX, shared, max_cond)
    if (is.null(fpY)) fpY <- gcm_ci_fingerprint(gY, shared, max_cond)

    # This captures shared CI where the co variate set Z is not shared    
    fpXfull.noZ <- gcm_ci_fingerprint(gX, max_cond = max_cond) |> 
        filter(x %in% shared, y %in% shared) |> 
        distinct(x,y)
    fpYfull.noZ <- gcm_ci_fingerprint(gY, max_cond = max_cond) |> 
        filter(x %in% shared, y %in% shared) |> 
        distinct(x,y)
    ci_from_shared_in_either <- full_join(fpXfull.noZ, fpYfull.noZ, by = c("x", "y"))
    ci_from_shared_in_both <- inner_join(fpXfull.noZ, fpYfull.noZ, by = c("x", "y"))
    
    sXY <- L1_score_direction(fpX, fpY)
    sYX <- L1_score_direction(fpY, fpX)
    
    # TODO: Is NA_real_ the correct response for n_union = 0?
    kx <- gcm_claim_keys(fpX); ky <- gcm_claim_keys(fpY)
    n_union <- length(union(kx, ky)) + nrow(ci_from_shared_in_either)
    n_intersect <- length(intersect(kx, ky)) + nrow(ci_from_shared_in_both)
    jaccard <- if (n_union == 0) NA_real_ else n_intersect / n_union
    # n_union == 0: both saturated -> identical (empty) claim sets
    
    list(
        d = 1 - jaccard,
        dXY = sXY$d, dYX = sYX$d,
        L1XY = 1 - sXY$d, L1YX = 1 - sYX$d,
        jaccard = jaccard,                    # symmetric summary
        markov_equiv_on_margin = setequal(kx, ky),
        n_claims_X = sXY$n_claims, n_claims_Y = sYX$n_claims,
        only_in_X = sXY$lost, only_in_Y = sYX$lost,
        fpX = fpX, fpY = fpY,
        shared_nodes = shared,
        ci_from_shared_in_either = ci_from_shared_in_either,
        ci_from_shared_in_both = ci_from_shared_in_both
    )
}


# ---- 1. Adjustment set for one pair --------------------------------

gcm_adj_set <- function(g, exposure, outcome,
                        adj_type = c("all" ,"optimal", "canonical")) {
    # Returns: NULL           -> not identifiable by adjustment
    #          character(0)   -> identifiable, empty set
    #          <chr vector>   -> identifiable, the set
    # NB "minimal" deliberately not offered: it returns multiple sets,
    # and unlist() would union them into a set that is not itself an
    # adjustment set.
    adj_type <- match.arg(adj_type)
    s <- switch(adj_type,
                canonical = adjustmentSets(g, exposure, outcome,
                                           type = "canonical"),
                all = adjustmentSets(g, exposure, outcome,
                                           type = "all"),
                optimal   = gcm_optimal_adjustment_set(g, exposure, outcome))
    if (length(s) == 0) return(NULL)
    sort(as.character(unname(unlist(s))))
}

# ---- 2. Per-graph fingerprint ---------------------------------------

gcm_adj_fingerprint <- function(g, pairs = NULL,
                                adj_type = c("all" ,"optimal", "canonical")) {
    # tibble: x, y, z (list of sets/NULLs), id, causal
    # `causal` = x is an ancestor of y in THIS graph
    adj_type <- match.arg(adj_type)
    if (!inherits(g, "dagitty")) g <- tidygraph_to_dagitty(g)
    if (is.null(pairs)) {
        pairs <- expand_grid(x = names(g), y = names(g)) |>
            filter(x != y)
    }
    an <- map(set_names(names(g)), \(v) ancestors(g, v))  # precompute
    
    pairs |>
        mutate(
            Z      = map2(x, y, \(x, y) gcm_adj_set(g, x, y, adj_type)),
            id     = !map_lgl(Z, is.null),
            causal = map2_lgl(x, y, \(x, y) x %in% an[[y]])
        )
}



L2_consistency_pairwise_score <- function(
        # each row input corresponds to a node-pair and their L2 statements
        id_ref, # identifiable in reference graph
        id_oth, # identifiable in comparison graph
        z_ref,  # adjustment sets in reference graph
        z_oth,  # adjustment sets in comparison graph
        use = TRUE # this node-pair should be used in the score
        ) {
    err  <- map2_int(z_ref, z_oth, \(a, b) length(setdiff(a, b)))
    size <- lengths(z_ref)                 # NULL -> 0
    unit <- pmax(size, 1)                  # max(|Z~_X|, 1)
    
    excluded <- !use | (!id_ref & id_oth)  # iota_X < iota_Y, or filtered
    
    n_neg     <- sum(!id_ref & use)             # X's negative claims
    n_refuted <- sum(!id_ref & use & id_oth)    # ...that Y refutes
    r <- if (n_neg > 0) n_refuted / n_neg else NA_real_
    
    delta <- case_when(
        excluded           ~ 0,
        id_ref & !id_oth   ~ as.numeric(unit),
        !id_ref & !id_oth  ~ 0,
        .default = as.numeric(err)         # both identify: raw omission count
    )
    w <- if_else(excluded, 0, as.numeric(unit))
    
    W <- sum(w)
    list(
        d_pair = if_else(w > 0, delta / w, NA_real_),  # per-pair, diagnostic
        w_pair = w,
        delta  = delta,                                # raw numerator terms
        d = if (W > 0) sum(delta) / W else NA_real_,
        n_effective  = sum(w > 0),
        r = r,
        n_negative_claims = n_neg,
        n_refuted         = n_refuted,
        total_weight = W
    )
}


# ---- 4. L2 consistency Orchestrator ------

gcm_L2_consistency <- function(gX, gY,
                                     adj_type = c("optimal", "canonical"),
                                     causal_pairs_only = FALSE,
                                     # empty_agreement = c("vacuous", "agree"),
                                     fpX = NULL, fpY = NULL) {
    # causal_pairs_only: restrict to pairs where x is an ancestor of y
    #   in AT LEAST ONE graph (union, not intersection: one-sided causal
    #   claims are the most diagnostic disagreements).
    # empty_agreement: how to weight pairs where both graphs identify
    #   with an empty shared-restricted set. "vacuous" = weight 0
    #   (current behaviour); "agree" = weight 1 when both sets are
    #   empty (recommended with adj_type = "optimal", see notes).
    # fpX, fpY: optionally pass precomputed full fingerprints from
    #   gcm_adj_fingerprint() to skip all dagitty calls here.
    adj_type <- match.arg(adj_type)
    
    if (!inherits(gX, "dagitty")) gX <- tidygraph_to_dagitty(gX)
    if (!inherits(gY, "dagitty")) gY <- tidygraph_to_dagitty(gY)
    
    shared <- intersect(names(gX), names(gY))
    pairs  <- expand_grid(x = shared, y = shared) |> filter(x != y)
    
    if (is.null(fpX)) fpX <- gcm_adj_fingerprint(gX, pairs, adj_type)
    else              fpX <- semi_join(fpX, pairs, by = c("x", "y"))
    if (is.null(fpY)) fpY <- gcm_adj_fingerprint(gY, pairs, adj_type)
    else              fpY <- semi_join(fpY, pairs, by = c("x", "y"))
    stopifnot(nrow(fpX) == nrow(pairs), nrow(fpY) == nrow(pairs))
    fpX <- arrange(fpX, x, y); fpY <- arrange(fpY, x, y)
    pairs <- arrange(pairs, x, y)
    
    zx_sh <- map(fpX$z, \(z) if (is.null(z)) NULL else intersect(z, shared))
    zy_sh <- map(fpY$z, \(z) if (is.null(z)) NULL else intersect(z, shared))
    use   <- if (causal_pairs_only) fpX$causal | fpY$causal else TRUE
    
    sXY <- L2_consistency_pairwise_score(fpX$id, fpY$id, zx_sh, zy_sh,
                               use)
    sYX <- L2_consistency_pairwise_score(fpY$id, fpX$id, zy_sh, zx_sh,
                               use)
    
    df.out <- pairs |>
        mutate(
            dXY = sXY$d_pair, dYX = sYX$d_pair,
            wXY = sXY$w_pair, wYX = sYX$w_pair,
            idX = fpX$id, idY = fpY$id,
            causalX = fpX$causal, causalY = fpY$causal,
            zx_sh = zx_sh, zy_sh = zy_sh
        )
    
    list(
        dXY  = sXY$d,     dYX  = sYX$d,
        # rXY: share of gX's non-identifiability claims that gY refutes
        rXY = sXY$r, rYX = sYX$r,
        # NaN/NA diagnosis: if d is NA, these say why (no usable evidence)
        n_effective_XY  = sXY$n_effective,
        n_effective_YX  = sYX$n_effective,
        total_weight_XY = sXY$total_weight,
        total_weight_YX = sYX$total_weight,
        df = df.out, shared_nodes = shared,
        settings = list(adj_type = adj_type,
                        causal_pairs_only = causal_pairs_only)
    )
}

# L1_distance_from_ci_fingerprints <- function(fpX, fpY, vocab, return_list = FALSE) {
#     # Takes a list of L1-statements (fingerprints) from 
#     # two graphs from the gcm_ci_fingerprint function
#     # and their shared nodes (the vocab)
#     
#     # Restricting to statements only between shared nodes 
#     # (on x, y, allowing any conditional nodes in Z)
#     fpXinV <- fpX |> 
#         filter(x %in% vocab, y %in% vocab) 
#     fpYinV <- fpY |> 
#         filter(x %in% vocab, y %in% vocab) 
#     
#     # All pairs already interdependent on graph
#     fpXind <- fpXinV |> filter(Z == "") |> distinct(x,y)
#     fpYind <- fpXinV |> filter(Z == "") |> distinct(x,y)
#     
#     # All possible pairs 
#     fpXpossible <- fpXinV |> distinct(x, y)
#     fpYpossible <- fpYinV |> distinct(x, y)
#     
#     trim_ci_fp_to_non_trivial_shared <- function(fp) {
#         fp |> 
#             rowwise() |> 
#             mutate(Zlist = str_split(Z, pattern = ",")) |> 
#             mutate(Z = str_c(sort(intersect(Zlist,vocab)), collapse = ",")) |> 
#             select(-Zlist) |> 
#             filter(Z != "") |> 
#             distinct() |> 
#             ungroup()
#     }
#     
#     # CI relationships only with sahred nodes
#     fpXsharedV <- 
#         fpXinV |> 
#         trim_ci_fp_to_non_trivial_shared()
#     fpYsharedV <- 
#         fpYinV |> 
#         trim_ci_fp_to_non_trivial_shared()
#     
#     # We have 3 sets of statements about pairs (x,y) 
#     # of nodes in the shared vocab
#     # 1. pairs that are already independent
#     # 2. pairs that can be conditionally independent
#     # 3. pairs that are conditionally independent on some nodes from the shared
#     # vocabulary
#     
#     independent_intersection <- inner_join(fpXind, fpYind, by = c("x", "y"))
#     independent_union <- full_join(fpXind, fpYind, by = c("x", "y"))
#     ci_possible_intersection <- inner_join(fpXpossible, fpYpossible, by = c("x", "y"))
#     ci_possible_union <- full_join(fpXpossible, fpYpossible, by = c("x", "y"))
#     ci_exact_intersection <- inner_join(fpXsharedV, fpYsharedV, by = c("x", "y", "Z"))
#     ci_exact_union <- full_join(fpXsharedV, fpYsharedV, by = c("x", "y", "Z"))
#     
#     L1_Jaccard <- (
#         nrow(independent_intersection) + 
#             nrow(ci_possible_intersection) +
#             nrow(ci_exact_intersection)) / (
#                 nrow(independent_union) +
#                     nrow(ci_possible_union) +
#                     nrow(ci_exact_union)
#             )
#     
#     L1_dist <- 1 - L1_Jaccard
#     
#     if (return_list) {
#         list(
#             d = L1_dist,
#             J = L1_Jaccard,
#             fpXind = fpXind,
#             fpXpossible = fpXpossible,
#             fpXsharedV = fpXsharedV,
#             fpYind = fpYind,
#             fpYpossible = fpYpossible,
#             fpYsharedV = fpYsharedV
#         )
#     } else {
#         return(L1_dist)
#     }
# }

L1_distance_via_latent_projection <- function(gX, gY, return_list = FALSE,
                                              max_cond = Inf) {
    S <- intersect(gcm_nodelist(gX),gcm_nodelist(gY)) |> pull(name)
    
    gXproj <- gcm_latent_projection(gX, keep_nodes = S)
    gYproj <- gcm_latent_projection(gY, keep_nodes = S)
    
    fpX <- gcm_ci_fingerprint(gXproj, max_cond = max_cond)
    fpY <- gcm_ci_fingerprint(gYproj, max_cond = max_cond)
    
    gamma_diff <- dplyr::symdiff(fpX, fpY)
    gamma_union <- dplyr::union(fpX, fpY)
    
    if (nrow(gamma_union) == 0) {
        d <- 0
    } else {
        d <- nrow(gamma_diff) / nrow(gamma_union)
    }
    
    if (!return_list) {
        return(d)
    } else {
        return(list(
            d = d, gammaX = fpX, gammaY = fpY
        ))
    }
}

L2_distance_via_latent_projection <- function(gX, gY, return_list = FALSE,
                                              adj_type = "optimal") {
    S <- intersect(gcm_nodelist(gX),gcm_nodelist(gY)) |> pull(name)
    
    gXproj <- gcm_latent_projection(gX, keep_nodes = S)
    gYproj <- gcm_latent_projection(gY, keep_nodes = S)
    
    fpX <- gcm_adj_fingerprint(gXproj, adj_type = adj_type)
    fpY <- gcm_adj_fingerprint(gYproj, adj_type = adj_type)
    
    gamma_diff <- dplyr::symdiff(fpX, fpY)
    gamma_union <- dplyr::union(fpX, fpY)
    
    if (nrow(gamma_union) == 0) {
        d <- 0
    } else {
        d <- nrow(gamma_diff) / nrow(gamma_union)
    }
    
    if (!return_list) {
        return(d)
    } else {
        return(list(
            d = d, gammaX = fpX, gammaY = fpY
        ))
    }
}

L1_distance_from_ci_fingerprints <- function(fpX, fpY, vocab, return_list = FALSE) {
    
    gamma_L1 <- function(fp) {
        fp_in <- fp |>
            filter(x %in% vocab, y %in% vocab) |>
            mutate(Zset = map(Z, \(z) setdiff(str_split(z, ",")[[1]], ""))) |>
            filter(map_lgl(Zset, \(s) all(s %in% vocab))) |>
            mutate(Zkey = map_chr(Zset, \(s)
                                  if (length(s) == 0) "<marginal>" else str_c(sort(s), collapse = ","))
            )
        
        fp_in |>
            distinct(x, y, Zkey) |>
            rename(Z = Zkey) |>
            mutate(kind = "ci")
    }
    
    gX <- gamma_L1(fpX)
    gY <- gamma_L1(fpY)
    
    inter <- inner_join(gX, gY, by = c("kind", "x", "y", "Z"))
    union <- full_join( gX, gY, by = c("kind", "x", "y", "Z"))
    
    J <- if (nrow(union) == 0) 1 else nrow(inter) / nrow(union)
    d <- 1 - J
    
    if (return_list) list(d = d, J = J, gX = gX, gY = gY) else d
}

# L2_distance_from_adj_fingerprints <- function(fpX, fpY, vocab, return_list = FALSE) {
#     # Takes a list of 
#     
#     fpXinV <- fpX |> filter(x %in% vocab, y %in% vocab)
#     fpYinV <- fpY |> filter(x %in% vocab, y %in% vocab)
#     
#     # Identifiable nodes
#     fpXid <- fpXinV |> filter(id) |> distinct(x,y)
#     fpYid <- fpYinV |> filter(id) |> distinct(x,y)
#     
#     # From identifiable, how? Sort into shared adjustmen list Z
#     adj_fp_from_shared <- function(fp) {
#         fp |> 
#             filter(id) |>
#             rowwise() |> 
#             mutate(z_in_v = list(sort(intersect(unlist(z), vocab)))) |> 
#             mutate(Z = str_c(z_in_v, collapse = ",")) |> 
#             mutate(Z = if_else(length(z) == 0, "No adjusment", Z)) |> 
#             select(x, y, Z)
#     }
#     
#     fpXsharedZ <- fpXinV |>
#         adj_fp_from_shared()
#     fpYsharedZ <- fpYinV |>
#         adj_fp_from_shared()
#     
#     identifiable_intersection <- inner_join(fpXid, fpYid, by = c("x", "y"))
#     identifiable_union <- full_join(fpXid, fpYid, by = c("x", "y"))
#     adj_set_intersection <- inner_join(fpXsharedZ, fpYsharedZ, by = c("x", "y", "Z"))
#     adj_set_union <- full_join(fpXsharedZ, fpYsharedZ, by = c("x", "y", "Z"))
#     
#     L2_Jaccard <- (
#         nrow(identifiable_intersection) +
#             nrow(adj_set_intersection)
#     ) / 
#         (nrow(identifiable_union) +
#              nrow(adj_set_union))
#     
#     L2_dist <- 1 - L2_Jaccard
#     
#     if (return_list) {
#         list(
#             d = L2_dist,
#             J = L2_Jaccard,
#             fpXid = fpXid,
#             fpXsharedZ = fpXsharedZ,
#             fpYid = fpYid,
#             fpYsharedZ = fpYsharedZ
#         )
#     } else {
#         return(L2_dist)
#     }
# }

L2_distance_from_adj_fingerprints <- function(fpX, fpY, vocab, return_list = FALSE) {
    
    gamma_L2 <- function(fp) {
        fp_in <- fp |>
            filter(x %in% vocab, y %in% vocab, id) |>
            rowwise() |>
            filter(all(unlist(z) %in% vocab)) |>
            ungroup()
        
        id_stmts <- fp_in |>
            distinct(x, y) |>
            mutate(kind = "id", Z = NA_character_)
        
        adj_stmts <- fp_in |>
            rowwise() |>
            mutate(Z = if (length(unlist(z)) == 0) "<empty>"
                   else str_c(sort(unlist(z)), collapse = ",")) |>
            ungroup() |>
            distinct(x, y, Z) |>
            mutate(kind = "adj")
        
        bind_rows(id_stmts, adj_stmts)
    }
    
    gX <- gamma_L2(fpX)
    gY <- gamma_L2(fpY)
    
    inter <- inner_join(gX, gY, by = c("kind", "x", "y", "Z"))
    union <- full_join( gX, gY, by = c("kind", "x", "y", "Z"))
    
    J <- if (nrow(union) == 0) 1 else nrow(inter) / nrow(union)
    d <- 1 - J
    
    if (return_list) list(d = d, J = J, gX = gX, gY = gY) else d
}


#### Clustering, revisted

# ===================================================================
# Coarsening degradation via the pairwise L2 machinery
#
# Principle (commuting square / approximate abstraction, Beckers,
# Eberhardt & Halpern 2019):
#
#   degradation(g, .f) = D( abstract-then-answer , answer-then-abstract )
#                      = D( fingerprint(quotient(g, .f)) ,
#                           lift(fingerprint)(g, .f) )
#
# where lift = query the FINE graph at cluster resolution (set-valued
# exposures/outcomes = cluster members) and map the resulting covariate
# sets through .f into cluster vocabulary. Both fingerprints share the
# (x, y, z, id, causal) schema, so L2_consistency_pairwise_score / the orchestrator
# consume them unchanged.
#
# Theorems inherited from the soundness results (Anand et al. 2023,
# adjustment-level version = the old recipe_sound check):
#   * coarse_id  =>  fine_id, hence
#   * r_{fine -> coarse} == 0        (refutation rate is a bug detector)
#   * the contradiction cell of d_{fine -> coarse} == old `lost_id`
# ===================================================================

# ---- Lifted fingerprint: fine graph at cluster resolution -----------

gcm_cluster_adj_fingerprint <- function(g, .f, pairs = NULL, map_sets = TRUE,
                                        adj_type = c("canonical",
                                                     "optimal")) {
    # Same schema as gcm_adj_fingerprint(), but:
    #   x, y   : cluster names
    #   z      : fine adjustment set for do(members(x)) on members(y),
    #            mapped through .f into cluster vocabulary
    #   causal : any member of x is a fine-graph ancestor of any
    #            member of y
    adj_type <- match.arg(adj_type)
    if (adj_type == "optimal") {
        # O-set completeness for set-valued exposures needs amenability
        # conditions (Henckel et al. 2022, multi-exposure case); the
        # canonical set is complete for arbitrary X, Y sets (Perkovic
        # et al. 2018), so it is the safe default here.
        warning("optimal sets with set-valued exposures are not ",
                "guaranteed complete; using canonical instead.")
        adj_type <- "canonical"
    }
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    gd  <- tidygraph_to_dagitty(g)
    mem <- gcm_cluster_memberships(g, .f)
    cl  <- names(mem)
    if (is.null(pairs))
        pairs <- expand_grid(x = cl, y = cl) |> filter(x != y)
    
    an <- map(set_names(names(gd)), \(v) ancestors(gd, v))
    
    pairs |>
        mutate(
            Z = map2(x, y, \(x, y) {
                s <- tryCatch(
                    adjustmentSets(gd,
                                   exposure = mem[[x]],
                                   outcome  = mem[[y]],
                                   type = "canonical"),
                    error = function(e) list())
                if (length(s) == 0) return(NULL)
                # in gcm_cluster_adj_fingerprint(), replace the final line of z's map2 with:
                out <- sort(as.character(unname(unlist(s))))
                if (map_sets) sort(unique(.f(out))) else out
            }),
            id = !map_lgl(Z, is.null),
            causal = map2_lgl(x, y, \(x, y)
                              any(mem[[x]] %in% unlist(an[mem[[y]]])))
        )
}


gcm_cluster_L1_degradation <- function(g, .f, noisy = TRUE, max_cond = Inf) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    
    qd  <- gcm_cluster_graph(g, .f)
    
    fp_fine   <- gcm_ci_fingerprint(
        g, 
        z_as_list = TRUE, max_cond = max_cond)      
    fp_fine_coarse_vocab <- # lifted claims
        fp_fine |>
        coarsen_vocab(.f = .f) |> 
        filter(x != y) |>
        rowwise() |> 
        filter(!(x %in% Z)) |> 
        filter(!(x %in% Z)) |> 
        ungroup() |> 
        distinct() 
    
    fp_coarse <- gcm_ci_fingerprint(
        qd, z_as_list = TRUE, max_cond = max_cond) # coarse graph claims
    
    lost_claims <- dplyr::setdiff(fp_coarse, fp_fine_coarse_vocab)
    all_claims <- dplyr::union(fp_fine_coarse_vocab, fp_coarse)
    
    degradation <- nrow(lost_claims) / nrow(all_claims)
    if (noisy) message(str_c("L1 degradation of ", degradation, " computed."))
    list(
        degradation = degradation,
        claims_in_coarse = fp_coarse,
        claims_in_fine = fp_fine_coarse_vocab
    )
}

coarsen_vocab <- function(d, .f, sep = ",") {
    
    .fl <- function(col) {
        map_chr(str_split(col, sep), \(v)
                v |> setdiff("") |> .f() |> unique() |> sort() |> str_c(collapse = sep))
    }
    
    .fll <- function(col) {
        map(col, \(v)
            v |> as.character() |> .f() |> unique() |> sort())
    }
    
    d |>
        mutate(
            across(where(is.character), .fl),
            across(where(is.list),      .fll)
        )
}

gcm_cluster_L2_degradation <- function(g, .f, noisy = TRUE,
                                       adj_type = "all") {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    
    qd  <- tidygraph_to_dagitty(gcm_cluster_graph(g, .f))
    
    fp_fine   <- gcm_adj_fingerprint(g, adj_type = adj_type) |>
    filter(causal) |> select(-causal)# lifted claims
    fp_coarse <- gcm_adj_fingerprint(qd, adj_type = adj_type) |>
    filter(causal) |> select(-causal)# quotient claims
    fp_fine_coarse_vocab <- fp_fine |> 
        coarsen_vocab(.f = .f) |> 
        filter( # exclude if x,y now merged
            x != y) |> 
        rowwise() |>
        filter( # exclude if x now a member of z
            !(x %in% Z)) |>
        filter( # exclude if y now a member of z
            !(y %in% Z)) |>
        ungroup() |>
        distinct()
    
    lost_claims <- dplyr::setdiff(fp_fine_coarse_vocab, fp_coarse)
    all_claims <- dplyr::union(fp_fine_coarse_vocab, fp_coarse)
    
    degradation <- nrow(lost_claims) / nrow(all_claims)
    if (noisy) message(str_c("L2 degradation of ", degradation, " computed."))
    list(
        degradation = degradation,
        claims_in_coarse = fp_coarse,
        claims_in_fine = fp_fine
    )
}


# Comparing two graphs ----------------

align_adjacency_matrices <- function(g1, g2, common = TRUE, nodes = NULL) {
    
    nodes1 <- g1 %N>% pull(name)
    nodes2 <- g2 %N>% pull(name)
    
    if (is.null(nodes)) {
        nodes <- if (common) intersect(nodes1, nodes2) else union(nodes1, nodes2)
    }
    
    adj1 <- as_adjacency_matrix(g1, sparse = FALSE)
    adj2 <- as_adjacency_matrix(g2, sparse = FALSE)
    
    mat1 <- matrix(0, length(nodes), 
                   length(nodes), 
                   dimnames = list(nodes, nodes))
    mat2 <- matrix(0, length(nodes), 
                   length(nodes), 
                   dimnames = list(nodes, nodes))
    
    # Restrict to nodes that are both in the matrix AND in the target node set
    shared1 <- intersect(rownames(adj1), nodes)
    shared2 <- intersect(rownames(adj2), nodes)
    
    mat1[shared1, shared1] <- adj1[shared1, shared1]
    mat2[shared2, shared2] <- adj2[shared2, shared2]
    
    list(m1 = mat1, m2 = mat2, nodes = nodes)
}






gcm_latent_projection <- function(g, 
                                  latent_nodes = NULL, 
                                  keep_nodes = NULL,
                                  noisy = FALSE) {
    
    ig <- as.igraph(g)
    all_nodes <- igraph::V(ig)$name
    
    if (is.null(keep_nodes)) {
        include_nodes <- setdiff(all_nodes, latent_nodes)
    } else {
        include_nodes <- intersect(all_nodes, keep_nodes)
    }
    latent_nodes <- setdiff(all_nodes, include_nodes)
    
    if (length(latent_nodes) == 0) {
        if (noisy) message("No latent nodes in graph, returning original graph")
        return(g)
    }
    
    edge_exists <- function(ig, from, to) {
        paste(from, to) %in% 
            paste(igraph::tail_of(ig, igraph::E(ig))$name,
                  igraph::head_of(ig, igraph::E(ig))$name)
    }
    
    ig_proj <- ig
    
    for (node in latent_nodes) {
        parents  <- igraph::neighbors(ig_proj, node, mode = "in")$name
        children <- igraph::neighbors(ig_proj, node, mode = "out")$name
        
        causal_edges <- expand_grid(from = parents, to = children) |> 
            filter(from != to,
                   !edge_exists(ig_proj, from, to))
        
        confound_edges <- expand_grid(from = children, to = children) |> 
            filter(from != to,
                   !edge_exists(ig_proj, from, to))
        
        new_edges <- bind_rows(causal_edges |> 
                                   mutate(path_type = "causal"), 
                               confound_edges |> 
                                   mutate(path_type = "confound"))
        
        if (nrow(new_edges) > 0) {
            ig_proj <- igraph::add_edges(
                ig_proj,
                as.vector(rbind(new_edges$from, new_edges$to)),
                attr = list(path_type = new_edges$path_type)
            )
        }
    }
    
    if (any(duplicated(t(apply(as.matrix(gcm_edgelist(g)), 1, sort)))))
        warning("graph contains reciprocal edges — not a DAG")
    
    igraph::delete_vertices(ig_proj, latent_nodes) |> 
        as_tbl_graph()
}

g_to_cpdag <- function(g) {
    m <- as_adjacency_matrix(g, sparse = FALSE)
    cpdag_m <- pcalg::dag2cpdag(m)
    # dag2cpdag drops dimnames, so restore from original
    dimnames(cpdag_m) <- dimnames(m)
    igraph::graph_from_adjacency_matrix(cpdag_m, mode = "directed") |>
        as_tbl_graph()
}


gcm_hamming_dist <- function(g1, g2, common = TRUE, nodes = NULL, normalise = FALSE, 
                           all_mistakes_as_one = FALSE, use_cp_dag = TRUE) {
    if (use_cp_dag) {
        g1 <- g_to_cpdag(g1)
        g2 <- g_to_cpdag(g2)
    }
    M <- align_adjacency_matrices(g1, g2, common = common, nodes = nodes)
    hd <- SID::hammingDist(M$m1, M$m2, allMistakesOne = all_mistakes_as_one)
    
    if (normalise) {
        n <- length(M$nodes)
        hd <- hd / (n * (n - 1))
    }
    
    return(hd)    
}

gcm_L2_consistency_SID <- function(g1, g2,
                                           common = TRUE, nodes = NULL) {
    if (!(class(gX)[[1]] == "dagitty")) gX <- tidygraph_to_dagitty(gX)
    if (!(class(gY)[[1]] == "dagitty")) gY <- tidygraph_to_dagitty(gY)
    
    M <- align_adjacency_matrices(g1, g2, common = common, nodes = nodes)
    sid_raw <- SID::structIntervDist(M$m1, M$m2)
    
    n <- length(M$nodes)
    max_sid <- n * (n - 1)
    
    result <- list_flatten(list(sid_raw, nodes = M$nodes))
    result$sid.normalised           <- result$sid / max_sid
    result$max_sid       <- max_sid
    
    result
}







