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


#####################################################
# graph metric functions ------------
#####################################################

gcm_node_metrics <- function(g) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    pagerank_no_outcome <- 
        g %N>%
        filter(name != "Grade") |> 
        mutate(pagerank_no_outcome = centrality_pagerank(damping = 0.85)) |> 
        as_tibble() |> 
        select(name, pagerank_no_outcome)
    
    g %N>%
        mutate(
            node_w = 1,
            node_group = code_node_group(name),
            node_type = code_node_type(name),
            measurability = code_node_measurability(name),
            degree = centrality_degree(),
            btw = centrality_betweenness(),
            btw_rank = rank(-btw, ties.method = "min"),   # rank (highest = 1)
            btw_rel_rank = (btw_rank - 1) / (n() - 1),
            btw_rel = 1 - btw_rel_rank,
            eigen = centrality_eigen(),
            pagerank = centrality_pagerank(damping = 0.5) # lowers the influence of sinks
        ) |> 
        as_tibble() |> 
        left_join(pagerank_no_outcome, by = "name") |> 
        mutate(pagerank_no_outcome = replace_na(pagerank_no_outcome, 0))
}

gcm_summary_metrics <- function(g){
    if (class(g)[[1]] == "dagitty") g <- dag_to_tg(g)
    pagerank_no_outcome <- 
        g %N>%
        filter(name != "Grade") |> 
        mutate(pagerank_no_outcome = centrality_pagerank(damping = 0.5)) |> 
        as_tibble() |> 
        select(name, pagerank_no_outcome)
    
    n.df <- g %N>%
        mutate(
            pgrnk = centrality_pagerank(damping = 0.5),
            ng = code_node_group(name),
            m = code_node_measurability(name)) |> 
        as_tibble() |> 
        left_join(pagerank_no_outcome, by = "name")
    
    # p_traced <- mean(n.df$m == "Traced")
    # p_untraced <- mean(n.df$m == "Untraced")
    # pgrnk_traced <- sum(n.df |> filter(m == "Traced") |> pull(pgrnk))
    # pgrnk_untraced <- sum(n.df |> filter(m == "Untraced") |> pull(pgrnk))
    # pgrnk_no_outcome_traced <- sum(n.df |> filter(m == "Traced") |> pull(pagerank_no_outcome))
    # pgrnk_no_outcome_untraced <- sum(n.df |> filter(m == "Untraced") |> pull(pagerank_no_outcome))
    
    p_dw <- mean(n.df$m == "DW")
    p_o <- mean(n.df$m == "O")
    p_po <- mean(n.df$m == "PO")
    p_l <- mean(n.df$m == "L")
    measure_metric <- (3 * p_dw + 2 * p_o + 1 * p_po ) / 3
    
    graph_metrics <- list(
        order     = gorder(g), # number of nodes
        size      = gsize(g), # number of edges
        density   = edge_density(g),
        diameter  = diameter(g), # longest shortest path
        mean_path = mean_distance(g),
        clustering= transitivity(g, type = "global"),
        acyclic = igraph::is_acyclic(as.igraph(g)),
        measurability = measure_metric,
        p_dw = p_dw,
        p_o = p_o,
        p_po = p_po,
        p_l = p_l
    )
    graph_metrics
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


switches_exist <- function(d, X, Y) {
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



# Causality helpers -------------------------

has_directed_path <- function(g, from, to) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    ig   <- as.igraph(g)
    from_idx <- which(igraph::V(ig)$name == from)
    to_idx   <- which(igraph::V(ig)$name == to)
    igraph::distances(ig, v = from_idx, to = to_idx, mode = "out") |> 
        is.finite() 
}

gcm_get_adjustment_sets_from_node_list <- function(g, node_list) {
    # Returns as a data frame with exposure, outcome, list of 
    # adjustment sets (list of lists)
    if (class(g)[[1]] != "dagitty") g <- tidygraph_to_dagitty(g)
    df_of_tests <-
        cross_join(
            tibble(exposure = node_list),
            tibble(outcome = node_list)
        ) |> 
        filter(exposure != outcome) |> 
        rowwise() |> 
        mutate(has_path = map2_lgl(
            exposure, 
            outcome, \(e, o) has_directed_path(g, e, o)))
    
    
    df_of_tests |> 
        mutate(
            adjustment_sets = pmap(
                list(exposure, outcome),
                \(exp, out) g_get_adjustment_sets(g, exp, out)
            ))
}


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

# ---------------------------------------
# Cluster metrics, L1, L2 consistency

# -------------------------------------------------------------------
# Gate 1: admissibility (acyclic quotient)
# -------------------------------------------------------------------

#' A coarsening is admissible as a C-DAG iff the quotient is acyclic.
#' Returns the specific obstructions: 2-cycles (the common case; these
#' are what would otherwise masquerade as bidirected edges) and any
#' larger strongly connected components.
gcm_cluster_check_admissible <- function(g, .f) {
    nm <- g |> activate(nodes) |> pull(name)
    el <- gcm_edgelist(g) |>
        mutate(across(everything(), .f)) |>
        filter(from != to) |>
        distinct(from, to)
    
    two_cycles <- el |>
        inner_join(el, by = c(from = "to", to = "from")) |>
        filter(from < to) |>
        transmute(a = from, b = to)
    
    q <- graph_from_data_frame(el, vertices = unique(.f(nm)))
    acyclic <- is_dag(q)
    og.acyclic <- is_dag(g)
    
    cyclic_clusters <- if (!acyclic) {
        comp <- components(q, mode = "strong")
        names(comp$membership)[
            comp$membership %in% which(comp$csize > 1)
        ]
    } else character(0)
    
    og.cyclic_clusters <- if (!og.acyclic) {
        og.comp <- components(g, mode = "strong")
        names(og.comp$membership)[
            og.comp$membership %in% which(og.comp$csize > 1)
        ]
    } else character(0)
    
    induced_cycles <- 
        if (og.acyclic) {
            cyclic_clusters
        } else {
            setdiff(
                .f(og.cyclic_clusters), 
                cyclic_clusters)
        }
    
    admissible <- length(induced_cycles) == 0
    
    list(
        admissible = admissible,
        og.acyclic = og.acyclic,
        acyclic      = acyclic,
        og.cyclic_clusters = og.cyclic_clusters,
        two_cycles      = two_cycles,
        cyclic_clusters = cyclic_clusters
    )
}


# -------------------------------------------------------------------
# Gate 2 + L1 cost: cluster-level CI audit
# -------------------------------------------------------------------
# For every cluster pair (A, B) and conditioning set of clusters Z:
#   coarse_indep : d-separation in the quotient
#   fine_indep   : SET-level d-separation of members in the fine DAG
#
# Theory (Anand et al. 2023, soundness of C-DAG d-separation):
#   admissible  =>  coarse_indep implies fine_indep,
# so status == "VIOLATION" should be IMPOSSIBLE; treat any occurrence
# as a bug (unit test). "lost" rows are the genuine information cost.

gcm_cluster_L1_consistency_ci_audit <- function(g, .f, max_cond = Inf) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidy(g)
    gd  <- tidy_graph_to_dag(g)
    mem <- gcm_cluster_memberships(g, .f)
    qd <- tidy_graph_to_dag(gcm_cluster_graph(g, .f))
    cl  <- names(mem)
    
    rows <- list()
    for (p in combn(cl, 2, simplify = FALSE)) {
        rest <- setdiff(cl, p)
        for (k in 0:min(max_cond, length(rest))) {
            Zs <- if (k == 0) list(character(0))
            else combn(rest, k, simplify = FALSE)
            for (Z in Zs) {
                coarse <- dseparated(qd, p[1], p[2], Z)
                # above works, but error below
                fine   <- dseparated(gd, 
                                     as.character(mem[[p[1]]]),
                                     as.character(mem[[p[2]]]),
                                     as.character(unlist(mem[Z]))
                )
                rows[[length(rows) + 1]] <- tibble(
                    A = p[1], B = p[2],
                    Z = paste(sort(Z), collapse = ","),
                    coarse_indep = coarse, fine_indep = fine)
            }
        }
    }
    bind_rows(rows) |>
        mutate(status = case_when(
            coarse_indep  & fine_indep  ~ "agree_indep",
            !coarse_indep & !fine_indep ~ "agree_dep",
            !coarse_indep & fine_indep  ~ "lost",      # information cost
            coarse_indep  & !fine_indep ~ "VIOLATION"  # impossible if admissible
        ))
}

gcm_cluster_L1_consistency_ci_loss <- function(audit) {
    stopifnot(!any(audit$status == "VIOLATION"))   # unit test on theory
    mean(audit$status == "lost")
}

# -------------------------------------------------------------------
# Cost: ancestral inflation
# -------------------------------------------------------------------
# Fine ancestry always survives coarsening; the quotient can only ADD
# ancestral relations (via concatenated between-cluster paths that have
# no fine counterpart). One-directional by construction.
# Should not function with cycles well. 

gcm_cluster_L1_consistency_ancestral_inflation <- function(g, .f) {
    gd  <- tidygraph_to_dagitty(g)
    mem <- gcm_cluster_memberships(g, .f)
    qd <- tidy_graph_to_dag(gcm_cluster_graph(g, .f))
    cl  <- names(mem)
    
    expand_grid(an = cl, de = cl) |>
        filter(an != de) |>
        rowwise() |>
        mutate(
            coarse_anc = an %in% setdiff(ancestors(qd, de), de),
            fine_anc   = length(intersect(
                mem[[an]],
                setdiff(unlist(map(mem[[de]], ~ ancestors(gd, .x))),
                        mem[[de]]))) > 0
        ) |>
        ungroup() |>
        mutate(spurious = coarse_anc & !fine_anc)
}

# wrapping all into one
gcm_cluster_L1_consistency <- function(g, .f, mx = 3) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    
    admiss <- gcm_cluster_check_admissible(g, .f)
    audit <- gcm_cluster_L1_consistency_ci_audit(g, .f, max_cond = mx)
    ci_loss <- gcm_cluster_L1_consistency_ci_loss(audit)
    anc_inflation <- gcm_cluster_L1_consistency_ancestral_inflation(g, .f)
    
    clustering_L1_consistency <- tibble(
        ci_loss = ci_loss,
        spurious_ancestors = mean(anc_inflation$spurious, na.rm = TRUE)
    )
    
    return(list(
        clustering_L1_consistency = clustering_L1_consistency,
        audit = audit,
        admissible_results = admiss,
        ancestral_inflation_results = anc_inflation
        
    ))
}

# -------------------------------------------------------------------
# L2 cost: identifiability / adjustment audit
# -------------------------------------------------------------------
# For each ordered cluster pair (X, Y):
#   fine_id    : effect of do(members(X)) on members(Y) adjustment-
#                identifiable in the fine graph
#   coarse_id  : same at cluster level
#   lost_id    : fine_id & !coarse_id  -- the coarsening destroyed
#                identifiability (e.g. mediator merged with confounder).
#                Transit-cluster theory (Tikka et al. 2023) characterizes
#                the clusterings that never do this.
#   recipe_sound : do the coarse minimal adjustment sets, expanded to
#                their members, satisfy the adjustment criterion in the
#                fine graph? (Sound by Anand et al.; kept as unit test.)

gcm_cluster_L2_consistency_adjustment_audit <- function(g, .f, type = "canonical") {
    gd  <- tidygraph_to_dagitty(g)
    mem <- gcm_cluster_memberships(g, .f)
    qd <- tidygraph_to_dagitty(gcm_cluster_graph(g, .f))
    cl  <- names(mem)
    
    expand_grid(x = cl, y = cl) |>
        filter(x != y) |>
        pmap_dfr(function(x, y) {
            coarse_sets <- tryCatch(
                adjustmentSets(qd, exposure = x, outcome = y, type = type),
                error = function(e) list())
            fine_sets <- tryCatch(
                adjustmentSets(gd, exposure = mem[[x]], outcome = mem[[y]],
                               type = type),
                error = function(e) list())
            coarse_id <- length(coarse_sets) > 0
            fine_id   <- length(fine_sets)   > 0
            sound <- if (coarse_id) {
                all(map_lgl(coarse_sets, function(Z)
                    isAdjustmentSet(gd, unlist(mem[unlist(Z)]),
                                    exposure = mem[[x]], outcome = mem[[y]])))
            } else NA
            tibble(x, y, fine_id, coarse_id,
                   lost_id      = fine_id & !coarse_id,
                   recipe_sound = sound)
        })
}

gcm_cluster_L2_consistency <- function(g, .f) {
    if (class(g)[[1]] == "dagitty") g <- dagitty_to_tidygraph(g)
    audit <- gcm_cluster_L2_consistency_adjustment_audit(g, .f)
    return(tibble(
        ErrLostId = sum(audit$lost_id, na.rm = TRUE) / sum(audit$fine_id)
    ))
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

gcm_fuzzy_L1_consistency <- function(gX, gY) {
    # if X _||_ Y | Z and Z in g1 subset of 
    
    shared <- intersect(gcm_nodelist(gX), gcm_nodelist(gY)) |> 
        pull(name)
    
    if (length(shared) < 2) {
        message()
    }
    
    imp_ci_gX <- gcm_implied_conditional_independencies(gX) |> 
        mutate(ZXinY = map(Z, \(x) x[x %in% shared])) |> 
        rowwise() |> 
        mutate(
            ZXlen = length(Z),
            ZXinYlen = length(ZXinY)) |> 
        rename(ZX = Z) |> 
        ungroup() |> 
        summarise(
            indX = any(ZXlen == 0),
            ZXinY = list(ZXinY[ZXinYlen > 0]),
            .by = c(X, Y)
        )
    
    imp_ci_gY <- gcm_implied_conditional_independencies(gY) |> 
        mutate(ZYinX = map(Z, \(x) x[x %in% shared])) |> 
        rowwise() |> 
        mutate(
            ZYlen = length(Z),
            ZYinXlen = length(ZYinX)) |> 
        rename(ZY = Z) |> 
        ungroup() |> 
        summarise(
            indY = any(ZYlen == 0),
            ZYinX = list(ZYinX[ZYinXlen > 0]),
            .by = c(X, Y)
        )
    
    imp_ci <- inner_join(
        imp_ci_gX,
        imp_ci_gY,
        by = c("X", "Y")
    ) |>
        mutate(
            .a = map(ZXinY, ~ unique(map(.x, sort))),
            .b = map(ZYinX, ~ unique(map(.x, sort))),
            Zint   = map2(
                .a, .b, 
                ~ .x[map_lgl(.x, \(s) any(map_lgl(.y, \(t) identical(s, t))))]),
            Zuni   = map2(.a, .b, ~ unique(c(.x, .y))),
            n_int  = map_int(Zint, length),
            n_uni  = map_int(Zuni, length),
            jaccard = n_int / n_uni
        ) |> 
        mutate(
            distance = case_when(
                # no overlapping conditioning sets, but agree on independence
                n_int == 0 & n_uni == 0 & indX == indY ~ 0,
                # no overlapping sets, disagree on independence 
                n_int == 0 & n_uni == 0 & indX != indY~ 0.5, 
                # overlapping sets, use jaccard distance
                TRUE ~ jaccard
            )
        ) |> 
        select(-.a, -.b) 
        
    if (nrow(imp_ci) == 0) {
        message("No shared independence relationships")
        return(
            list(
                L1soft = 0,
                L1hard = 0,
                ci_compasisons = imp_ci
            )
        )
    }
    
    # mean jaccard for soft L1, mean jaccard > 0 for hard
    L1_soft <- mean(imp_ci$distance, na.rm = T)
    L1_hard <- mean(imp_ci$distance > 0, na.rm = T)
    
    return(
        list(
        L1soft = L1_soft,
        L1hard = L1_hard,
        ci_comparisions = imp_ci)
    )
}



# adj_fingerprint 
# getting the ZX, idX etc from function below

# Needs to work better for M5, D5, D7, ...
# Many NaNs

gcm_fuzzy_L2_consistency <- function(gX, gY) {
    # Takes in two tidygraph / dagitty objects
    # Examines them for L2 causal consistency through checking if, for each pair,
    # (x,y) where x is an ancestor of y, and (x.y) are in both graphs, 
    # the adjustment sets for P(y|do(x)) meet the following criteria:
    # - If P(y|do(x)) is identifiable in g1 then it is identifiable in g2
    # - If P(y|do(x)) is identifiable in g2 then it is identifiable in g1
    # - For the canonical adjustment set Z1 = z1, z2, ..., zn in g1 
    #   and the canonical adjustment set Z2 = w1, w2, ..., wm then 
    #   (a) Z1 not in Z2 are not in V(g2)
    #   (b) Z2 not in Z1 are not in V(g1)
    #   The count of the nodes invalidating (a) and (b) is ErrXY, ErrYX
    if (!(class(gX)[[1]] == "dagitty")) gX <- tidygraph_to_dagitty(gX)
    if (!(class(gY)[[1]] == "dagitty")) gY <- tidygraph_to_dagitty(gY)
    
    shared <- intersect(names(gX), names(gY))
    
    pairs <- expand_grid(x = shared, y = shared) |> 
        filter(x != y)
    
    canAdjSet <- function(g, exposure, outcome) {
        adjustmentSets(x = g, 
                       exposure = exposure,
                       outcome = outcome,
                       type = "canonical") |> 
            unlist() |> 
            unname() |> 
            as.character() |> 
            sort()
    }
    
    df.out <- 
        pairs |> 
        rowwise() |> 
        mutate(
            canonical_ZX = pmap(
                list(x,y),
                \(x,y)
                canAdjSet(gX, x, y)
            ),
            canonical_ZY = pmap(
                list(x,y),
                \(x,y)
                canAdjSet(gY, x, y)
            )
        ) |> 
        mutate(
            idX = length(canonical_ZX) > 0,
            idY = length(canonical_ZX) > 0
        ) |> 
        mutate(
            idXY = !(idY & !idX),
            idYX = !(idX & !idY)
        ) |> 
        mutate(
            # nodes in adjustment sets in both models
            ZXshared = list(canonical_ZX[canonical_ZX %in% shared]),
            ZYshared = list(canonical_ZY[canonical_ZY %in% shared])
        ) |> 
        mutate( # Err is count of nodes in shared 
            ZXsharedMax = length(ZXshared),
            ZYsharedMax = length(ZYshared),
        ) |> 
        mutate(
            ErrXY = sum(!(ZXshared %in% ZYshared)),
            ErrYX = sum(!(ZYshared %in% ZXshared))
        ) |> 
        select(x, y, idXY, idYX, ErrXY, ErrYX, everything())
    
    idXY <-  mean(df.out$idXY)
    idYX <-  mean(df.out$idYX)
    if (sum(df.out$ZXsharedMax) > 0) {
        ErrXY <- sum(df.out$ErrXY) / sum(df.out$ZXsharedMax)
    } else {
        # mean of id agreement
        ErrXY <- mean(df.out$idXY == df.out$idYX)
    }
    if (sum(df.out$ZYsharedMax) > 0) {
        ErrYX <- sum(df.out$ErrYX) / sum(df.out$ZYsharedMax)
    } else {
        ErrYX <- mean(df.out$idYX == df.out$idXY)
    }
    
    return(list(
        idXY = idXY, idYX = idYX,
        ErrXY = ErrXY, ErrYX = ErrYX,
        df = df.out, shared_nodes = shared
    ))
    
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

gcm_latent_projection <- function(g, latent_nodes) {
    
    ig <- as.igraph(g)
    all_nodes <- igraph::V(ig)$name
    
    include_nodes <- setdiff(all_nodes, latent_nodes)
    latent_nodes <- setdiff(all_nodes, include_nodes)
    
    if (length(latent_nodes) == 0) {
        message("No latent nodes in graph, returning original graph")
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
        
        new_edges <- bind_rows(causal_edges, confound_edges)
        
        if (nrow(new_edges) > 0) {
            ig_proj <- igraph::add_edges(
                ig_proj,
                as.vector(rbind(new_edges$from, new_edges$to))
            )
        }
    }
    
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

gcm_L2_consistency <- function(gX, gY) {
    fuzzyL2 <- gcm_fuzzy_L2_consistency(gX, gY)
    sidXYL2 <- gcm_L2_consistency_SID(gX, gY)
    sidYXL2 <- gcm_L2_consistency_SID(gY, gX)
    df.out <- 
        tibble(
        sidXY = sidXYL2$sid.normalised,
        sidYX = sidYXL2$sid.normalised,
        idXY = fuzzyL2$idXY,
        idYX = fuzzyL2$idYX,
        ErrXY = fuzzyL2$ErrXY,
        ErrYX = fuzzyL2$ErrYX
    )
    return(list(
        L2_consistincy = df.out,
        fuzzy = fuzzyL2,
        sidXY = sidXYL2,
        sidYX = sidYXL2
    ))
}






