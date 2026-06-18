# graph-helpers

## Functions to help manipulate dagitty and tidygraph objects

# exposures(dag_s1)
library(igraph)
library(dagitty)

cluster_DAGitty <- function(
        dag, 
        # A function that takes a node name and returns the cluster
        phi = NULL, 
        # Alternatively, a named list of the partitions
        partition = NULL) {
    
}

g_nodelist <- function(g) {
    # tidygraph to a tibble of nodes
    g %N>%
        as_tibble() |> 
        select(name) 
}

g_edgelist <- function(g){ 
    # tidygraph to tibble of edges
    g %E>%
        as_tibble() %>%
        mutate(
            from = g %N>% pull(name) %>% .[from],
            to   = g %N>% pull(name) %>% .[to]
        ) %>%
        select(from, to)
}

g_possible_edges <- function(g) {
    # m <- as.character(substitute(g)) |>  
    #     str_remove("^tdag_")
    g_nodelist(g) |> 
        rename(from = name) |> 
        cross_join(g_nodelist(g) |> 
                       rename(to = name))
}

g_all_descendants_edge_list <- function(g) {
    # Get descendant edges
    expanded <- igraph::distances(g, mode = "out") < Inf
    closure_edges <- which(expanded, arr.ind = TRUE)
    
    expanded_edges <- data.frame(
        from = V(g)[closure_edges[,1]]$name,
        to   = V(g)[closure_edges[,2]]$name
    ) |> dplyr::filter(from != to)
    return(expanded_edges)
}

# Superseeded by g_cluster_graph
g_merge_nodes <- function(g, regex_from, string_to, no.self.loops = TRUE) {
    # merges nodes in a graph by renaming edgelist
    g_el <- g_edgelist(g)
    g_el_m <- g_el |> 
        mutate(from = str_replace(from, regex_from, string_to),
               to = str_replace(to, regex_from, string_to)) |> 
        distinct()
    tbl_graph(edges = g_el_m)
}

# same but for dags
dag_merge_nodes <- function(dag, mapping) {
    # mapping is a named list, like: c("Old.var" = "Merged.var", "Old.var2" = "Merged.var")
    g <- graph_from_data_frame(dagitty::edges(dag))
    old_names <- V(g)$name
    # Replace names using mapping if present, else keep original
    new_names <- ifelse(old_names %in% names(mapping),
                        mapping[old_names],
                        old_names)
    V(g)$name <- unname(new_names)
    new_g <- simplify(g, remove.loops = TRUE, remove.multiple = TRUE)
    
    # Convert back to dagitty
    edges_m <- igraph::as_data_frame(new_g, what="edges")
    dag_str <- paste0("dag { ", paste0(edges_m$from, " -> ", edges_m$to, collapse="; "), " }")
    g_m <- dagitty(dag_str)
    g_m
}

# Causal consistency measures

# Getting L1 level conditions
d_implied_conditional_independencies <- function(dag) {
    # Takes in a DAGitty object and exports a data frame
    # X, Y are the conditionally independent variables (order not important)
    # Z is the set of conditions under which they are independent, as a list
    dag |> 
        dagitty::impliedConditionalIndependencies() |> 
        map(\(x) tibble(X = x$X, Y = x$Y, Z = list(x$Z))) |> 
        list_rbind() |> 
        arrange(X, Y)
}

g_implied_conditional_independencies <- function(g) {
    tidy_graph_to_dag(g) |> 
        d_implied_conditional_independencies()
}


# For L2 consistency
g_get_adjustment_sets <- function(g, exposure = "S.Kn", outcome = "Grade") {
    # Compute adjustment sets
    sets <- adjustmentSets(
        tidy_graph_to_dag(g), 
        exposure = exposure, 
        outcome = outcome)
    as.list(sets)
}

has_directed_path <- function(g, from, to) {
    ig   <- as.igraph(g)
    from_idx <- which(igraph::V(ig)$name == from)
    to_idx   <- which(igraph::V(ig)$name == to)
    igraph::distances(ig, v = from_idx, to = to_idx, mode = "out") |> 
        is.finite() 
}

g_get_adjustment_sets_from_node_list <- function(g, node_list) {
    # Returns as a data frame with exposure, outcome, list of 
    # adjustment sets (list of lists)
    df_of_tests <-
        cross_join(
            tibble(exposure = node_list),
            tibble(outcome = node_list)
        ) |> 
        filter(exposure != outcome) |> 
        rowwise() |> 
        mutate(has_path = map2_lgl(exposure, outcome, \(e, o) has_directed_path(g, e, o)))

    
    
    
    df_of_tests |> 
        mutate(
            adjustment_sets = pmap(
                list(exposure, outcome),
                \(exp, out) g_get_adjustment_sets(g, exp, out)
            ))
}


# graph metric functions ------------

g_node_metrics <- function(g) {
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

g_summary_metrics <- function(g){
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


d_count_paths <- function(dag) {
    tibble(
        total     = nrow(as_tibble(dagitty::paths(dag, limit = 1e4, 
                                                  from = from, to = to))),
        backdoor  = nrow(as_tibble(dagitty::paths(dagitty::backDoorGraph(dag), 
                                                  from = from, to = to,
                                                  limit = 1e4))),
        frontdoor = nrow(as_tibble(dagitty::paths(dag, directed = TRUE, 
                                                  from = from, to = to,
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

g_merge <- function(g1, g2) {
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

g_latent_projection <- function(g, latent_nodes) {
    
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


g_hamming_dist <- function(g1, g2, common = TRUE, nodes = NULL, normalise = FALSE, 
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

g_structural_intervention_dist <- function(g1, g2,
                                           common = TRUE, nodes = NULL) {
    M <- align_adjacency_matrices(g1, g2, common = common, nodes = nodes)
    sid_raw <- SID::structIntervDist(M$m1, M$m2)
    
    n <- length(M$nodes)
    max_sid <- n * (n - 1)
    
    result <- list_flatten(list(sid_raw, nodes = M$nodes))
    result$sid.normalised           <- result$sid / max_sid
    result$max_sid       <- max_sid
    
    result
}
    
    


# Graph clustering / coarsening ---------------

g_cluster_graph <- function(g, .f,
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
        g_edgelist() |> 
        mutate(across(from:to, .f)) |> 
        count(from, to) |> 
        filter(from != to) |> 
        rename(edge_w = n)
    
    tbl_graph(nodes = new_nodes, edges = new_edges)
}
