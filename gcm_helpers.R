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

gcm_nodelist <- function(g, as_tibble = TRUE) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    # tidygraph to a tibble of nodes
    if (as_tibble) {
        g %N>%
            as_tibble() |> 
            select(name) 
    } else {
        g %N>%
            as_tibble() |> 
            pull(name)
    }
}

gcm_edgelist <- function(g){ 
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    # tidygraph to tibble of edges
    g %E>%
        as_tibble() %>%
        mutate(
            from = g %N>% pull(name) %>% .[from],
            to   = g %N>% pull(name) %>% .[to]
        ) %>%
        select(from, to)
}

gcm_descendants_fingerprint <- function(g, vocab = NULL,
                                        full_table = FALSE) {
    if (inherits(g, "tbl_graph")) g <- tidygraph_to_dagitty(g)
    if (is.null(vocab)) {
        S <- gcm_nodelist(g, as_tibble = FALSE)
    } else {
        S <- vocab
    }
    
    df.full <- 
        expand_grid(x = S, y = S) |> 
        filter(x != y) |> 
        mutate(
            y_ancestors = map(
                y, 
                \(v)
                dagitty::ancestors(g, v)
            )
        ) |> 
        mutate(causal_path = pmap_lgl(list(x,y_ancestors),
                                  \(xx, yy)
                                  xx %in% yy
                                  )) 
    if (full_table) {
        return(df.full)
    } else {
        df.full |> 
            select(-y_ancestors) |> 
            distinct()
    }
    
}

gcm_possible_edges <- function(g) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    # m <- as.character(substitute(g)) |>  
    #     str_remove("^tdag_")
    gcm_nodelist(g) |> 
        rename(from = name) |> 
        cross_join(gcm_nodelist(g) |> 
                       rename(to = name))
}

gcm_all_descendants_edge_list <- function(g) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
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
    if (inherits(g, "tbl_graph")) g <- dagitty_to_tidygraph(g)
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

gcm_merge <- function(g1, g2, sep = ";",
                      expand_edge_types = TRUE,
                      edge_weight_var = NULL) {
    if (inherits(g1, "dagitty")) g1 <- dagitty_to_tidygraph(g1)
    if (inherits(g2, "dagitty")) g2 <- dagitty_to_tidygraph(g2)
    
    nodes <- bind_rows(g1 %N>% as_tibble(), g2 %N>% as_tibble()) |>
        distinct(name, .keep_all = TRUE)
    
    named_edge_df <- function(g) {
        nm <- igraph::V(g)$name
        ne <- 
            g %E>% 
            as_tibble() |> 
            mutate(across(from:to, \(i) nm[i]))
        
        if (!("type" %in% names(ne))) {
            ne <- ne |> 
                mutate(type = "->")
        }
        
        if (!is.null(edge_weight_var)) {
            if (!(edge_weight_var %in% names(ne))) {
                ne <- ne |> 
                    mutate({{edge_weight_var}} := 1)
                message(paste(
                    "Imputed edge weights, check if", edge_weight_var, 
                    "is the correct variable for wedge weight."))
            }
        }
        return(ne)
    }
    
    first_eq_last <- function(x) str_sub(x, 1, 1) == str_sub(x, -1, -1)
    
    # Need to impute 'type' and 'edge_w' if needed.
    edges <- bind_rows(
        named_edge_df(g1), 
        named_edge_df(g2)) |>
        mutate( 
            # This deals with edges like X -- Y in one graph and Y -- X in
            # the other. We want these described the same. 
            .sym = first_eq_last(type), # checking if edge type is symmetric
            .a   = if_else(.sym & to < from, to, from),
            .b   = if_else(.sym & to < from, from, to),
            from = .a, to = .b
        ) |>
        select(-.sym, -.a, -.b) |>
        group_by(from, to, type) |>
        summarise(
            across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
            across(where(is.character) | where(is.factor),
                   \(x) str_flatten(unique(as.character(x)), collapse = sep,
                                    na.rm = TRUE)),
            .groups = "drop"
        )
    
    tbl_graph(nodes = nodes, edges = edges)
}

gcm_split_types <- function(g, sep = ";") {
    tbl_graph(
        nodes = g %N>% as_tibble(),
        edges = g %E>% as_tibble() |>
            separate_longer_delim(type, delim = sep) |>
            distinct()
    )
}


# Quick plot of a GCM
#' Quick plot of a GCM / DAG / PDAG / PAG with mixed edge types
#'
#' Endpoint marks are parsed from the first and last character of `type`,
#' so both the PAG notation ("o->", "o-o", "-->", "<-o") and the dagitty
#' notation ("@->", "@-@") work, alongside "->", "<->", "--".
#'
#'   first char:  "<" arrowhead   "-" tail   "o"/"@" circle
#'   last  char:  ">" arrowhead   "-" tail   "o"/"@" circle
#'
#' Circles are drawn as open points just inside the node, since ggraph has
#' no circle arrow type. Their offset is in data units, so `mark_offset`
#' may need tuning once per figure size.
#'
#' @param edge_weight_var name of a numeric edge column, as a string. NULL
#'   (default) draws every edge at `flat_alpha`. Otherwise the column is
#'   mapped to edge alpha over `alpha_range`, so weak edges fade rather
#'   than disappear.

gcm_qplot <- function(g,
                      layout          = "sugiyama",
                      node_size       = 6,
                      label_size      = 3,
                      arrow_mm        = 2.5,
                      arc             = 0.25,
                      mark_size       = 1.9,
                      mark_offset     = NULL,
                      edge_weight_var = NULL,
                      alpha_range     = c(0.3, 0.95),
                      flat_alpha      = 0.9,
                      title           = NULL) {
    
    g <- tidygraph::as_tbl_graph(g)
    node_df <- tibble::as_tibble(g, active = "nodes")
    
    if (!"name" %in% names(node_df)) {
        g <- g |> activate(nodes) |> mutate(name = as.character(dplyr::row_number()))
        node_df <- tibble::as_tibble(g, active = "nodes")
    }
    if (!"type" %in% names(tibble::as_tibble(g, active = "edges")))
        g <- g |> activate(edges) |> mutate(type = "->")
    
    ## ---- resolve the weight column ------------------------------------------
    e_tbl <- tibble::as_tibble(g, active = "edges")
    use_w <- FALSE
    
    if (!is.null(edge_weight_var)) {
        if (!edge_weight_var %in% names(e_tbl)) {
            warning("edge_weight_var '", edge_weight_var,
                    "' not found in edge data; using flat alpha")
        } else if (!is.numeric(e_tbl[[edge_weight_var]])) {
            warning("edge_weight_var '", edge_weight_var,
                    "' is not numeric; using flat alpha")
        } else {
            use_w <- TRUE
        }
    }
    
    g <- g |> activate(edges) |> mutate(
        .w = if (use_w) .data[[edge_weight_var]] else flat_alpha
    )
    
    if (use_w && anyNA(tibble::as_tibble(g, active = "edges")$.w)) {
        warning("NA values in '", edge_weight_var, "' set to the observed minimum")
        g <- g |> activate(edges) |>
            mutate(.w = tidyr::replace_na(.w, min(.w, na.rm = TRUE)))
    }
    
    ## ---- parse endpoint marks ------------------------------------------------
    mark_of <- function(ch) dplyr::case_when(
        ch %in% c("<", ">") ~ "arrow",
        ch %in% c("o", "@") ~ "circle",
        TRUE                ~ "tail"
    )
    
    g <- g |> activate(edges) |> mutate(
        .m_from = mark_of(substr(type, 1, 1)),
        .m_to   = mark_of(substr(type, nchar(type), nchar(type))),
        .ends   = dplyr::case_when(
            .m_from == "arrow" & .m_to == "arrow" ~ "both",
            .m_to   == "arrow"                    ~ "last",
            .m_from == "arrow"                    ~ "first",
            TRUE                                  ~ "none"
        ),
        .style  = dplyr::case_when(
            .m_from == "circle" | .m_to == "circle" ~ "partial",
            .ends == "both"                         ~ "bidirected",
            .ends == "none"                         ~ "undirected",
            TRUE                                    ~ "directed"
        )
    )
    
    edge_df <- tibble::as_tibble(g, active = "edges")
    
    if (any(grepl(";", edge_df$type, fixed = TRUE)))
        warning("collapsed `type` values found - run gcm_split_types() first")
    
    ## ---- layout --------------------------------------------------------------
    lay <- if (all(c("x", "y") %in% names(node_df))) {
        ggraph::create_layout(g, layout = "manual", x = node_df$x, y = node_df$y)
    } else {
        ggraph::create_layout(g, layout = layout)
    }
    
    ## ---- circle marker positions, in data units ------------------------------
    seg <- edge_df |>
        mutate(x1 = lay$x[from], y1 = lay$y[from],
               x2 = lay$x[to],   y2 = lay$y[to],
               len = sqrt((x2 - x1)^2 + (y2 - y1)^2)) |>
        filter(len > 0)
    
    off <- if (is.null(mark_offset)) 0.15 * stats::median(seg$len) else mark_offset
    
    marks <- bind_rows(
        seg |> filter(.m_from == "circle") |>
            mutate(mx = x1 + off / len * (x2 - x1),
                   my = y1 + off / len * (y2 - y1)),
        seg |> filter(.m_to == "circle") |>
            mutate(mx = x2 - off / len * (x2 - x1),
                   my = y2 - off / len * (y2 - y1))
    ) |> select(mx, my, .w)
    
    ## ---- layers --------------------------------------------------------------
    cap <- ggraph::circle(node_size / 2 + 1.5, "mm")
    ar  <- function(e) grid::arrow(type = "closed",
                                   length = grid::unit(arrow_mm, "mm"), ends = e)
    
    link <- function(e) geom_edge_link(
        aes(filter = .ends == e & .style != "bidirected",
            edge_colour = .style, edge_linetype = .style, edge_alpha = .w),
        arrow = if (e == "none") NULL else ar(e),
        start_cap = cap, end_cap = cap, edge_width = 0.5
    )
    
    ## identity when flat, so a constant column isn't rescaled to a midpoint
    alpha_edge <- if (use_w) {
        scale_edge_alpha_continuous(range = alpha_range, name = edge_weight_var)
    } else {
        scale_edge_alpha_identity()
    }
    alpha_mark <- if (use_w) {
        scale_alpha_continuous(range = alpha_range, guide = "none")
    } else {
        scale_alpha_identity()
    }
    
    ggraph(lay) +
        link("none") + link("last") + link("first") +
        ## bidirected arced, so it separates from a co-existing directed edge
        geom_edge_arc(
            aes(filter = .style == "bidirected",
                edge_colour = .style, edge_linetype = .style, edge_alpha = .w),
            strength = arc, arrow = ar("both"),
            start_cap = cap, end_cap = cap, edge_width = 0.5
        ) +
        geom_point(data = marks, aes(x = mx, y = my, alpha = .w),
                   shape = 21, fill = "white", colour = "grey30",
                   size = mark_size, stroke = 0.6, inherit.aes = FALSE) +
        geom_node_point(size = node_size, colour = "grey85") +
        geom_node_text(aes(label = name), size = label_size) +
        alpha_edge + alpha_mark +
        scale_edge_colour_manual(values = c(directed = "grey10", bidirected = "grey30",
                                            undirected = "grey30", partial = "grey30"),
                                 guide = "none") +
        scale_edge_linetype_manual(values = c(directed = "solid", bidirected = "dashed",
                                              undirected = "dotted", partial = "dotted"),
                                   guide = "none") +
        labs(title = title) +
        theme_graph(base_family = "sans") +
        theme(plot.margin = margin(2, 2, 2, 2))
}


## --- checks -----------------------------------------------------------------
## g <- tbl_graph(nodes = tibble(name = c("CD.Str","Grade","S.Att","S.Kn")),
##                edges = tibble(from = c(1,1,1,3,4), to = c(2,3,4,2,2),
##                               type = c("o->","o-o","o-o","o->","o->"),
##                               n = c(8, 1, 4, 6, 2)))
## gcm_qplot(g)                            # all edges at flat_alpha
## gcm_qplot(g, edge_weight_var = "n")     # n=1 at 0.3, n=8 at 0.95
## gcm_qplot(g, edge_weight_var = "nope")  # warns, falls back to flat
## --- checks -----------------------------------------------------------------
## tbl_graph(nodes = tibble(name = c("CD.Str","Grade","S.Att","S.Kn")),
##           edges = tibble(from = c(1,1,1,3,4), to = c(2,3,4,2,2),
##                          type = c("o->","o-o","o-o","o->","o->"))) |>
##   gcm_qplot()
##   -> 5 dotted grey edges; open circles at CD.Str on all three of its edges,
##      at S.Att and S.Kn on their own ends, arrowheads into Grade

## --- checks -----------------------------------------------------------------
## g <- gcm_projection(dagitty("dag{L->X L->Y X->Y}"), latent_nodes = "L")
## gcm_plot(g)
##   both X->Y (straight) and X<->Y (dashed arc) must be visible and distinct
##
## gcm_latent_projection(dagitty("dag{L1->L2 L2->X L2->Y L1->Z}"),
##                       latent_nodes = c("L1","L2")) |> gcm_plot()
##   three dashed arcs, no solid lines

# ===================================== #
# Graph projections ---------------------
# ===================================== #



#' Latent projection of a DAG onto a set of nodes
#'
#' Verma-Pearl projection. Returns an ADMG over the retained nodes:
#'   x -> y   iff a directed path x ~> y has all intermediate nodes latent
#'   x <-> y  iff x and y are both reachable from some latent node via
#'            paths whose intermediate nodes are all latent
#'
#' Both edges may hold for the same pair; that is correct, not a duplicate.
#' Edge type is carried in the `type` column ("->" / "<->"), since igraph
#' stores every edge directionally.

gcm_projection <- function(
        g,
        # Either use latent (latent projection) or keep (projection onto)
        latent_nodes = NULL,
        keep_nodes   = NULL,
        noisy        = FALSE) {
    
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    ig <- igraph::as.igraph(g)
    
    all_nodes <- igraph::V(g)$name
    
    if (is.null(keep_nodes)) {
        observed <- setdiff(all_nodes, latent_nodes)
    } else {
        absent <- setdiff(keep_nodes, all_nodes)
        if (length(absent) && noisy)
            message("not in graph, ignored: ", paste(absent, collapse = ", "))
        observed <- intersect(all_nodes, keep_nodes)
    }
    latent <- setdiff(all_nodes, observed)
    
    if (length(latent) == 0) {
        if (noisy) message("no latent nodes; returning original graph")
        if (!("type" %in% edge_attr_names(g))) {
            g <- g %E>%
                mutate(type = "->")
        }
        return(g)
    }
    
    ## observed nodes reachable from `v` through all-latent intermediates.
    ## traversal stops on reaching an observed node, so intermediates stay latent.
    reach_obs <- function(v) {
        frontier <- igraph::neighbors(g, v, mode = "out")$name
        found    <- character(0)
        seen     <- character(0)
        while (length(frontier)) {
            found    <- union(found, intersect(frontier, observed))
            step     <- setdiff(intersect(frontier, latent), seen)
            seen     <- union(seen, step)
            frontier <- unique(unlist(lapply(
                step, function(u) igraph::neighbors(g, u, mode = "out")$name
            )))
        }
        setdiff(found, v)
    }
    
    dir_edges <- observed |>
        map(\(x) tibble(from = x, to = reach_obs(x), type = "->")) |>
        list_rbind()
    
    bi_edges <- latent |>
        map(\(l) {
            s <- sort(reach_obs(l))
            if (length(s) < 2) return(NULL)
            as_tibble(t(combn(s, 2)), .name_repair = ~ c("from", "to")) |>
                mutate(type = "<->")
        }) |>
        list_rbind()
    
    edges <- unique(rbind(dir_edges, bi_edges))
    nodes <- data.frame(name = observed, stringsAsFactors = FALSE)
    
    if (is.null(edges) || nrow(edges) == 0) {
        return(tidygraph::tbl_graph(nodes = nodes, directed = TRUE))
    }
    
    tidygraph::tbl_graph(
        nodes    = nodes,
        edges    = data.frame(from = match(edges$from, observed),
                              to   = match(edges$to,   observed),
                              type = edges$type,
                              stringsAsFactors = FALSE),
        directed = TRUE
    )
}


## --- checks -----------------------------------------------------------------
# gcm_projection(dagitty("dag{L1->L2 L2->X L2->Y L1->Z}"),
#                       latent_nodes = c("L1","L2"))
##   -> X<->Y, X<->Z, Y<->Z            (3 edges, none directed)
##
# gcm_projection(dagitty("dag{L->X L->Y X->Y}"), latent_nodes = "L")
#   -> X->Y and X<->Y                 (both, on the same pair)
##
## gcm_latent_projection(dagitty("dag{X->L L->Y}"), latent_nodes = "L")
##   -> X->Y                           (chain, no confounding)
##
## gcm_latent_projection(dagitty("dag{X->L Y->L L->Z}"), latent_nodes = "L")
##   -> X->Z, Y->Z                     (collider: X,Y NOT confounded)


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



gcm_to_cpdag <- function(g) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    m <- as_adjacency_matrix(g, sparse = FALSE)
    cpdag_m <- pcalg::dag2cpdag(m)
    # dag2cpdag drops dimnames, so restore from original
    dimnames(cpdag_m) <- dimnames(m)
    igraph::graph_from_adjacency_matrix(cpdag_m, mode = "directed") |>
        as_tbl_graph()
}




# =============================== #
# ===== Clustering DAGs ==========
# =============================== #

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
    
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    
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



# ========================================= #
# ======== Causal fingerprints ===========
# ========================================= #


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


gcm_ci_fingerprint <- function(g, vocab = NULL, max_cond = Inf) {
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
                    rows[[length(rows) + 1]] <- tibble(
                        x = p[1], y = p[2],
                        # It is critical that Z is sorted here!!!
                        Z = list(sort(Z)))
                }
            }
        }
    }

    if (length(rows) == 0) {
        tibble(x = character(), y = character(), Z = list())
    } else {
        return(bind_rows(rows))
    }
}

# ---- 1. Adjustment set for one pair --------------------------------

gcm_adj_set <- function(g, exposure, outcome,
                        adj_type = c("canonical", "all" ,"optimal" )) {
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

gcm_adj_fingerprint <- function(g, vocab = NULL,
                                adj_type = c("canonical", "all" ,"optimal" )) {
    # tibble: x, y, z (list of sets/NULLs), id, causal
    # `causal` = x is an ancestor of y in THIS graph
    adj_type <- match.arg(adj_type)
    if (!inherits(g, "dagitty")) g <- tidygraph_to_dagitty(g)
    
    if (is.null(vocab)) vocab <- names(g)
    stopifnot(all(vocab %in% names(g)))
    vocab <- sort(vocab)
    
    pairs <- expand_grid(x = vocab, y = vocab) |>
        filter(x != y)
    
    an <- map(set_names(names(g)), \(v) ancestors(g, v))  # precompute
    
    pairs |>
        mutate(
            Z      = map2(x, y, \(x, y) gcm_adj_set(g, x, y, adj_type)),
            id     = !map_lgl(Z, is.null),
            causal = map2_lgl(x, y, \(x, y) x %in% an[[y]])
        )
}

# Unconditional (marginal) independence, from a graph -------------
gcm_marginal_independence <- function(g, S = NULL) {
    if (inherits(g, "tbl_graph")) g <- tidygraph_to_dagitty(g)
    if (is.null(S)) {
        S <- gcm_nodelist(g, as_tibble = FALSE)
    }
    
    pairs <- expand_grid(x = S, y = S) |> 
        filter(x != y)
    
    pairs |>
        mutate(
            independent = map2_lgl(
                x, y,
                \(xx, yy) dagitty::dseparated(g, xx, yy, character(0))
            )
        )
}


# ========================================= #
# ======== Causal dissimilarities ===========
# ========================================= #
# These are for examining causal consistency between
# any two graphs. 

# L0 dissimilarities ----------
L0_dissimilarity_vocab <- function(gX, gY, return_list = FALSE, kappa = 1) {
    if (inherits(gX, "tbl_graph")) gX <- tidygraph_to_dagitty(gX)
    if (inherits(gY, "tbl_graph")) gY <- tidygraph_to_dagitty(gY)
    
    vocab_shared <- intersect(names(gX), names(gY))
    vocab_diff <- symdiff(names(gX), names(gY))
    vocab_union <- union(names(gX), names(gY))
    
    d <- length(vocab_diff) / (length(vocab_union) + kappa)
    
    if (return_list) {
        list(
            d = d,
            S = vocab_shared,
            full_vocab = vocab_union
        )
    } else {
        return(d)
    }
}


L1_dissimilarity_ci <- function(
        gX = NULL, gY = NULL,
        fpX = NULL, fpY = NULL,
        return_list = FALSE,
        max_cond = Inf,
        restriction_to_shared = TRUE,
        kappa = 1,
        include_separability = TRUE,
        return_marginal_independence = FALSE) {
    
    S <- intersect(gcm_nodelist(gX), gcm_nodelist(gY)) |> pull(name)
    
    if (is.null(fpX)) fpX <- gcm_ci_fingerprint(gX, vocab = S, max_cond = max_cond)
    if (is.null(fpY)) fpY <- gcm_ci_fingerprint(gY, vocab = S, max_cond = max_cond)
    
    if (include_separability) {
        fpX <- bind_rows(fpX, fpX |> distinct(x, y) |> mutate(Z = list("CI")))
        fpY <- bind_rows(fpY, fpY |> distinct(x, y) |> mutate(Z = list("CI")))
    }
    
    if (restriction_to_shared) {
        Splus <- union(S, "CI")
        fpX <- fpX |> filter(x %in% S, y %in% S, map_lgl(Z, \(z) all(z %in% Splus)))
        fpY <- fpY |> filter(x %in% S, y %in% S, map_lgl(Z, \(z) all(z %in% Splus)))
    }
    
    gamma_diff  <- dplyr::symdiff(fpX, fpY)
    gamma_union <- dplyr::union(fpX, fpY)
    d <- nrow(gamma_diff) / (nrow(gamma_union) + kappa)
    
    marg_X <- marg_Y <- NULL
    if (return_marginal_independence) {
        marg_X <- gcm_marginal_independence(fpX, S)
        marg_Y <- gcm_marginal_independence(fpY, S)
    }
    
    if (!return_list) {
        if (return_marginal_independence) {
            return(list(d = d, marginal_independence_X = marg_X, marginal_independence_Y = marg_Y))
        }
        return(d)
    } else {
        out <- list(d = d, gammaX = fpX, gammaY = fpY, S = S)
        if (return_marginal_independence) {
            out$marginal_independence_X <- marg_X
            out$marginal_independence_Y <- marg_Y
        }
        return(out)
    }
}

# L1 dissimilarities ----------------
L1_dissimilarity_ci <- function(
        gX = NULL, gY = NULL,
        fpX = NULL, fpY = NULL,
        return_list = FALSE,
        max_cond = Inf,
        restriction_to_shared = TRUE,
        kappa = 1,
        # include the binary as well - maybe change to indep
        include_separability = TRUE,
        include_marginal_dependence = TRUE) {

    S <- intersect(gcm_nodelist(gX),gcm_nodelist(gY)) |> pull(name)

    if (is.null(fpX)) {
        fpX <- gcm_ci_fingerprint(gX, vocab = S, max_cond = max_cond)
    }

    if (is.null(fpY)) {
        fpY <- gcm_ci_fingerprint(gY, vocab = S, max_cond = max_cond)
    }


    if (include_separability) {
        fpX <-
            bind_rows(
                fpX,
                fpX |>
                    distinct(x,y) |>
                    mutate(Z = list("CI"))
            )
        fpY <-
            bind_rows(
                fpY,
                fpY |>
                    distinct(x,y) |>
                    mutate(Z = list("CI"))
            )
    }
    
    if (include_marginal_dependence) {
        marg_X <- gcm_marginal_independence(gX, S) |>
            filter(!independent) |> 
            distinct(x,y) |> 
            mutate(Z = list("Marginally dependent"))
        marg_Y <- gcm_marginal_independence(gY, S) |> 
            filter(!independent) |> 
            distinct(x,y) |> 
            mutate(Z = list("Marginally dependent"))
        
        fpX <- bind_rows(fpX, marg_X)
        fpY <- bind_rows(fpY, marg_Y)
    }

    if (restriction_to_shared) {
        Splus <- union(S, "CI")
        fpX <- fpX |>
            filter(
                x %in% S,
                y %in% S,
                map_lgl(Z, \(z) all(z %in% Splus))
            )
        fpY <- fpY |>
            filter(
                x %in% S,
                y %in% S,
                map_lgl(Z, \(z) all(z %in% Splus))
            )
    }


    gamma_diff <- dplyr::symdiff(fpX, fpY)
    gamma_union <- dplyr::union(fpX, fpY)

    d <- nrow(gamma_diff) / (nrow(gamma_union) + kappa)

    if (!return_list) {
        return(d)
    } else {
        return(list(
            d = d, gammaX = fpX, gammaY = fpY, S = S
        ))
    }
}

L1_dissimilarity_ci_latent_projection <- 
    function(
        gX, gY, return_list = FALSE,
        max_cond = Inf, kappa = 1,
        include_separability = TRUE) {
        S <- intersect(gcm_nodelist(gX),gcm_nodelist(gY)) |> pull(name)
        
        gXproj <- gcm_latent_projection(gX, keep_nodes = S)
        gYproj <- gcm_latent_projection(gY, keep_nodes = S)
        
        fpX <- gcm_ci_fingerprint(gXproj, max_cond = max_cond)
        fpY <- gcm_ci_fingerprint(gYproj, max_cond = max_cond)
        
        
        if (include_separability) {
            fpX <- 
                bind_rows(
                    fpX,
                    fpX |> 
                        distinct(x,y) |> 
                        mutate(Z = list("CI"))
                )
            fpY <- 
                bind_rows(
                    fpY,
                    fpY |> 
                        distinct(x,y) |> 
                        mutate(Z = list("CI"))
                )
        }
        
        gamma_diff <- dplyr::symdiff(fpX, fpY)
        gamma_union <- dplyr::union(fpX, fpY)
        
        # +1 includes empty set ("no claim" is available to both)
        d <- nrow(gamma_diff) / (nrow(gamma_union) + kappa)
        
        if (!return_list) {
            return(d)
        } else {
            return(list(
                d = d, gammaX = fpX, gammaY = fpY, S = S
            ))
        }
    }

L2_dissimilarity_adj_set <- function(
        gX = NULL, gY = NULL,
        fpX = NULL, fpY = NULL,
        return_list = FALSE,
        kappa = 1,
        restriction_to_shared = TRUE,
        restriction_to_causal = FALSE,
        include_identifiability = TRUE,
        include_ancestral = TRUE,
        adj_type = "canonical") {
    
    S <- intersect(gcm_nodelist(gX),gcm_nodelist(gY)) |> pull(name)
    
    
    if (is.null(fpX)) {
        fpX <- gcm_adj_fingerprint(gX, adj_type = adj_type, vocab = S)
    }    
        
    if (is.null(fpY)) {
        fpY <- gcm_adj_fingerprint(gY, adj_type = adj_type, vocab = S)
    }
    
    if (include_identifiability) {
        fpX <- 
            bind_rows(
                fpX,
                fpX |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Identifiable"))
            )
        fpY <- 
            bind_rows(
                fpY,
                fpY |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Identifiable"))
            )
    }
    
    if (include_ancestral) {
        fpX <- 
            bind_rows(
                fpX,
                gcm_descendants_fingerprint(gX) |> 
                    filter(causal_path) |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Ancestral"))
            )
        fpY <- 
            bind_rows(
                fpY,
                gcm_descendants_fingerprint(gY) |> 
                    filter(causal_path) |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Ancestral"))
            )
    }
    
    # Might not be needed now vocab is being used
    if (restriction_to_shared) {
        Splus <- union(S, "Identifiable")
        fpX <- fpX |>
            filter(
                x %in% S,
                y %in% S,
                map_lgl(Z, \(z) all(z %in% Splus))
            )
        fpY <- fpY |>
            filter(
                x %in% S,
                y %in% S,
                map_lgl(Z, \(z) all(z %in% Splus))
            )
    }
    
        
    gamma_diff <- dplyr::symdiff(fpX, fpY)
    gamma_union <- dplyr::union(fpX, fpY)
    
    
    d <- nrow(gamma_diff) / (nrow(gamma_union) + kappa)
    
    if (!return_list) {
        return(d)
    } else {
        return(list(
            d = d, gammaX = fpX, gammaY = fpY
        ))
    }
}

L2_dissimilarity_adj_set_latent_projection <- function(
        gX, gY, 
        return_list = FALSE,
        include_identifiability = TRUE,
        kappa = 1,
        adj_type = "canonical") {
    S <- intersect(gcm_nodelist(gX),gcm_nodelist(gY)) |> pull(name)
    
    gXproj <- gcm_latent_projection(gX, keep_nodes = S)
    gYproj <- gcm_latent_projection(gY, keep_nodes = S)
    
    fpX <- gcm_adj_fingerprint(gXproj, adj_type = adj_type)
    fpY <- gcm_adj_fingerprint(gYproj, adj_type = adj_type)
    
    if (include_identifiability) {
        fpX <- 
            bind_rows(
                fpX,
                fpX |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Identifiable"))
            )
        fpY <- 
            bind_rows(
                fpY,
                fpY |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Identifiable"))
            )
    }
    
    gamma_diff <- dplyr::symdiff(fpX, fpY)  
    gamma_union <- dplyr::union(fpX, fpY)  
    
    
    d <- nrow(gamma_diff) / (nrow(gamma_union) + kappa) 
    
    if (!return_list) {
        return(d)
    } else {
        return(list(
            d = d, gammaX = fpX, gammaY = fpY
        ))
    }
}


# ======================================= #
# ======= Cluster degradation ===== 
# ======================================= #


coarsen_vocab <- function(d, .f, sep = ",") {
    # This words on a tibble with x,y,Z (char, char, list)   
    # which is the output format of L1 and L2 fingerprints
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
        ) |> filter(
            # Drop statements that become internal to a cluster
            x != y,
            # Drop statements whose Z set reaches into either endpoint
            !map2_lgl(x, Z, \(a, s) a %in% s),
            !map2_lgl(y, Z, \(b, s) b %in% s)
        ) |>
        distinct()
}


gcm_cluster_L1_degradation <- function(
        g, .f, fp = NULL,
        noisy = TRUE, max_cond = Inf,
        include_separability = TRUE,
        kappa = 1) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    
    CDAG  <- gcm_cluster_graph(g, .f)
    
    if (is.null(fp)) {
        fp_fine   <- gcm_ci_fingerprint(
            g, 
            max_cond = max_cond)      
    } else {
        fp_fine <- fp
    }
    
    fp_fine_coarse_vocab <- # lifted claims
        fp_fine |>
        coarsen_vocab(.f = .f)
    
    fp_coarse <- gcm_ci_fingerprint(
        CDAG, max_cond = max_cond) # coarse graph claims
    
    if (include_separability) {
        fp_fine_coarse_vocab <- 
            bind_rows(
                fp_fine_coarse_vocab,
                fp_fine_coarse_vocab |> 
                    distinct(x,y) |> 
                    mutate(Z = list("CI"))
            )
        fp_coarse <- 
            bind_rows(
                fp_coarse,
                fp_coarse |> 
                    distinct(x,y) |> 
                    mutate(Z = list("CI"))
            )
    }
    
    diff_claims <- dplyr::symdiff(fp_coarse, fp_fine_coarse_vocab)
    all_claims <- dplyr::union(fp_fine_coarse_vocab, fp_coarse)
    
    degradation <- nrow(diff_claims) / (nrow(all_claims) + kappa)
    
    if (noisy) message(str_c("L1 degradation of ", degradation, " computed."))
    
    list(
        degradation = degradation,
        symdiff = diff_claims,
        all_claims = all_claims,
        claims_in_coarse = fp_coarse,
        claims_in_fine = fp_fine_coarse_vocab,
        full_claims_in_fine = fp_fine
    )
}


gcm_cluster_L2_degradation <- function(
        g, .f, noisy = TRUE,
        fp = NULL,
        adj_type = "canonical",
        include_separability = TRUE,
        kappa = 1) {
    if (inherits(g, "dagitty")) g <- dagitty_to_tidygraph(g)
    
    CDAG  <- gcm_cluster_graph(g, .f)
    
    if (is.null(fp)) {
        fp   <- gcm_adj_fingerprint(g, adj_type = adj_type)
    } 
    fp_fine <- fp |> 
        filter(causal) |> select(-causal) 
    fp_coarse <- gcm_adj_fingerprint(CDAG, adj_type = adj_type) |>
        filter(causal) |> select(-causal) # quotient claims
    fp_fine_coarse_vocab <- fp_fine |> 
        coarsen_vocab(.f = .f) # lifted claims
    
    if (include_separability) {
        fp_fine_coarse_vocab <- 
            bind_rows(
                fp_fine_coarse_vocab,
                fp_fine_coarse_vocab |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Identifiable"))
            )
        fp_coarse <- 
            bind_rows(
                fp_coarse,
                fp_coarse |> 
                    distinct(x,y) |> 
                    mutate(Z = list("Identifiable"))
            )
    }
    
    diff_claims <- dplyr::symdiff(fp_fine_coarse_vocab, fp_coarse)
    all_claims <- dplyr::union(fp_fine_coarse_vocab, fp_coarse)
    
    degradation <- nrow(diff_claims) / (nrow(all_claims) + kappa)
    if (noisy) message(str_c("L2 degradation of ", degradation, " computed."))
    list(
        degradation = degradation,
        symdiff = diff_claims,
        all_claims = all_claims,
        claims_in_coarse = fp_coarse,
        claims_in_fine = fp_fine_coarse_vocab,
        full_claims_in_fine = fp_fine
    )
}







