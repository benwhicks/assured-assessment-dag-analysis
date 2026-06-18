# ===================================================================
# Validating and selecting graph coarsenings (abstraction maps)
# for comparing DAGs at L1/L2 across granularities.
#
# Companion to compare_dags_L1_L2.R -- source that first for
# tidy_to_dagitty(), project_to(), ci_compare(), l2_summary(), etc.
#
# Design:
#   Gate 1  check_admissible()        quotient must be acyclic (C-DAG)
#   Gate 2  cluster_ci_audit()        soundness is guaranteed by theory
#                                     (Anand et al. 2023); the audit
#                                     verifies it and quantifies LOSS
#   Costs   ancestral_inflation()     spurious coarse ancestral relations
#           cluster_adjustment_audit() lost adjustment-identifiability
#   Search  score_map_pair()          loss of a candidate (f1, f2)
#           search_coarsening()       constrained greedy merge search
#
# Key references:
#   Anand, Ribeiro, Tian & Bareinboim (2023, AAAI)  C-DAGs: soundness of
#       d-separation and identification at cluster level
#   Tikka, Helske & Karvanen (2023, JMLR)  transit clusters: which
#       clusterings preserve identifiability
#   Beckers & Halpern (2019, AAAI); Beckers, Eberhardt & Halpern (2019,
#       UAI)  constructive / approximate abstraction
#   Rischel & Weichwald (2021, UAI)  compositional abstraction error
#   Zennaro et al. (2023, CLeaR)  learning abstraction maps from data
# ===================================================================

library(tidygraph)
library(igraph)
library(dagitty)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

# source("compare_dags_L1_L2.R")

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

g_edgelist <- function(g) {
  nm <- g |> activate(nodes) |> pull(name)
  g |> 
      activate(edges) |> 
      as_tibble() |>
      mutate(
          from = nm[from], 
          to = nm[to], 
          .keep = "none")
}

#' Build a vectorized string->string map from a named lookup vector
#' (names = fine nodes, values = coarse names). Unmapped names pass
#' through unchanged -- convenient for partial / incremental maps.
map_from_lookup <- function(lookup) {
  force(lookup)
  function(x) ifelse(x %in% names(lookup), unname(lookup[x]), x)
}

#' Cluster membership induced by .f on a graph's nodes
memberships <- function(g, .f) {
  nm <- g |> activate(nodes) |> pull(name)
  split(nm, .f(nm))
}

# -------------------------------------------------------------------
# Gate 1: admissibility (acyclic quotient)
# -------------------------------------------------------------------

#' A coarsening is admissible as a C-DAG iff the quotient is acyclic.
#' Returns the specific obstructions: 2-cycles (the common case; these
#' are what would otherwise masquerade as bidirected edges) and any
#' larger strongly connected components.
check_admissible <- function(g, .f) {
  nm <- g |> activate(nodes) |> pull(name)
  el <- g_edgelist(g) |>
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

#' Quotient graph as a dagitty DAG (errors if inadmissible).
quotient_dagitty <- function(g, .f) {
  adm <- check_admissible(g, .f)
  if (!adm$admissible) {
    stop("Inadmissible coarsening. 2-cycles: ",
         paste(adm$two_cycles$a, "<->", adm$two_cycles$b, collapse = "; "),
         " | cyclic clusters: ",
         paste(adm$cyclic_clusters, collapse = ", "))
  }
  nm <- g |> activate(nodes) |> pull(name)
  cl <- unique(.f(nm))
  el <- g_edgelist(g) |>
    mutate(across(everything(), .f)) |>
    filter(from != to) |> distinct(from, to)
  dagitty(paste0(
    "dag {\n",
    paste(c(cl, paste(el$from, "->", el$to)), collapse = "\n"),
    "\n}"))
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

cluster_ci_audit <- function(g, .f, max_cond = Inf) {
  gd  <- tidy_to_dagitty(g)
  mem <- memberships(g, .f)
  qd  <- quotient_dagitty(g, .f) # not sure if this is right
  qd <- tidy_graph_to_dag(g_cluster_graph(g, .f))
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

ci_loss <- function(audit) {
  stopifnot(!any(audit$status == "VIOLATION"))   # unit test on theory
  mean(audit$status == "lost")
}

# -------------------------------------------------------------------
# Cost: ancestral inflation
# -------------------------------------------------------------------
# Fine ancestry always survives coarsening; the quotient can only ADD
# ancestral relations (via concatenated between-cluster paths that have
# no fine counterpart). One-directional by construction.

ancestral_inflation <- function(g, .f) {
  gd  <- tidy_to_dagitty(g)
  mem <- memberships(g, .f)
  # qd  <- quotient_dagitty(g, .f) # not sure about this one
  qd <- tidy_graph_to_dag(g_cluster_graph(g, .f))
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

cluster_adjustment_audit <- function(g, .f, type = "minimal") {
  gd  <- tidy_to_dagitty(g)
  mem <- memberships(g, .f)
  # qd  <- quotient_dagitty(g, .f)
  qd <- tidy_graph_to_dag(g_cluster_graph(g, .f))
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

# -------------------------------------------------------------------
# Two-graph comparison under candidate maps
# -------------------------------------------------------------------
# Pipeline per candidate (f1, f2):
#   coarsen each graph -> quotient DAGs over (partly) shared vocabulary
#   latent-project each quotient onto the shared cluster names
#   L1 + L2 comparison (functions from compare_dags_L1_L2.R)
#   penalize map complexity (number of merges)
#
# NB: coarsen-then-project and project-then-coarsen do not commute in
# general. Coarsening first keeps the abstraction map defined on actual
# model variables, which is the defensible order for your argument.

n_merges <- function(g, .f) {
  nm <- g |> activate(nodes) |> pull(name)
  length(nm) - length(unique(.f(nm)))
}

score_map_pair <- function(g1, g2, f1 = identity, f2 = identity,
                           lambda = 0.02, max_cond = Inf) {
  a1 <- check_admissible(g1, f1)
  a2 <- check_admissible(g2, f2)
  if (!a1$admissible || !a2$admissible) {
    return(tibble(admissible = FALSE, n_shared = NA, ci_jaccard = NA,
                  recipe_fail_12 = NA, recipe_fail_21 = NA,
                  merges = NA, loss = Inf))
  }
  q1 <- quotient_dagitty(g1, f1)
  q2 <- quotient_dagitty(g2, f2)
  shared <- intersect(names(q1), names(q2))
  if (length(shared) < 2) {
    return(tibble(admissible = TRUE, n_shared = length(shared),
                  ci_jaccard = NA, recipe_fail_12 = NA,
                  recipe_fail_21 = NA, merges = NA, loss = Inf))
  }
  m1 <- project_to(q1, shared)
  m2 <- project_to(q2, shared)

  ci  <- ci_compare(m1, m2, max_cond = min(max_cond, length(shared) - 2))
  s12 <- l2_summary(m1, m2)
  s21 <- l2_summary(m2, m1)
  k   <- n_merges(g1, f1) + n_merges(g2, f2)

  tibble(
    admissible     = TRUE,
    n_shared       = length(shared),
    ci_jaccard     = ci$jaccard,
    recipe_fail_12 = s12$prop_recipe_fails,
    recipe_fail_21 = s21$prop_recipe_fails,
    merges         = k,
    loss = (1 - ci$jaccard) +
           mean(c(s12$prop_recipe_fails, s21$prop_recipe_fails),
                na.rm = TRUE) +
           lambda * k
  )
}

# -------------------------------------------------------------------
# Constrained greedy search over coarsenings of g1
# -------------------------------------------------------------------
# candidates: a list of character vectors -- each vector is a set of
# fine nodes that are ALLOWED to share a cluster (semantic constraint,
# e.g. items belonging to the same survey construct). The search only
# proposes merges within a candidate group, into a cluster named after
# the group. This is what keeps the optimization from laundering
# genuine disagreement into "granularity difference".
#
# Returns the loss trajectory so you can report the decomposition:
#   loss(identity)           = naive disagreement
#   loss(identity) - min loss = granularity-attributable
#   min loss (net of penalty) = genuine residual disagreement

search_coarsening <- function(g1, g2, candidates,
                              f2 = identity, lambda = 0.02,
                              max_cond = Inf, quietly = FALSE) {
  stopifnot(is.list(candidates))
  if (is.null(names(candidates)))
    names(candidates) <- map_chr(candidates, ~ paste(.x, collapse = "_"))

  lookup <- character(0)                 # current (partial) map
  best <- score_map_pair(g1, g2, map_from_lookup(lookup), f2,
                         lambda = lambda, max_cond = max_cond)
  trajectory <- list(best |> mutate(step = 0, move = "identity"))

  repeat {
    # propose: add one not-yet-mapped node from some group into that
    # group's cluster (clusters named after the group)
    proposals <- imap(candidates, function(members, cluster) {
      pending <- setdiff(members, names(lookup))
      map(pending, function(node) {
        lk <- c(lookup, set_names(cluster, node))
        # ensure already-merged partners exist: first merge in a group
        # is a rename, which is fine (renaming = trivial merge)
        lk
      })
    }) |> flatten()

    if (length(proposals) == 0) break

    scores <- map_dfr(proposals, function(lk)
      score_map_pair(g1, g2, map_from_lookup(lk), f2,
                     lambda = lambda, max_cond = max_cond))
    i <- which.min(scores$loss)

    if (scores$loss[i] >= best$loss) break   # no improvement: stop

    lookup <- proposals[[i]]
    best   <- scores[i, ]
    moved  <- setdiff(names(lookup), names(trajectory)) # last added node
    trajectory[[length(trajectory) + 1]] <-
      best |> mutate(step = length(trajectory),
                     move = paste(tail(names(lookup), 1), "->",
                                  tail(unname(lookup), 1)))
    if (!quietly) message("step ", length(trajectory) - 1, ": ",
                          tail(names(lookup), 1), " -> ",
                          tail(unname(lookup), 1),
                          "  loss = ", round(best$loss, 4))
  }

  list(map        = map_from_lookup(lookup),
       lookup     = lookup,
       best       = best,
       trajectory = bind_rows(trajectory))
}

# ===================================================================
# Worked example
# ===================================================================

run_coarsening_example <- function() {
  # Fine graph: survey items i1..i3 measuring construct "Capability",
  # plus PD (professional development), Practice, Outcome
  g1 <- tbl_graph(
    nodes = tibble(name = c("PD", "i1", "i2", "i3",
                            "Practice", "Outcome")),
    edges = tribble(
      ~from, ~to,
      "PD", "i1", "PD", "i2", "PD", "i3",
      "i1", "Practice", "i2", "Practice", "i3", "Practice",
      "Practice", "Outcome"
    ) |> mutate(across(everything(),
            ~ match(.x, c("PD", "i1", "i2", "i3",
                          "Practice", "Outcome")))))

  # Coarse graph: someone else's model at construct level,
  # with an extra confounder "School"
  g2 <- tbl_graph(
    nodes = tibble(name = c("PD", "Capability", "Practice",
                            "Outcome", "School")),
    edges = tribble(
      ~from, ~to,
      "PD", "Capability", "Capability", "Practice",
      "Practice", "Outcome",
      "School", "PD", "School", "Outcome"
    ) |> mutate(across(everything(),
            ~ match(.x, c("PD", "Capability", "Practice",
                          "Outcome", "School")))))

  # Semantic constraint: only the items may merge, into "Capability"
  candidates <- list(Capability = c("i1", "i2", "i3"))

  # Audit the intended map on g1 first
  f_cap <- map_from_lookup(set_names(rep("Capability", 3),
                                     c("i1", "i2", "i3")))
  print(check_admissible(g1, f_cap))
  audit <- cluster_ci_audit(g1, f_cap)
  print(count(audit, status))
  print(cluster_adjustment_audit(g1, f_cap) |> filter(lost_id))

  # Search: how much disagreement is granularity vs genuine?
  res <- search_coarsening(g1, g2, candidates)
  print(res$trajectory)
  message("naive loss (identity): ",
          round(res$trajectory$loss[1], 4),
          " | best loss: ", round(res$best$loss, 4),
          " | granularity-attributable: ",
          round(res$trajectory$loss[1] - res$best$loss, 4))
  invisible(res)
}

# run_coarsening_example()
