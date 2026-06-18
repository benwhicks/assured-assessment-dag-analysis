# ===================================================================
# Comparing causal DAGs with partially overlapping node sets
# at L1 (observational/Markov) and L2 (interventional) levels
# of Pearl's causal hierarchy.
#
# Disclaimer: Used Claude Fable to generate initial code and ref, 
# before checking and adjusting. 
#
# Strategy: latent-project each DAG onto the shared variables
# (yielding MAGs), then compare the projections. Induced subgraphs
# do NOT represent the marginal models and must not be used.
#
# Key references:
#   Verma & Pearl (1990)            latent projection
#   Richardson & Spirtes (2002)     ancestral graph Markov models
#   Zhang (2008, JMLR)              MAGs/PAGs, FCI completeness
#   Ali, Richardson & Spirtes (2009) Markov equivalence of MAGs
#   Peters & Buhlmann (2015)        Structural Intervention Distance
#   Perkovic et al. (2018, JMLR)    complete generalized adjustment criterion
#   van der Zander et al. (2019)    adjustment-set algorithms (dagitty)
#
# Requires: dagitty (>= 0.3-1 for toMAG), tidygraph, igraph,
#           dplyr, purrr, tidyr; optionally pcalg, RBGL, graph, SID
# ===================================================================

library(dagitty)
library(tidygraph)
library(igraph)
library(dplyr)
library(purrr)
library(tidyr)

# -------------------------------------------------------------------
# 0. Interop: tidygraph <-> dagitty
# -------------------------------------------------------------------

#' tbl_graph (directed, with a `name` node attribute) -> dagitty DAG
tidy_to_dagitty <- function(tg) {
  stopifnot(inherits(tg, "tbl_graph"))
  nm <- tg |> activate(nodes) |> pull(name)
  el <- tg |> activate(edges) |> as_tibble()
  lines <- c(nm, if (nrow(el) > 0) paste(nm[el$from], "->", nm[el$to]))
  dagitty(paste0("dag {\n", paste(lines, collapse = "\n"), "\n}"))
}

#' dagitty graph (DAG or MAG) -> tbl_graph, edge `type` preserved
dagitty_to_tidy <- function(g) {
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

# -------------------------------------------------------------------
# 1. Latent projection onto the shared margin
# -------------------------------------------------------------------

#' Project a dagitty DAG onto `keep`, marginalizing everything else.
#' Returns a MAG over `keep` (bidirected edges = confounding through
#' the dropped variables).
project_to <- function(g, keep) {
  stopifnot(all(keep %in% names(g)))
  latents(g) <- setdiff(names(g), keep)
  toMAG(g)
}

# -------------------------------------------------------------------
# 2. L1: comparing the marginal independence models
# -------------------------------------------------------------------

# ---- 2a. Edge-mark Hamming distance between two MAGs --------------
# Syntactic distance in the spirit of SHD (Tsamardinos et al. 2006),
# extended to edge marks. NB: can be nonzero for Markov-equivalent
# MAGs; for an equivalence-class distance use the PAG route (2c) or
# the CI fingerprint (2b).

mark_matrix <- function(m) {
  v <- names(m); n <- length(v)
  M <- matrix("", n, n, dimnames = list(v, v))
  e <- dagitty::edges(m)
  for (k in seq_len(nrow(e))) {
    a <- e$v[k]; b <- e$w[k]
    marks <- switch(e$e[k],
      "->"  = c("-", ">"),   # tail at a, arrow at b
      "<->" = c(">", ">"),
      "--"  = c("-", "-"),
      stop("unhandled edge type: ", e$e[k]))
    M[a, b] <- marks[2]      # mark at b, as seen from a
    M[b, a] <- marks[1]      # mark at a, as seen from b
  }
  M
}

#' Hamming distance over unordered node pairs.
#' count_marks = TRUE: same skeleton but different orientation counts 1.
shd_mag <- function(m1, m2, 
                    count_marks = TRUE,
                    project_to_shared = TRUE) {
  v1 <- names(m1)
  v2 <- names(m2)
  v.shared <- intersect(v1, v2)
  if (!project_to_shared) stopifnot(setequal(v1, v2))
  if (project_to_shared) {
      m1 <- project_to(m1, v.shared)
      m2 <- project_to(m2, v.shared)
  }
  v <- names(m1)
  A <- mark_matrix(m1)[v, v]
  B <- mark_matrix(m2)[v, v]
  d <- 0L
  for (i in seq_along(v)) for (j in seq_along(v)) if (i < j) {
    adjA <- A[i, j] != ""; adjB <- B[i, j] != ""
    if (adjA != adjB) {
      d <- d + 1L
    } else if (adjA && count_marks &&
               (A[i, j] != B[i, j] || A[j, i] != B[j, i])) {
      d <- d + 1L
    }
  }
  d
}

# ---- 2b. Conditional-independence fingerprints ---------------------
# Enumerate ALL m-separation statements x _||_ y | Z over the shared
# variables. This characterizes the marginal Markov equivalence class:
# identical fingerprints <=> Markov-equivalent projections.
# Exponential in |shared| -- fine up to ~12-14 shared nodes.

# Note that it is important to use the dseparated version, the 
# adjustment sets dagitty function abbreviates the node names, 
# which is quite annoying for this kind of code!
ci_fingerprint <- function(m, vars = names(m),
                           max_cond = length(vars) - 2) {
  out <- tibble(
      X = "X", Y = "Y", Z = list("Z")
  ) |> 
      slice(0)
  
  for (p in combn(vars, 2, simplify = FALSE)) {
    rest <- setdiff(vars, p)
    for (k in 0:min(max_cond, length(rest))) {
      Zs <- if (k == 0) list(character(0))
            else combn(rest, k, simplify = FALSE)
      for (Z in Zs) {
        if (dseparated(m, p[1], p[2], Z)) {   # m-separation for MAGs
          out <- 
              bind_rows(out,
                        tibble(X = p[1], Y = p[2], Z = list(Z)))
              # c(out, sprintf("%s _||_ %s | {%s}",
              #                   p[1], p[2],
              #                   paste(sort(Z), collapse = ",")))
        }
      }
    }
  }
  distinct(out)
}

# TODO: Use CI fingerprint to do a similar version of fuzzy L2 but for L1
#.      Worth noting that this will be symmetric, whereas L2 is not. 
#.      Desiridata - needs to work with shared nodes in a similar way
#       to fuzzy_L2_consistency

fuzzy_L1_consistency <- function(gX, gY) {
    # Compares two graphs for L1 causal consistency
    # Graphs may have only some shared nodes.
    # worth returning the number of shared nodes!
}


ci_compare <- function(m1, m2, ...) {
    # should this be on the projection instead?
    # that would be more akin to the shd_mag
    # or save that for the sih_mag?
  f1 <- ci_fingerprint(m1, ...)
  f2 <- ci_fingerprint(m2, ...)
  list(
    n_g1        = length(f1),
    n_g2        = length(f2),
    n_common    = length(intersect(f1, f2)),
    jaccard     = length(intersect(f1, f2)) /
                  max(1, length(union(f1, f2))),
    markov_equiv = all_equal(f1, f2), 
    only_in_g1  = anti_join(f1, f2,
                            by = c("X","Y","Z")),
    only_in_g2  = anti_join(f2, f1,
                            by = c("X","Y","Z"))
  )
}

# ---- 2c. (Optional) PAG-level comparison via oracle FCI -------------
# Run FCI with a d-separation oracle on each FULL DAG, restricted to
# the shared variables. The output PAG is the marginal Markov
# equivalence class itself; mark-Hamming on PAGs is then an
# equivalence-class-level distance.
# Requires pcalg, graph, RBGL (Bioconductor).

dag_amat <- function(g) {
  v <- names(g)
  A <- matrix(0, length(v), length(v), dimnames = list(v, v))
  e <- dagitty::edges(g)
  stopifnot(all(e$e == "->"))        # full DAG only
  for (k in seq_len(nrow(e))) A[e$v[k], e$w[k]] <- 1
  A
}

oracle_pag <- function(g_full, observed) {
  if (!requireNamespace("pcalg", quietly = TRUE) ||
      !requireNamespace("RBGL",  quietly = TRUE)) {
    stop("oracle_pag() needs pcalg and RBGL")
  }
  v_all <- names(g_full)
  gnel  <- as(dag_amat(g_full), "graphNEL")
  suff  <- list(g = gnel, jp = RBGL::johnson.all.pairs.sp(gnel))
  idx   <- match(observed, v_all)
  itest <- function(x, y, S, suffStat)
    pcalg::dsepTest(idx[x], idx[y], idx[S], suffStat)
  pcalg::fci(suffStat = suff, indepTest = itest,
             alpha = 0.5, labels = observed)   # alpha irrelevant: oracle
}

pag_shd <- function(p1, p2) {
  A <- p1@amat; v <- colnames(A)
  B <- p2@amat[v, v]
  d <- 0L
  for (i in seq_along(v)) for (j in seq_along(v)) if (i < j) {
    if (A[i, j] != B[i, j] || A[j, i] != B[j, i]) d <- d + 1L
  }
  d
}

# -------------------------------------------------------------------
# 3. L2: comparing interventional implications
# -------------------------------------------------------------------

# ---- 3a. SID, when applicable ---------------------------------------
# SID (Peters & Buhlmann 2015) is defined for DAGs (and DAG vs CPDAG,
# with bounds) over a COMMON node set. After projection it applies
# only if both MAGs are bidirected-edge-free. Asymmetric by design.

is_pure_dag <- function(m) {
  e <- dagitty::edges(m)
  nrow(e) == 0 || all(e$e == "->")
}

sid_if_possible <- function(m1, m2) {
  if (!requireNamespace("SID", quietly = TRUE)) return(NA_real_)
  if (!is_pure_dag(m1) || !is_pure_dag(m2))      return(NA_real_)
  v  <- names(m1)
  A1 <- dag_amat(m1)[v, v]
  A2 <- dag_amat(m2)[v, v]
  c(sid_1_to_2 = SID::structIntervDist(A1, A2)$sid,
    sid_2_to_1 = SID::structIntervDist(A2, A1)$sid)
}

# ---- 3b. Adjustment-set cross-validity (SID generalized to MAGs) ----
# For each ordered pair (x, y) of shared nodes:
#   * does g_from identify P(y | do(x)) via covariate adjustment?
#   * if so, do its (minimal) adjustment sets remain valid in g_to
#     under the generalized adjustment criterion (Perkovic et al. 2018)?
# dagitty's adjustmentSets()/isAdjustmentSet() handle MAGs natively
# (van der Zander, Liskiewicz & Textor 2019).
#
# NB: adjustment-identifiability is not all of L2 identifiability
# (front-door / ID-algorithm constructions exist without adjustment
# sets). Frame this metric as agreement of adjustment recipes.

adjustment_cross_validity <- function(
        g_from, 
        g_to,
        vars = intersect(names(g_from),
                         names(g_to)),
        type = "minimal") {
    pairs <- expand_grid(x = vars, y = vars) |> 
        filter(x != y)
    pmap_dfr(pairs, function(x, y) {
    Z_from <- tryCatch(adjustmentSets(g_from, exposure = x, outcome = y,
                                      type = type),
                       error = function(e) list())
    Z_to   <- tryCatch(adjustmentSets(g_to, exposure = x, outcome = y,
                                      type = type),
                       error = function(e) list())
    id_from <- length(Z_from) > 0
    id_to   <- length(Z_to)   > 0
    if (!id_from) {
      return(tibble(x, y, identified_from = FALSE, identified_to = id_to,
                    n_sets_from = 0L, n_valid_in_to = NA_integer_,
                    any_valid = NA, set_jaccard = NA_real_))
    }
    valid <- map_lgl(Z_from, function(Z) {
        Z <- as.character(unlist(Z))
        if (length(Z) == 0) {
            return(dseparated(g_to, x, y, c()))
        }
        isAdjustmentSet(g_to, Z = Z, exposure = x, outcome = y)
    })
    # valid <- map_lgl(Z_from, function(Z)
    #   isAdjustmentSet(g_to, Z, exposure = x, outcome = y))
    canon <- function(L) map_chr(L, ~ paste(sort(unlist(.x)),
                                            collapse = ","))
    sj <- if (id_to) {
      a <- canon(Z_from); b <- canon(Z_to)
      length(intersect(a, b)) / length(union(a, b))
    } else NA_real_
    tibble(x, y,
           identified_from = TRUE, identified_to = id_to,
           n_sets_from   = length(Z_from),
           n_valid_in_to = sum(valid),
           any_valid     = any(valid),
           set_jaccard   = sj)
  })
}

#' Headline asymmetric scores, SID-flavoured:
#' proportion of ordered pairs where g_from's recipe FAILS in g_to.
l2_summary <- function(g_from, g_to, ...) {
  res <- adjustment_cross_validity(g_from, g_to, ...)
  res |>
    summarise(
      n_pairs            = n(),
      prop_id_agree      = mean(identified_from == identified_to),
      prop_recipe_fails  = mean(identified_from & !coalesce(any_valid, TRUE)),
      mean_set_jaccard   = mean(set_jaccard, na.rm = TRUE)
    )
}

# manual code

adj_fingerprint 
# getting the ZX, idX etc from function below

fuzzy_L2_consistency <- function(gX, gY) {
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
    if (!class(gX) == "dagitty") gX <- tidy_graph_to_dag(gX)
    if (!class(gY) == "dagitty") gY <- tidy_graph_to_dag(gY)
    
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
                canAdjSet(g1, x, y)
            ),
            canonical_ZY = pmap(
                list(x,y),
                \(x,y)
                canAdjSet(g2, x, y)
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
        select(x, y, idXY, idYX, ErrXY, everything())
    
    idXY <-  mean(df.out$idXY)
    idYX <-  mean(df.out$idYX)
    ErrXY <- mean(df.out$ErrXY)
    ErrYX <- mean(df.out$ErrYX)
    
    return(list(
        idXY = idXY, idYX = idYX,
        ErrXY = ErrXY, ErrYX = ErrYX,
        df = df.out, shared_nodes = shared
    ))
    
}

# ===================================================================
# 4. Worked example: partially overlapping node sets
# ===================================================================

run_example <- function() {
  # Graph 1: nodes X, M, Y, C, U   (U, M not in graph 2)
  g1 <- dagitty("dag {
    U -> X ; U -> M
    X -> M ; M -> Y
    X -> C ; C -> Y
  }")

  # Graph 2: nodes X, Y, C, W      (W not in graph 1)
  g2 <- dagitty("dag {
    X -> Y
    X -> C ; C -> Y
    W -> C ; W -> Y
  }")

  shared <- intersect(names(g1), names(g2))
  message("Shared nodes: ", paste(shared, collapse = ", "))

  m1 <- project_to(g1, shared)   # X -> Y (via M), X -> C, C -> Y; U adds X<->? (depends)
  m2 <- project_to(g2, shared)   # W induces C <-> Y
  print(m1); print(m2)

  # ---- L1 ----
  message("\nMAG edge-mark Hamming: ", shd_mag(m1, m2))
  ci <- ci_compare(m1, m2)
  message("CI Jaccard: ", round(ci$jaccard, 3),
          " | Markov equivalent: ", ci$markov_equiv)
  if (length(ci$only_in_g1))
    message("Claims only in g1:\n  ", paste(ci$only_in_g1, collapse = "\n  "))
  if (length(ci$only_in_g2))
    message("Claims only in g2:\n  ", paste(ci$only_in_g2, collapse = "\n  "))

  # ---- L2 ----
  s <- sid_if_possible(m1, m2)
  message("\nSID (if both projections are DAGs): ",
          if (all(is.na(s))) "not applicable (bidirected edges present)"
          else paste(names(s), s, collapse = "; "))

  cv12 <- adjustment_cross_validity(m1, m2)
  cv21 <- adjustment_cross_validity(m2, m1)
  print(cv12)
  message("\nL2 summary, g1 -> g2:"); print(l2_summary(m1, m2))
  message("L2 summary, g2 -> g1:");   print(l2_summary(m2, m1))

  invisible(list(m1 = m1, m2 = m2, ci = ci, cv12 = cv12, cv21 = cv21))
}

# run_example()
