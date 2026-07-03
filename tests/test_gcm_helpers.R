# ===================================================================
# test-gcm-consistency.R
#
# Regression suite for the L1/L2 consistency framework.
# Run with: testthat::test_file("test-gcm-consistency.R")
#
# Each expectation is hand-derived; derivations are in the comments so
# the tests double as worked examples for the paper.
#
# Assumed API (adjust names if yours differ):
#   gcm_L1_consistency(gX, gY, max_cond)      -> $dXY $dYX $jaccard
#                                                $markov_equiv_on_margin
#   gcm_L2_consistency(gX, gY, ...)           -> $dXY $dYX $rXY $rYX
#   gcm_cluster_L1_degradation(g, .f)         -> $d_degradation $n_lost
#                                                $lost_claims $audit
#   gcm_cluster_L2_degradation(g, .f)         -> $d_degradation $ErrLostId
#                                                $d_coarse_to_fine
#   gcm_cluster_check_admissible(g, .f)       -> $admissible $two_cycles
# ===================================================================

library(testthat)
library(dagitty)
# source("gcm_functions.R")   # <- your consolidated function file(s)

tol <- 1e-10

# -------------------------------------------------------------------
# 1. Identity: identical graphs are distance zero at both rungs
# -------------------------------------------------------------------

test_that("identical graphs: L1 and L2 distances are zero", {
    g <- dagitty("dag{ Z -> X ; Z -> Y ; X -> M ; M -> Y }")
    
    l1 <- gcm_L1_consistency(g, g)
    expect_equal(l1$dXY, 0, tolerance = tol)
    expect_equal(l1$dYX, 0, tolerance = tol)
    expect_equal(l1$jaccard, 1, tolerance = tol)
    expect_true(l1$markov_equiv_on_margin)
    
    l2 <- gcm_L2_consistency(g, g)
    expect_equal(l2$dXY, 0, tolerance = tol)
    expect_equal(l2$dYX, 0, tolerance = tol)
    # negative claims exist (e.g. reversed pairs) but none are refuted:
    expect_true(is.na(l2$rXY) || l2$rXY == 0)
    expect_true(is.na(l2$rYX) || l2$rYX == 0)
})

# -------------------------------------------------------------------
# 2. Minimal Markov-equivalence separator: A -> B  vs  B -> A
# -------------------------------------------------------------------
# L1: the only claim over {A, B} is A _||_ B | {}, denied by both
#     (adjacent nodes). Empty fingerprints = agreed saturation. d = 0,
#     Markov equivalent (same skeleton, no colliders).
# L2: gX identifies (A, B) with the empty recipe; (B, A) is NOT
#     adjustment-identifiable in gX (back-door B <- A blockable only
#     by conditioning on the outcome). Mirror image in gY. Each
#     direction: one contradiction out of one scoreable pair -> d = 1;
#     one negative claim, refuted -> r = 1.
# This is the minimal L1 = 0 / L2 = 1 witness.

test_that("edge reversal: L1 identical, L2 maximal", {
    gX <- dagitty("dag{ A -> B }")
    gY <- dagitty("dag{ B -> A }")
    
    l1 <- gcm_L1_consistency(gX, gY)
    expect_equal(l1$dXY, 0, tolerance = tol)
    expect_equal(l1$dYX, 0, tolerance = tol)
    expect_true(l1$markov_equiv_on_margin)
    
    l2 <- gcm_L2_consistency(gX, gY)
    expect_equal(l2$dXY, 1, tolerance = tol)
    expect_equal(l2$dYX, 1, tolerance = tol)
    expect_equal(l2$rXY, 1, tolerance = tol)
    expect_equal(l2$rYX, 1, tolerance = tol)
})

# -------------------------------------------------------------------
# 3. The triangle pair: Markov-equivalent, graded L2 disagreement
# -------------------------------------------------------------------
# gX: X -> Y with confounder Z; gY: Y -> X with confounder Z.
# Complete graphs: no CIs at all -> L1 = 0 via agreed saturation.
# L2 hand-derivation (canonical sets, pooled scorer, X -> Y direction):
#   (X,Y): gX id with {Z}, gY not id      -> contradiction, delta=1, w=1
#   (Y,X): gX not id, gY id               -> excluded (feeds r)
#   (Z,X), (Z,Y): both id, empty sets     -> 0/1 each
#   (X,Z), (Y,Z): both not id             -> 0/1 each
# unrestricted: d = 1/5 = 0.2 ; causal filter drops (X,Z),(Y,Z): d = 1/3
# r: three negative claims in gX ((Y,X),(X,Z),(Y,Z)), one refuted -> 1/3

test_that("triangle pair: L1 zero, L2 graded, filter changes denominator", {
    gX <- dagitty("dag{ X -> Y ; Z -> X ; Z -> Y }")
    gY <- dagitty("dag{ Y -> X ; Z -> X ; Z -> Y }")
    
    l1 <- gcm_L1_consistency(gX, gY)
    expect_equal(l1$dXY, 0, tolerance = tol)
    expect_true(l1$markov_equiv_on_margin)
    
    l2_all <- gcm_L2_consistency(gX, gY, causal_pairs_only = FALSE)
    expect_equal(l2_all$dXY, 0.2, tolerance = tol)
    expect_equal(l2_all$dYX, 0.2, tolerance = tol)
    expect_equal(l2_all$rXY, 1/3, tolerance = tol)
    
    l2_c <- gcm_L2_consistency(gX, gY, causal_pairs_only = TRUE)
    expect_equal(l2_c$dXY, 1/3, tolerance = tol)
    expect_equal(l2_c$dYX, 1/3, tolerance = tol)
})

# -------------------------------------------------------------------
# 4. Degradation under the identity coarsening is zero
# -------------------------------------------------------------------
# The oracle that would have caught the tibble-masking bug: with
# .f = identity the quotient IS the graph, so the audit's two answer
# columns must be identical row-by-row, and every metric must be 0.

test_that("identity coarsening: zero degradation at L1 and L2", {
    g <- tidygraph::as_tbl_graph(
        igraph::graph_from_literal(Z -+ X, Z -+ Y, X -+ M, M -+ Y))
    
    d1 <- gcm_cluster_L1_degradation(g, identity)
    expect_equal(d1$d_degradation, 0, tolerance = tol)
    expect_equal(d1$n_lost, 0L)
    expect_identical(d1$audit$fine_indep, d1$audit$coarse_indep)
    expect_false(any(d1$audit$status == "VIOLATION"))
    
    d2 <- gcm_cluster_L2_degradation(g, identity)
    expect_equal(d2$d_degradation, 0, tolerance = tol)
    expect_equal(d2$ErrLostId, 0, tolerance = tol)
    expect_true(is.na(d2$d_coarse_to_fine) || d2$d_coarse_to_fine == 0)
})

# -------------------------------------------------------------------
# 5. Spurious-path merge: hand-computed lossy coarsening
# -------------------------------------------------------------------
# Fine graph: A -> B1 ; B2 -> C   (two disconnected chains)
# Coarsening: {B1, B2} -> B       Quotient: A -> B -> C
#
# Cluster-level L1 claims (clusters A, B, C):
#   fine (set-level d-separation):  A _||_ C | {}   and   A _||_ C | B
#   coarse (quotient):              A _||_ C | B    only
#     (A -> B -> C is open marginally: the spurious path)
# => 1 of 2 fine claims lost: d_degradation = 1/2, lost claim is the
#    marginal one. Also the canonical spurious-ancestry example:
#    A is a coarse ancestor of C with no fine counterpart.

test_that("spurious-path merge: d1_deg = 1/2, marginal claim lost", {
    g  <- tidygraph::as_tbl_graph(
        igraph::graph_from_literal(A -+ B1, B2 -+ C))
    .f <- function(x) ifelse(x %in% c("B1", "B2"), "B", x)
    
    expect_true(gcm_cluster_check_admissible(g, .f)$admissible)
    
    d1 <- gcm_cluster_L1_degradation(g, .f)
    expect_equal(d1$d_degradation, 0.5, tolerance = tol)
    expect_equal(d1$n_lost, 1L)
    expect_match(d1$lost_claims, "A.*C")          # the marginal claim
    expect_false(any(d1$audit$status == "VIOLATION"))
    
    # L2 signature: no identifiability is destroyed (ErrLostId = 0),
    # but the quotient INVENTS a needed adjustment for (C, A) -- fine
    # says empty recipe (null effect, no paths), quotient demands {B}.
    # Pushforward fine->coarse is clean; the discrepancy lives in the
    # coarse-as-claimant direction.
    d2 <- gcm_cluster_L2_degradation(g, .f)
    expect_equal(d2$ErrLostId, 0, tolerance = tol)
    expect_equal(d2$d_degradation, 1/6, tolerance = tol)
    expect_gt(d2$d_degradation_pullback, 0)
})

# -------------------------------------------------------------------
# 6. Admissibility gate: the mediator-confounder merge announces itself
# -------------------------------------------------------------------
# Fine: X -> M -> Y ; C -> X ; C -> Y.  Merging {M, C} induces the
# 2-cycle X <-> MC (X -> M internally forward, C -> X internally
# backward). The classically "lossy" merge is caught at the gate, not
# in the metrics.

test_that("mediator-confounder merge is inadmissible with a 2-cycle", {
    g  <- tidygraph::as_tbl_graph(
        igraph::graph_from_literal(X -+ M, M -+ Y, C -+ X, C -+ Y))
    .f <- function(x) ifelse(x %in% c("M", "C"), "MC", x)
    
    adm <- gcm_cluster_check_admissible(g, .f)
    expect_false(adm$admissible)
    expect_gt(nrow(adm$two_cycles), 0)
    expect_true(any(
        (adm$two_cycles$a == "MC" & adm$two_cycles$b == "X") |
            (adm$two_cycles$a == "X"  & adm$two_cycles$b == "MC")))
})

# -------------------------------------------------------------------
# 7. Input guards: malformed maps fail loudly, not silently
# -------------------------------------------------------------------

test_that("a map producing empty-string names errors legibly", {
    g  <- tidygraph::as_tbl_graph(igraph::graph_from_literal(A -+ B))
    .f <- function(x) stringr::str_replace(x, "^B$", "")   # eats B
    expect_error(gcm_cluster_L1_degradation(g, .f))
    expect_error(gcm_cluster_L2_degradation(g, .f))
})
