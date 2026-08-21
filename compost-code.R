# arxode

# adj_fingerprint 
# getting the ZX, idX etc from function below

# Needs to work better for M5, D5, D7, ...
# Many NaNs

gcm_fuzzy_L2_consistency <- function(gX, gY,
                                     adj_type = "canonical") {
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
    # The distance is then weighted according to the number of nodes shared
    # An exception to the distance is returning the maximum value (1) if
    # gX deems the nodes identifiable and gY thinks they are not, and if
    # both models agree the node pair is unidentifiable then the min distance 
    # (0) is returned. 
    if (!(class(gX)[[1]] == "dagitty")) gX <- tidygraph_to_dagitty(gX)
    if (!(class(gY)[[1]] == "dagitty")) gY <- tidygraph_to_dagitty(gY)
    
    shared <- intersect(names(gX), names(gY))
    
    pairs <- expand_grid(x = shared, y = shared) |> 
        filter(x != y)
    
    canAdjSet <- function(g, exposure, outcome) {
        if (str_detect(adj_type, "^[Cc]")) {
            s <- adjustmentSets(g, exposure = exposure, outcome = outcome, 
                                type = adj_type)
        } else if (str_detect(adj_type, "^[Oo]")) {
            s <- gcm_optimal_adjustment_set(g, exposure, outcome)
        } else {
            warning(str_c(adj_type, ": adjustment method not recognised."))
            return(NULL)
        }
        if (length(s) == 0) return(NULL)              # not identifiable
        sort(as.character(unname(unlist(s))))         # identifiable (possibly empty)
    }
    
    df.out <- pairs |>
        rowwise() |>
        mutate(
            zx = list(canAdjSet(gX, x, y)),
            zy = list(canAdjSet(gY, x, y)),
            idX = !is.null(zx),
            idY = !is.null(zy),
            idXY = !(idX & !idY),     # identifiable in X ⇒ identifiable in Y
            idYX = !(idY & !idX),
            zx_sh = list(intersect(zx, shared)),   
            zy_sh = list(intersect(zy, shared)),
            # Err on the adjustment sets
            ErrXY = length(setdiff(zx_sh, zy_sh)),
            ErrYX = length(setdiff(zy_sh, zx_sh)),
            ZXsharedMax = length(zx_sh),
            ZYsharedMax = length(zy_sh)
        ) |>
        mutate(
            # per-pair X->Y distance
            dXY = case_when(
                # identifiability violation: max distance
                idX & !idY            ~ 1,   
                # agree not identifiable: min distance (and no adj sets)
                !idX & !idY           ~ 0,   
                # graded set error
                idX &  idY            ~ ErrXY / pmax(ZXsharedMax, 1),  
                # antecedent false: excluded. X => Y uninformative for X => !Y
                TRUE                  ~ NA_real_                   
            ),
            dYX = case_when(
                idY & !idX            ~ 1,
                !idY & !idX           ~ 0,
                idY &  idX            ~ ErrYX / pmax(ZYsharedMax, 1),
                TRUE                  ~ NA_real_
            )
        ) |> 
        ungroup() |> 
        mutate(
            # weight: how much comparable evidence underlies each pair's score
            w_XY = case_when(
                # excluded, won't enter mean (this is when X => !Y)
                is.na(dXY)                          ~ 0,            
                # VACUOUS zero: not comparable. This is where the denominator
                # is forced to be 1 in dXY above
                idX & idY & ZXsharedMax == 0        ~ 0,            
                # genuine mutual-agreement on non-identifiability (a strong statement)
                !idX & !idY                         ~ 1,            
                # weight by comparable nodes using the set error from dXY
                TRUE                                ~ ZXsharedMax   
            ),
            w_YX = case_when(
                is.na(dYX)                          ~ 0,            
                idY & idX & ZYsharedMax == 0        ~ 0,            
                !idY & !idX                         ~ 1,            
                TRUE                                ~ ZYsharedMax   
            )
        ) |> 
        select(x, y, dXY, dYX, everything())
    
    dXY <- weighted.mean(df.out$dXY, w = df.out$w_XY, na.rm = TRUE)
    dYX <- weighted.mean(df.out$dYX, w = df.out$w_YX, na.rm = TRUE)
    
    L2XY <- 1 - dXY
    L2YX <- 1 - dYX
    
    # idXY <-  mean(df.out$idXY)
    # idYX <-  mean(df.out$idYX)
    # if (sum(df.out$ZXsharedMax) > 0) {
    #     ErrXY <- sum(df.out$ErrXY) / sum(df.out$ZXsharedMax)
    # } else {
    #     # mean of id agreement
    #     ErrXY <- mean(df.out$idXY == df.out$idYX)
    # }
    # if (sum(df.out$ZYsharedMax) > 0) {
    #     ErrYX <- sum(df.out$ErrYX) / sum(df.out$ZYsharedMax)
    # } else {
    #     ErrYX <- mean(df.out$idYX == df.out$idXY)
    # }
    
    return(list(
        # causal consistency (asymmetric)
        L2XY = L2XY, L2YX = L2YX,
        # L2-consistency distance between models (asymmetric)
        dXY = dXY, dYX = dYX,
        # artefacts for diagnosis
        df = df.out, shared_nodes = shared
    ))
    
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
                n_int == 0 & n_uni == 0 & indX != indY~ 1, 
                # overlapping sets, use jaccard distance
                TRUE ~ 1 - jaccard
            )
        ) |> 
        select(-.a, -.b) 
    
    if (nrow(imp_ci) == 0) {
        # message("No shared independence relationships")
        return(
            list(
                distance = 1,
                L1_consistency = 0,
                ci_compasisons = imp_ci
            )
        )
    }
    
    distance <- mean(imp_ci$distance, na.rm = T)
    L1_consistency <- 1 - distance
    # no longer used
    L1_hard <- 1 - mean(imp_ci$distance == 0, na.rm = T)
    
    return(
        list(
            distance = distance,
            L1_consistency = L1_consistency,
            ci_comparisions = imp_ci)
    )
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

# ---- Lifted fingerprint: fine graph at cluster resolution -----------

gcm_cluster_adj_fingerprint <- function(
        g, .f, pairs = NULL, map_sets = TRUE,
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


gcm_claim_keys <- function(fp) {
    if (nrow(fp) == 0) return(character(0))
    a  <- pmin(fp$x, fp$y)          # elementwise lexicographic min/max
    b  <- pmax(fp$x, fp$y)
    Zc <- map_chr(strsplit(fp$Z, ","), 
                  \(z) paste(sort(z, method = "radix"), collapse = ","))
    paste(a, "_||_", b, "|", Zc)
}



## L0 consistency
# Not really a thing on the causal hierarchy, but instead just a measure
# of how many nodes (relatively) the two models share
# Here the "fingerprint" is simply the list of nodes
gcm_L0_consistency <- function(gX, gY, kappa = 1) {
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


# Getting L1 level conditions

gcm_implied_conditional_independencies <- function(g) {
    if ("tbl_graph" %in% class(g)) g <- tidygraph_to_dagitty(g)
    g |> 
        dagitty::impliedConditionalIndependencies() |> 
        map(\(x) tibble(X = x$X, Y = x$Y, Z = list(x$Z))) |> 
        list_rbind() |> 
        arrange(X, Y)
}


