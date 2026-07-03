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
