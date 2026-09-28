#' HIV1 substitution model fitted via IQ-TREE2 (TVM+F+R5) from real HIV1
#' CRF01_AE sequences (Shenzhen 2006-2015; Yang et al. 2021,
#' https://doi.org/10.1093/ve/veab094). Default model for simulate.sequences().
#' @export
HIV1.SUBSTITUTION.MODEL <- list(
  gtr.rates = c(1.50852, 8.64133, 0.65614, 0.51625, 8.64133, 1.00000), # A-C,A-G,A-T,C-G,C-T,G-T
  bf        = c(0.390, 0.172, 0.218, 0.221),                          # A, C, G, T
  site.rate.categories = list(
    list(proportion = 0.602, rate = 0.031),
    list(proportion = 0.170, rate = 0.882),
    list(proportion = 0.157, rate = 2.185),
    list(proportion = 0.053, rate = 5.258),
    list(proportion = 0.018, rate = 11.941)
  )
)

#' simulate.sequences
#'
#' Simulate DNA sequences forward along resolve.arg()'s per-segment trees
#' under a substitution model with among-site rate variation, then stitch
#' the segments into one full sequence per tip. Requires phangorn.
#'
#' Strict clock only (all branches scaled by mu); no tip dates exported.
#'
#' @param resolved  output of resolve.arg()
#' @param mu        clock rate, substitutions/site per unit time
#' @param model     "GTR" (default), "HKY", or "JC"
#' @param gtr.rates length-6 relative rates (A-C,A-G,A-T,C-G,C-T,G-T), GTR only
#' @param bf        length-4 base frequencies (A,C,G,T); ignored for JC
#' @param kappa     transition/transversion ratio, HKY only
#' @param site.rate.categories  list of list(proportion=, rate=) for
#'                  among-site rate variation, must sum to 1. Defaults to
#'                  HIV1's fitted FreeRate profile; pass
#'                  list(list(proportion=1, rate=1)) to disable.
#'
#' @return named character vector, one sequence per tip
#' @export
simulate.sequences <- function(resolved, mu = 1e-3,
                                model = c("GTR", "HKY", "JC"),
                                gtr.rates = HIV1.SUBSTITUTION.MODEL$gtr.rates,
                                bf = HIV1.SUBSTITUTION.MODEL$bf,
                                kappa = 1,
                                site.rate.categories = HIV1.SUBSTITUTION.MODEL$site.rate.categories) {
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    stop("simulate.sequences() requires the 'phangorn' package. ",
         "Install it with install.packages(\"phangorn\").")
  }

  model <- match.arg(model)
  Q <- switch(model,
    GTR = gtr.rates,
    HKY = c(1, kappa, 1, 1, kappa, 1),
    JC  = rep(1, 6)
  )
  if (model == "JC") bf <- rep(0.25, 4)

  if (length(Q) != 6 || length(bf) != 4) {
    stop("simulate.sequences: gtr.rates must have length 6, bf length 4")
  }

  props <- vapply(site.rate.categories, `[[`, "proportion", FUN.VALUE = numeric(1))
  if (abs(sum(props) - 1) > 1e-6) {
    stop("simulate.sequences: site.rate.categories proportions must sum to 1 (got ",
         sum(props), ")")
  }

  category.lengths <- function(seg.len) {
    lens <- floor(seg.len * props)
    remainder <- seg.len - sum(lens)
    if (remainder > 0) {
      order.idx <- order(-props)
      lens[order.idx[seq_len(remainder)]] <- lens[order.idx[seq_len(remainder)]] + 1
    }
    lens
  }

  tip.seqs <- list()

  for (seg in resolved$local.trees) {
    if (is.null(seg$phylo)) {
      stop("simulate.sequences: segment ", seg$start, "-", seg$end,
           " has no valid tree (resolve.arg failed to parse it)")
    }

    tree <- seg$phylo
    tree$edge.length <- tree$edge.length * mu   # time -> substitutions/site
    seg.len <- seg$end - seg$start + 1L
    cat.lens <- category.lengths(seg.len)

    cat.sims <- list()
    for (k in seq_along(site.rate.categories)) {
      if (cat.lens[k] == 0) next
      cat.sims[[length(cat.sims) + 1]] <- phangorn::simSeq(
        tree, l = cat.lens[k], Q = Q, bf = bf, type = "DNA",
        rate = site.rate.categories[[k]]$rate
      )
    }
    sim   <- if (length(cat.sims) == 1) cat.sims[[1]] else do.call(c, cat.sims)
    align <- as.character(sim)

    for (tip in tree$tip.label) {
      tip.seqs[[tip]] <- c(tip.seqs[[tip]], align[tip, ])
    }
  }

  vapply(tip.seqs, paste, collapse = "", FUN.VALUE = character(1))
}
