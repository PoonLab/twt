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
                                site.rate.categories = HIV1.SUBSTITUTION.MODEL$site.rate.categories,
                                seed = NULL) {
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

  if (!is.null(seed)) set.seed(seed)

  segs <- resolved$local.trees
  seg.starts <- vapply(segs, function(s) as.numeric(s$start), numeric(1))
  seg.ends   <- vapply(segs, function(s) as.numeric(s$end), numeric(1))
  genome.len <- as.integer(max(seg.ends))
  if (sum(seg.ends - seg.starts + 1) != genome.len) {
    stop("simulate.sequences: segments do not tile sites 1-", genome.len,
         " exactly (total segment length ", sum(seg.ends - seg.starts + 1), ")")
  }

  # A site's rate category is a property of its genomic position, so draw it
  # once for the whole genome (exact proportions, random positions) and reuse
  # it for every segment, whichever local tree covers the site.
  site.cat <- sample(rep(seq_along(site.rate.categories),
                         times = category.lengths(genome.len)))

  tip.seqs <- list()

  for (seg in segs) {
    if (is.null(seg$phylo)) {
      stop("simulate.sequences: segment ", seg$start, "-", seg$end,
           " has no valid tree (resolve.arg failed to parse it)")
    }

    tree <- seg$phylo
    tree$edge.length <- tree$edge.length * mu   # time -> substitutions/site
    cols <- seq.int(seg$start, seg$end)
    seg.cat <- site.cat[cols]
    seg.align <- matrix(NA_character_, nrow = length(tree$tip.label), ncol = length(cols),
                        dimnames = list(tree$tip.label, NULL))

    for (k in seq_along(site.rate.categories)) {
      sel <- which(seg.cat == k)
      if (length(sel) == 0) next
      sim <- phangorn::simSeq(
        tree, l = length(sel), Q = Q, bf = bf, type = "DNA",
        rate = site.rate.categories[[k]]$rate
      )
      a <- as.character(sim)
      seg.align[rownames(a), sel] <- a
    }

    for (tip in tree$tip.label) {
      tip.seqs[[tip]] <- c(tip.seqs[[tip]], seg.align[tip, ])
    }
  }

  out <- vapply(tip.seqs, paste, collapse = "", FUN.VALUE = character(1))
  if (any(nchar(out) != genome.len)) {
    stop("simulate.sequences: not every tip appears in every segment tree ",
         "(expected ", genome.len, " sites per tip)")
  }
  out
}
