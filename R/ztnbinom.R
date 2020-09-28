## ztnbinom: Zero-truncated negative binomial

# from https://rdrr.io/rforge/countreg/src/R/ztnbinom.R
# maintainer: Achim Zeileis <Achim.Zeileis@R-project.org>
# released under GPL-2 / GPL-3

dztnbinom <- function(x, mu, theta, size, log = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- dnbinom(x, mu = mu, size = theta, log = TRUE) - pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[x < 1] <- -Inf
  rval[mu <= 0] <- 0
  if(log) rval else exp(rval)
}

pztnbinom <- function(q, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- log(pnbinom(q, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE) - dnbinom(0, mu = mu, size = theta)) -
    pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  rval[q < 1] <- if(lower.tail) -Inf else 0
  if(log.p) rval else exp(rval)
}

qztnbinom <- function(p, mu, theta, size, lower.tail = TRUE, log.p = FALSE) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  p_orig <- p
  p <- if(log.p) p else log(p)
  p <- p + pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
  p <- exp(p) + dnbinom(0, mu = mu, size = theta)
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = lower.tail, log.p = FALSE)
  if(lower.tail) rval[p_orig < dztnbinom(1, mu = mu, theta = theta, log = log.p)] <- 1
  rval
}

rztnbinom <- function(n, mu, theta, size) {
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  qztnbinom(runif(n), mu = mu, theta = theta)
}
