#' plot.OuterTree.with.dynamics
#'
#' Plots the outer transmission tree (left panel) alongside the epidemic
#' trajectory I(t) (right panel), with superinfection event times marked
#' as vertical dashed lines on both panels.
#'
#' @param ot       R6 object of class OuterTree
#' @param dynamics S3 object returned by sim.dynamics (with counted=TRUE, 
#'                 the default)
#' @param i.comp   character, name of the infected compartment to plot 
#'                 (default: "I")
#' @param pad      numeric, x-axis padding for tree panel (default: 1.05)
#' @export
plot.OuterTree.with.dynamics <- function(ot, dynamics, i.comp="I", pad=1.05) {

  # extract superinfection event times from outer log
  events <- ot$get.log()
  events$time <- as.numeric(events$time)
  trans <- events[events$event == "transmission", ]
  trans <- trans[order(trans$time), ]
  si.times <- trans$time[duplicated(trans$to.host)]

  # extract I(t) from dynamics
  ev <- dynamics$events
  # dynamics may not have compartment columns if counted=FALSE
  if (!i.comp %in% names(ev)) {
    stop("Compartment '", i.comp, "' not found in dynamics$events. ",
         "Was sim.dynamics() called with counted=TRUE?")
  }
  ev$time <- as.numeric(ev$time)
  # prepend t=simTime row at initial sizes
  mod <- dynamics$model
  init <- mod$get.init.sizes()
  t0 <- as.numeric(mod$get.parameters()$simTime)
  row0 <- as.list(setNames(rep(0, ncol(ev)), names(ev)))
  row0$time <- t0
  for (cn in names(init)) if (cn %in% names(ev)) row0[[cn]] <- init[[cn]]
  ev <- rbind(as.data.frame(row0), ev)

  # layout: tree left, I(t) right
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  layout(matrix(c(1,2), nrow=1), widths=c(3,2))

  # left: outer tree
  par(mar=c(5,1,2,1))
  withCallingHandlers(
    plot.OuterTree(ot, pad=pad),
    warning = function(w) {
      if (grepl("no non-missing arguments to max", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )
  title(main="Transmission tree", cex.main=0.9)

  # mark SI times as vertical lines on tree
  if (length(si.times) > 0) {
    abline(v=si.times, col=adjustcolor("steelblue", alpha.f=0.4),
           lty=2, lwd=1)
  }

  # right: I(t) with SI times
  par(mar=c(5,4,2,1))
  xlim <- range(ev$time)
  ylim <- c(0, max(ev[[i.comp]], na.rm=TRUE) * 1.1)

  plot(ev$time, ev[[i.comp]], type="s",
       xlim=rev(xlim),   # reverse so t=0 (root) aligns with tree on left
       ylim=ylim,
       xlab="Time", ylab=paste0("N(", i.comp, ")"),
       col="black", lwd=1.5, bty="l",
       main="Epidemic trajectory", cex.main=0.9)

  # SI event times as vertical lines
  if (length(si.times) > 0) {
    abline(v=si.times, col=adjustcolor("steelblue", alpha.f=0.5),
           lty=2, lwd=1.2)
    # rug on x-axis
    rug(si.times, col="steelblue", lwd=1.5, ticksize=0.04)
  }

  # legend
  legend("topleft", bty="n", cex=0.8,
         legend=c(sprintf("I(t)  [%s]", i.comp),
                  sprintf("SI events (n=%d)", length(si.times))),
         col=c("black", "steelblue"), lty=c(1,2), lwd=c(1.5,1.2))

  invisible(list(si.times=si.times, n.si=length(si.times)))
}
