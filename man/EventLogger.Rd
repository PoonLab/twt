% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EventLogger.R
\name{EventLogger}
\alias{EventLogger}
\title{EventLogger}
\description{
\code{EventLogger} is an R6 class for an object that tracks migration, 
transmission, and coalescent events.  Note that bottleneck events are logged as 
coalescent events.
}
\examples{
# manually initialize an EventLog object
e <- EventLogger$new()
# note this log entry is not linked to existing Lineage or Compartment objects
e$add.event("transmission", time=1, line1="NA", comp1="host1", comp2="host2")
e$get.all.events()

}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{EventLogger$new()}}
\item \href{#method-get.all.events}{\code{EventLogger$get.all.events()}}
\item \href{#method-get.events}{\code{EventLogger$get.events()}}
\item \href{#method-add.event}{\code{EventLogger$add.event()}}
\item \href{#method-clear.events}{\code{EventLogger$clear.events()}}
\item \href{#method-blank.events}{\code{EventLogger$blank.events()}}
\item \href{#method-record.transmission}{\code{EventLogger$record.transmission()}}
\item \href{#method-record.migration}{\code{EventLogger$record.migration()}}
\item \href{#method-get.fixed.samplings}{\code{EventLogger$get.fixed.samplings()}}
\item \href{#method-store.fixed.samplings}{\code{EventLogger$store.fixed.samplings()}}
\item \href{#method-clone}{\code{EventLogger$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$new(events = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.all.events"></a>}}
\if{latex}{\out{\hypertarget{method-get.all.events}{}}}
\subsection{Method \code{get.all.events()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$get.all.events()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.events"></a>}}
\if{latex}{\out{\hypertarget{method-get.events}{}}}
\subsection{Method \code{get.events()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$get.events(event.type)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-add.event"></a>}}
\if{latex}{\out{\hypertarget{method-add.event}{}}}
\subsection{Method \code{add.event()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$add.event(
  type,
  time,
  line1 = NA,
  line2 = NA,
  comp1 = NA,
  comp2 = NA,
  type1 = NA,
  type2 = NA
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type:}}{event type, one of 'transmission', 'transition', 'migration', 
'coalescence' or 'bottleneck'.}

\item{\code{time:}}{CUMULATIVE time that event has occurred between two compartments 
       in a transmission/migration/coalescent event
Record which Lineages are transmitted from source to 
recipient (should only be one entry with Compartment as recipient)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clear.events"></a>}}
\if{latex}{\out{\hypertarget{method-clear.events}{}}}
\subsection{Method \code{clear.events()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$clear.events()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-blank.events"></a>}}
\if{latex}{\out{\hypertarget{method-blank.events}{}}}
\subsection{Method \code{blank.events()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$blank.events()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-record.transmission"></a>}}
\if{latex}{\out{\hypertarget{method-record.transmission}{}}}
\subsection{Method \code{record.transmission()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$record.transmission(recipient, lineages)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{recipient:}}{Compartment object}

\item{\code{lineages:}}{a vector of names of Lineages to transfer out
                  of recipient Compartment
Record which Lineages are transmitted from source to recipient
through a migration event}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-record.migration"></a>}}
\if{latex}{\out{\hypertarget{method-record.migration}{}}}
\subsection{Method \code{record.migration()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$record.migration(recipient, source, time, lineages)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{recipient:}}{Compartment object}

\item{\code{source:}}{Compartment object}

\item{\code{time:}}{double, time of migration event}

\item{\code{lineages:}}{a list of Lineage objects}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.fixed.samplings"></a>}}
\if{latex}{\out{\hypertarget{method-get.fixed.samplings}{}}}
\subsection{Method \code{get.fixed.samplings()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$get.fixed.samplings()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-store.fixed.samplings"></a>}}
\if{latex}{\out{\hypertarget{method-store.fixed.samplings}{}}}
\subsection{Method \code{store.fixed.samplings()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$store.fixed.samplings(model.fixed.samplings)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{EventLogger$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
