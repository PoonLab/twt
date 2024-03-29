% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{Lineage}
\alias{Lineage}
\title{Lineage}
\description{
\code{Lineage} is an R6 class for objects that represent pathogen lineages
that are carried by Compartments and which comprise the "inner" tree of the 
simulation.
}
\examples{
# load Compartments from a YAML object
path <- system.file('extdata', 'SI.yaml', package='twt')
settings <- yaml.load_file(path)
mod <- MODEL$new(settings)

# display first Lineage in first Compartment
comp <- mod$get.compartments()[[1]]
comp$get.lineages()  # display all 3 lineages

# manually add an unsampled Lineage
lin <- Lineage$new(name="L0", location=comp)
comp$add.lineage(lin)


}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Lineage$new()}}
\item \href{#method-copy}{\code{Lineage$copy()}}
\item \href{#method-get.name}{\code{Lineage$get.name()}}
\item \href{#method-get.type}{\code{Lineage$get.type()}}
\item \href{#method-get.sampling.time}{\code{Lineage$get.sampling.time()}}
\item \href{#method-get.location}{\code{Lineage$get.location()}}
\item \href{#method-set.location}{\code{Lineage$set.location()}}
\item \href{#method-set.location.by.name}{\code{Lineage$set.location.by.name()}}
\item \href{#method-clone}{\code{Lineage$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$new(name = NA, type = NA, sampling.time = NA, location = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-copy"></a>}}
\if{latex}{\out{\hypertarget{method-copy}{}}}
\subsection{Method \code{copy()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$copy(deep = FALSE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.name"></a>}}
\if{latex}{\out{\hypertarget{method-get.name}{}}}
\subsection{Method \code{get.name()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$get.name()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.type"></a>}}
\if{latex}{\out{\hypertarget{method-get.type}{}}}
\subsection{Method \code{get.type()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$get.type()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.sampling.time"></a>}}
\if{latex}{\out{\hypertarget{method-get.sampling.time}{}}}
\subsection{Method \code{get.sampling.time()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$get.sampling.time()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.location"></a>}}
\if{latex}{\out{\hypertarget{method-get.location}{}}}
\subsection{Method \code{get.location()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$get.location()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set.location"></a>}}
\if{latex}{\out{\hypertarget{method-set.location}{}}}
\subsection{Method \code{set.location()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$set.location(comp)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set.location.by.name"></a>}}
\if{latex}{\out{\hypertarget{method-set.location.by.name}{}}}
\subsection{Method \code{set.location.by.name()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$set.location.by.name(locationList, new.locationName)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Lineage$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
