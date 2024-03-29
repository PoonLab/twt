% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model.R
\name{Model}
\alias{Model}
\title{Model}
\description{
\code{Model} is an R6 class that defines an object that generates all of the 
objects of the various classes (e.g., Compartment, CompartmentType) that 
define a simulation model.  A \code{Model} object is immutable - it should
not change over the course of simulation.  Instead, we derive a \code{Run}
object class that inherits from a \code{Model}.
}
\examples{
require(twt)
# get path to example YAML file
path <- system.file('extdata', 'SI.yaml', package='twt')

# load file and parse to construct Model object
settings <- yaml.load_file(path)
mod <- Model$new(settings)

# display summary information (calls S3 print method)
mod

}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{initial.conds}}{Initial Conditions}

\item{\code{types}}{vector of CompartmentType objects}

\item{\code{compartments}}{vector of Compartment objects}

\item{\code{lineages}}{vector of Lineage objects}

\item{\code{fixed.samplings}}{list of Lineage names and sampling times for plotting}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{initial.conds}}{Initial Conditions}

\item{\code{types}}{vector of CompartmentType objects}

\item{\code{compartments}}{vector of Compartment objects}

\item{\code{lineages}}{vector of Lineage objects}

\item{\code{fixed.samplings}}{list of Lineage names and sampling times for plotting}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Model$new()}}
\item \href{#method-get.initial.conds}{\code{Model$get.initial.conds()}}
\item \href{#method-get.types}{\code{Model$get.types()}}
\item \href{#method-get.compartments}{\code{Model$get.compartments()}}
\item \href{#method-get.lineages}{\code{Model$get.lineages()}}
\item \href{#method-get.names}{\code{Model$get.names()}}
\item \href{#method-get.fixed.samplings}{\code{Model$get.fixed.samplings()}}
\item \href{#method-clone}{\code{Model$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$new(settings = NA, name = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.initial.conds"></a>}}
\if{latex}{\out{\hypertarget{method-get.initial.conds}{}}}
\subsection{Method \code{get.initial.conds()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$get.initial.conds()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.types"></a>}}
\if{latex}{\out{\hypertarget{method-get.types}{}}}
\subsection{Method \code{get.types()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$get.types()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.compartments"></a>}}
\if{latex}{\out{\hypertarget{method-get.compartments}{}}}
\subsection{Method \code{get.compartments()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$get.compartments()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.lineages"></a>}}
\if{latex}{\out{\hypertarget{method-get.lineages}{}}}
\subsection{Method \code{get.lineages()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$get.lineages()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.names"></a>}}
\if{latex}{\out{\hypertarget{method-get.names}{}}}
\subsection{Method \code{get.names()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$get.names(listR6obj)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get.fixed.samplings"></a>}}
\if{latex}{\out{\hypertarget{method-get.fixed.samplings}{}}}
\subsection{Method \code{get.fixed.samplings()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$get.fixed.samplings()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
