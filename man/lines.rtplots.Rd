% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{lines.rtplots}
\alias{lines.rtplots}
\title{Lines Method for 'rtplots' Object}
\usage{
\method{lines}{rtplots}(object, val = NULL, ind = 1, flip = F, ...)
}
\arguments{
\item{object}{a 'rtplots' object.}

\item{val}{the choice or accuracy value whose observations
are to be plotted.}

\item{ind}{a index to select observations for a specific
choice/accuracy value.}

\item{flip}{logical; if true, flips the y-axis values to
be negative.}

\item{...}{additional plotting parameters. See
\code{\link[graphics]{lines}}.}
}
\description{
Provides a method to 'rtplots' objects that allows
for drawing line segments on an existing plot.
}
