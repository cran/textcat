\name{textcat_xdist}
\alias{textcat_xdist}
\title{Cross-Distances Between \eqn{N}-Gram Profiles}
\description{
  Compute cross-distances between collections of \eqn{n}-gram profiles.
}
\usage{
textcat_xdist(x, p = NULL, method = "CT", ..., options = list())
}
\arguments{
  \item{x}{a textcat profile db (see \code{\link{textcat_profile_db}}),
    or an \R object of text documents extractable via
    \code{as.character}.}
  \item{p}{\code{NULL} (default), or as for \code{x}.
    The default is equivalent to taking \code{p} as \code{x} (but more
    efficient).}
  \item{method}{a character string specifying a built-in method, or a
    used-defined function for computing distances between \eqn{n}-gram
    profiles, or \code{NULL} (corresponding to the current value of
    \pkg{textcat} option \code{xdist_method} (see
    \code{\link{textcat_options}}).
    See \bold{Details} for available built-in methods.}
  \item{...}{options to be passed to the method for computing
    distances.}
  \item{options}{a list of such options.}
}
\details{
  If \code{x} (or \code{p}) is not a profile db, the \eqn{n}-gram
  profiles of the individual text documents extracted from it are
  computed using the profile method and options in \code{p} if this is a
  profile db, and using the current \pkg{textcat} profile method and
  options otherwise. 

  Currently, the following distance methods for \eqn{n}-gram profiles
  are available.
  \describe{
    \item{\code{"CT"}:}{the out-of-place measure of Cavnar and Trenkle.}
    \item{\code{"ranks"}:}{a variant of the Cavnar/Trenkle measure based
      on the aggregated absolute difference of the ranks of the combined
      \eqn{n}-grams in the two profiles.}
    \item{\code{"ALPD"}:}{the sum of the absolute differences in
      \eqn{n}-gram log frequencies.}
    \item{\code{"KLI"}:}{the Kullback-Leibler I-divergence
      \eqn{I(p, q) = \sum_i p_i \log(p_i/q_i)}{I(p, q) = sum_i p_i log(p_i/q_i)}
      of the \eqn{n}-gram frequency distributions \eqn{p} and \eqn{q} of the
      two profiles.}
    \item{\code{"KLJ"}:}{the Kullback-Leibler J-divergence
      \eqn{J(p, q) = \sum_i (p_i - q_i) \log(p_i/q_i)}{J(p, q) = sum_i (p_i - q_i) log(p_i/q_i)},
      the symmetrized variant \eqn{I(p, q) + I(q, p)} of the I-divergences.
    }
    \item{\code{"JS"}:}{the Jensen-Shannon divergence between the
      \eqn{n}-gram frequency distributions.}
    \item{\code{"cosine"}}{the cosine dissimilarity between the
      profiles, i.e., one minus the inner product of the frequency
      vectors normalized to Euclidean length one (and filled with zeros
      for entries missing in one of the vectors).}
      
  }
  For the measures based on distances of frequency distributions,
  \eqn{n}-grams of the two profiles are combined, and missing
  \eqn{n}-grams are given a small positive absolute frequency which can
  be controlled by option \code{eps}, and defaults to 1e-6.

  Options given in \code{...} and \code{options} are combined, and
  merged with the default xdist options specified by the \pkg{textcat}
  option \code{xdist_options} using \emph{exact} name matching.
}
\examples{
## Compute cross-distances between the TextCat byte profiles using the
## CT out-of-place measure.
d <- textcat_xdist(TC_byte_profiles)
## Visualize results of hierarchical cluster analysis on the distances.
plot(hclust(as.dist(d)), cex = 0.7)
}
