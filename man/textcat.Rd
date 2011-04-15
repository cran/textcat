\name{textcat}
\alias{textcat}
\title{N-Gram Based Text Categorization}
\description{
  Categorize texts by computing their n-gram profiles, and finding the
  closest category n-gram profile.
}
\usage{
textcat(x, p = TC_char_profiles, method = "CT")
}
\arguments{
  \item{x}{a character vector of texts, or an \R object which can be
    coerced to this using \code{as.character}.}
  \item{p}{a textcat profile db (see \code{\link{textcat_profile_db}}).}
  \item{method}{a character string specifying a built-in method, or a
    used-defined function for computing distances between n-gram
    profiles.  See \bold{Details} for available built-in methods.}
}
\details{
  For each given text, its n-gram profile is computed using the options
  in the category profile db.  Then, the distance between this profile
  and the category profiles is computed, and the text is categorized
  into the category of the closest profile (if this is not unique,
  \code{NA} is obtained).

  Currently, the following distance methods for n-gram profiles are
  available.
  \describe{
    \item{\code{"CT"}:}{the out-of-place measure of Cavnar and Trenkle.}
    \item{\code{"ranks"}:}{a variant of the Cavnar/Trenkle measure based
      on the aggregated absolute difference of the ranks of the combined
      n-grams in the two profiles.}
    \item{\code{"ALPD"}:}{the sum of the absolute differences in n-gram
      log frequencies.}
    \item{\code{"KLI"}:}{the Kullback-Leibler I-divergence
      \eqn{I(p, q) = \sum_i p_i \log(p_i/q_i)}{I(p, q) = sum_i p_i log(p_i/q_i)}
      of the n-gram frequency distributions \eqn{p} and \eqn{q} of the
      two profiles.}
    \item{\code{"KLJ"}:}{the Kullback-Leibler J-divergence
      \eqn{J(p, q) = \sum_i (p_i - q_i) \log(p_i/q_i)}{J(p, q) = sum_i (p_i - q_i) log(p_i/q_i)},
      the symmetrized variant \eqn{I(p, q) + I(q, p)} of the I-divergences.
    }
    \item{\code{"JS"}:}{the Jensen-Shannon divergence between the n-gram
      frequency distributions.}
  }
  For the measures based on distances of frequency distributions,
  n-grams of the two profiles are combined, and missing n-grams are
  given a small positive absolute frequency (currently, 1e-6).

  Unless the profile db uses bytes rather than characters, the texts in
  \code{x} should be encoded in UTF-8.
}
\references{
  W. B. Cavnar and J. M. Trenkle (1994),
  N-Gram-Based Text Categorization.
  In ``Proceedings of SDAIR-94, 3rd Annual Symposium on Document
  Analysis and Information Retrieval'', 161--175.
}
\examples{
textcat(c("This is an english sentence.",
          "Das ist ein deutscher satz."))
}  
