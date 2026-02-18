\name{textcat}
\alias{textcat}
\title{\eqn{N}-Gram Based Text Categorization}
\description{
  Categorize texts by computing their \eqn{n}-gram profiles, and finding
  the closest category \eqn{n}-gram profile.
}
\usage{
textcat(x, p = textcat::TC_char_profiles, method = "CT", ...,
        options = list())
}
\arguments{
  \item{x}{a character vector of texts, or an \R object which can be
    coerced to this using \code{as.character}, or a textcat profile db
    (see \code{\link{textcat_profile_db}}) created using the same method
    and options as \code{p}.}
  \item{p}{a textcat profile db.  By default, the TextCat character
    profiles are used (see \code{\link{TC_char_profiles}}).}
  \item{method}{a character string specifying a built-in method, or a
    user-defined function for computing distances between \eqn{n}-gram
    profiles.  See \code{\link{textcat_xdist}} for details.}
  \item{...}{options to be passed to the method for computing distances
    between profiles.}
  \item{options}{a list of such options.}
}
\details{
  For each given text, its \eqn{n}-gram profile is computed using the
  options in the category profile db.  Then, the distance between this
  profile and the category profiles is computed, and the text is
  categorized into the category of the closest profile (if this is not
  unique, \code{NA} is obtained).

  Unless the profile db uses bytes rather than characters, the texts in
  \code{x} should be encoded in UTF-8.
}
\references{
  \bibshow{textcat::Cavnar+Trenkle:1994,
    textcat::Hornik+Mair+Rauch:2013}
}
\examples{
textcat(c("This is an english sentence.",
          "Das ist ein deutscher satz."))
}  
