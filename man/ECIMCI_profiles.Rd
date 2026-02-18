\name{ECIMCI_profiles}
\alias{ECIMCI_profiles}
\title{ECI/MCI \eqn{N}-Gram Profiles}
\description{
  \eqn{N}-gram profile db for 26 languages based on the European Corpus
  Initiative Multilingual Corpus I.
}
\usage{ECIMCI_profiles}
\details{
  This profile db was built by Johannes Rauch, using the ECI/MCI corpus
  (\url{http://www.elsnet.org/eci.html}) and the default options
  employed by package \pkg{textcat}, with all text documents encoded in
  UTF-8.

  The category ids used for the db are the respective IETF language tags
  (see \link[NLP]{parse_IETF_language_tag} in package \pkg{NLP}), using the ISO 639-2
  Part B language subtags and, for Serbian, the script employed (i.e., 
  \code{"scc-Cyrl"} and \code{"scc-Latn"} for Serbian written in
  Cyrillic and Latin script, respectively; all other languages in the
  profile db are written in Latin script.)
}
\references{
  \bibshow{textcat::Armstrong-Warwick+Thompson+Mckelvie:1994}
}
\examples{
## Languages in the the ECI/MCI profile db:
names(ECIMCI_profiles)
## Key options used for the profile:
attr(ECIMCI_profiles, "options")[c("n", "size", "reduce", "useBytes")]
}
\keyword{datasets}
