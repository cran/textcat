\name{ECIMCI_profiles}
\alias{ECIMCI_profiles}
\title{ECI/MCI N-Gram Profiles}
\description{
  N-gram profile db for 26 languages based on the European Corpus
  Initiative Multilingual Corpus I.
}
\usage{ECIMCI_profiles}
\details{
  This profile db was built by Johannes Rauch, using the ECI/MCI corpus
  and the default options employed by package \pkg{textcat}, with all
  text documents encoded in UTF-8.

  The category ids used for the db are the respective IETF language tags
  (see \link[tau]{language} in package \pkg{tau}), using the ISO 639-2
  Part B language subtags and, for Serbian, the script employed (i.e., 
  \code{"scc-Cyrl"} and \code{"scc-Latn"} for Serbian written in
  Cyrillic and Latin script, respectively; all other languages in the
  profile are always written in Latin script.)
}
\keyword{datasets}
