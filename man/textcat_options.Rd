\name{textcat_options}
\alias{textcat_options}
\title{Textcat Options}
\description{
  Get and set options used for \eqn{n}-gram based text categorization.
}
\usage{
textcat_options(option, value)
}
\arguments{
  \item{option}{character string indicating the option to get or set
    (see \bold{Details}).  Can be abbreviated.  If missing, all options
    are returned as a list.}
  \item{value}{Value to be set.  If omitted, the current value of the
    given option is returned.}
}
\details{
  Currently, the following options are available:
  \describe{
    \item{\code{profile_method}:}{
      A character string or function specifying a method for computing
      \eqn{n}-gram profiles (see \code{\link{textcat_profile_db}}).

      Default: \code{"textcat"}.
    }
    \item{\code{profile_options}:}{
      A list of options to be passes to the method for computing
      profiles. 
    
      Default: none (empty list).
    }
    \item{\code{xdist_method}:}{
      A character string or function specifying a method for computing
      distances between \eqn{n}-gram profiles (see
      \code{\link{textcat_xdist}}).

      Default: \code{"CT"}, giving the Cavnar-Trenkle out of place
      measure.
    }
    \item{\code{xdist_options}:}{
      A list of options to be passes to the method for computing distances
      between profiles.
    
      Default: none (empty list).
    }
  }
}
