\name{textcat_options}
\alias{textcat_options}
\title{Textcat Options}
\description{
  Get and set options used for n-gram based text categorization.
}
\usage{
textcat_options(option, value)
}
\arguments{
  \item{option}{character string indicating the option to get or set
    (see \bold{Details}).  If missing, all options are returned as a
    list.}
  \item{value}{Value to be set.  If omitted, the current value of the
    given option is returned.}
}
\details{
  Currently, the following options are available:
  \describe{
    \item{\code{n}:}{
      the maximum number of character in the n-gram profiles.
      
      Default: \code{5L}.
    }
    \item{\code{split}:}{
      the regular expression pattern to be used in word splitting.
      
      Default: \code{"[[:space:][:punct:][:digit:]]+"}.
    }
    \item{\code{tolower}:}{
      A logical indicating whether to transform texts to lower case
      (after word splitting).

      Default: \code{TRUE}.
    }
    \item{\code{reduce}:}{
      A logical indicating whether a representation of n-grams more
      efficient than the one used by Cavnar and Trenkle should be
      employed.
      
      Default: \code{TRUE}.
    }
    \item{\code{useBytes}:}{
      A logical indicating whether to use byte n-grams rather than
      character n-grams.
      
      Default: \code{FALSE}.
    }
    \item{\code{ignore}:}{
      a character vector of n-grams to be ignored when computing n-gram
      profiles.
      
      Default: \code{"_"} (corresponding to a word boundary).
    }
    \item{\code{size}:}{
      The maximal number of n-grams used for a profile.
      
      Default: \code{1000L}.
    }
    \item{\code{method}:}{
      A character string or function specifying a method for computing
      distances between n-gram profiles (see \code{\link{textcat}}).
      
      Default: \code{"CT"}, giving the Cavnar-Trenkle out of place
      measure.
    }
  }
}
\seealso{
  \code{\link{textcat_profile_db}} for how the first 6 options are used
  when computing n-gram profiles.
  
  \code{\link[tau]{textcnt}} in package \pkg{tau} which provides the
  functionality for term or pattern counting of text documents employed
  by \pkg{textcat}.
}
