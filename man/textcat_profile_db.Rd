\name{textcat_profile_db}
\alias{textcat_profile_db}
\title{Textcat Profile Dbs}
\description{
  Create \eqn{n}-gram profile dbs for text categorization.
}
\usage{
textcat_profile_db(x, id = NULL, method = NULL, ..., options = list(),
                   profiles = NULL)
}
\arguments{
  \item{x}{a character vector of text documents, or an \R object of text
    documents extractable via \code{as.character}.
  }
  \item{id}{a character vector giving the categories of the texts to be 
    recycled to the length of \code{x}, or \code{NULL} (default),
    indicating to treat each text document separately.
  }
  \item{method}{a character string specifying a built-in method, or a
    user-defined function for computing distances between \eqn{n}-gram
    profiles, or \code{NULL} (default), corresponding to using the
    method and options used for creating \code{profiles} if this is not
    \code{NULL}, or otherwise the current value of \pkg{textcat} option
    \code{profile_method} (see \code{\link{textcat_options}}).
  }
  \item{...}{options to be passed to the method for creating profiles.}
  \item{options}{a list of such options.}
  \item{profiles}{a textcat profile db object.}
}
\details{
  The text documents are split according to the given categories, and
  \eqn{n}-gram profiles are computed using the specified method, with
  options either those used for creating \code{profiles} if this is not
  \code{NULL}, or by combining the options given in \code{...} and
  \code{options} and merging with the default profile options specified
  by the \pkg{textcat} option \code{profile_options} using \emph{exact}
  name matching.  The method and options employed for building the db
  are stored in the db as attributes \code{"method"} and
  \code{"options"}, respectively.

  There is a \code{\link{c}} method for combining profile dbs provided
  that these have identical options.  There are also a \code{[} method
  for subscripting and \code{\link{as.matrix}} and
  \code{\link[slam:matrix]{as.simple.triplet.matrix}} methods to
  \dQuote{export} the profiles to a dense matrix or the sparse simple
  triplet matrix representation provided by package \pkg{slam},
  respectively.

  Currently, the only available built-in method is \code{"textcnt"},
  which has the following options:
  \describe{
    \item{\code{n}:}{
      A numeric vector giving the numbers of characters or bytes in the
      \eqn{n}-gram profiles.
      
      Default: \code{1 : 5}.
    }
    \item{\code{split}:}{
      The regular expression pattern to be used in word splitting.
      
      Default: \code{"[[:space:][:punct:][:digit:]]+"}.
    }
    \item{\code{perl}:}{
      A logical indicating whether to use Perl-compatible regular
      expressions in word splitting.

      Default: \code{FALSE}.
    }
    \item{\code{tolower}:}{
      A logical indicating whether to transform texts to lower case
      (after word splitting).

      Default: \code{TRUE}.
    }
    \item{\code{reduce}:}{
      A logical indicating whether a representation of \eqn{n}-grams
      more efficient than the one used by Cavnar and Trenkle should be
      employed.
      
      Default: \code{TRUE}.
    }
    \item{\code{useBytes}:}{
      A logical indicating whether to use byte \eqn{n}-grams rather than
      character \eqn{n}-grams.
      
      Default: \code{FALSE}.
    }
    \item{\code{ignore}:}{
      a character vector of \eqn{n}-grams to be ignored when computing
      \eqn{n}-gram profiles.
      
      Default: \code{"_"} (corresponding to a word boundary).
    }
    \item{\code{size}:}{
      The maximal number of \eqn{n}-grams used for a profile.
      
      Default: \code{1000L}.
    }
  }
  This method uses \code{\link[tau]{textcnt}} in package \pkg{tau} for
  computing \eqn{n}-gram profiles, with \code{n}, \code{split},
  \code{perl} and \code{useBytes} corresponding to the respective
  \code{textcnt} arguments, and option \code{reduce} setting argument
  \code{marker} as needed.  \eqn{N}-grams listed in option \code{ignore}
  are removed, and only the most frequent remaining ones retained, with
  the maximal number given by option \code{size}.  

  Unless the profile db uses bytes rather than characters (i.e., option
  \code{useBytes} is \code{TRUE}), text documents in \code{x} containing
  non-ASCII characters must declare their encoding (see
  \code{\link{Encoding}}), and will be re-encoded to UTF-8.

  Note that option \code{n} specifies \emph{all} numbers of characters
  or bytes to be used in the profiles, and not just the maximal number:
  e.g., taking \code{n = 3} will create profiles only containing
  tri-grams.
}
\examples{
## Obtain the texts of the standard licenses shipped with R.
files <- dir(file.path(R.home("share"), "licenses"), "^[A-Z]",
             full.names = TRUE)
texts <- sapply(files,
                function(f) paste(readLines(f), collapse = "\n"))
names(texts) <- basename(files)
## Build a profile db using the same method and options as for building
## the ECIMCI character profiles.
profiles <- textcat_profile_db(texts, profiles = ECIMCI_profiles)
## Inspect the 10 most frequent n-grams in each profile.
lapply(profiles, head, 10L)
## Combine into one frequency table.
tab <- as.matrix(profiles)
tab[, 1 : 10]
## Determine languages.
textcat(profiles, ECIMCI_profiles)
}
