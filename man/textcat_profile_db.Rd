\name{textcat_profile_db}
\alias{textcat_profile_db}
\title{Textcat Profile Dbs}
\description{
  Create n-gram profile dbs for text categorization.
}
\usage{
textcat_profile_db(x, id, ...)
}
\arguments{
  \item{x}{a character vector of text documents, or an \R object of text
    documents extractable via \code{as.character}.
  }
  \item{id}{a character vector giving the categories of the texts.
    Recycled to the length of \code{x}.
  }
  \item{...}{further arguments specifying the options used for creating
    the n-gram profiles, see \code{\link{textcat_options}} for the
    (current) default options.  The names of the arguments are partially
    matched against the names of the defaults, and used for the options
    instead in case of unique matches.
  }
}
\details{
  The text documents are split according to the given categories, and
  n-gram profiles are computed via \code{\link[tau]{textcnt}} in package
  \pkg{tau}, with \pkg{textcat} options \code{n}, \code{split},
  \code{perl} and \code{useBytes} corresponding to the respective
  \code{textcnt} arguments, and option \code{reduce} setting argument
  \code{marker} as needed.  N-grams listed in option \code{ignore} are
  removed, and only the most frequent remaining ones retained, with the
  maximal number given by option \code{size}.  The options employed for
  building the db are stored in the db.

  There is a \code{\link{c}} method for combining profile dbs provided
  that these have identical options.

  Unless the profile db uses bytes rather than characters (i.e., option
  \code{useBytes} is \code{TRUE}), the text documents in \code{x} should
  be encoded in UTF-8.
}
