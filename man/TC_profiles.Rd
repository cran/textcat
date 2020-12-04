\name{TC_profiles}
\alias{TC_char_profiles}
\alias{TC_byte_profiles}
\title{TextCat \eqn{N}-Gram Profiles}
\description{
  TextCat \eqn{n}-gram byte and character profile dbs for language
  identification.
}
\usage{
TC_char_profiles
TC_byte_profiles
}
\details{
  TextCat (\url{https://www.let.rug.nl/vannoord/TextCat/}) is a Perl
  implementation of the Cavnar and Trenkle \dQuote{\eqn{N}-Gram-Based
    Text Categorization} technique by Gertjan van Noord which was
  subsequently integrated into SpamAssassin.  It provides byte
  \eqn{n}-gram profiles for 74 \dQuote{languages} (more precisely,
  language/encoding combinations).  The C library reimplementation
  \code{libtextcat} (\url{https://software.wise-guys.nl/libtextcat/})
  adds one more non-empty profile.

  \code{TC_byte_profiles} provides these byte profiles.

  \code{TC_char_profiles} provides a subset of 56 character profiles
  obtained by converting the byte sequences to UTF-8 strings where
  possible.

  The category ids are unchanged from the original, and give the full
  (English) name of the language, optionally combined the name of the
  encoding script.  Note that \sQuote{scots} indicates Scots, the
  Germanic language variety historically spoken in Lowland Scotland and
  parts of Ulster, to be distinguished from Scottish Gaelic (named
  \sQuote{scots_gaelic} in the profiles), the Celtic language variety
  spoken in most of the western Highlands and in the Hebrides (see
  \url{https://en.wikipedia.org/wiki/Scots_language}).
}
\examples{
## Languages in the TC byte profiles:
names(TC_byte_profiles)
## Languages only in the TC byte profiles:
setdiff(names(TC_byte_profiles), names(TC_char_profiles))
## Key options used for the profiles:
attr(TC_byte_profiles, "options")[c("n", "size", "reduce", "useBytes")]
attr(TC_char_profiles, "options")[c("n", "size", "reduce", "useBytes")]
}
\keyword{datasets}
