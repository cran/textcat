\name{TC_profiles}
\alias{TC_char_profiles}
\alias{TC_byte_profiles}
\title{TextCat N-Gram Profiles}
\description{
  TextCat N-gram byte and character profile dbs for language
  identification.
}
\usage{
TC_char_profiles
TC_byte_profiles
}
\details{
  TextCat is a Perl implementation of the Cavnar and Trenkle
  \dQuote{N-Gram-Based Text Categorization} technique by Gertjan van
  Noord \url{http://odur.let.rug.nl/~vannoord/TextCat/} which was
  subsequently integrated into SpamAssassin.  It provides byte N-gram
  profiles for 75 \dQuote{languages} (more precisely, language/encoding
  combinations).

  \code{TC_byte_profiles} provides these byte profiles.

  \code{TC_char_profiles} provides a subset of 56 character profiles
  obtained by converting the byte sequences to UTF-8 strings where
  possible.

  The category ids are unchanged from the original, and give the full
  (English) name of the language, optionally combined the name of the
  encoding script.  Note that \sQuote{scots} indicates Scots, the
  Germanic language variety historically spoken in Lowland Scotland and
  parts of Ulster, to be distinguished from Scottish Gaelic (named
  \sQuote{scots\_gaelic} in the profiles), the Celtic language variety
  spoken in most of the western Highlands and in the Hebrides (see
  \url{http://en.wikipedia.org/wiki/Scots_language}).
}
\keyword{datasets}
