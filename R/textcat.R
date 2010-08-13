### User interface.

## (Functions to be exported and methods to be registered).

textcat_options <-
local({
    options <-
        list(## options for computing fingerprints
             n = 5L, split = "[[:space:][:punct:][:digit:]]+",
             tolower = TRUE, reduce = TRUE, useBytes = FALSE,
             ignore = "_", size = 1000L,
             ## options for computing distances
             method = "CT")
    function(option, value) {
        if (missing(option)) return(options)
        if (missing(value))
            options[[option]]
        else
            options[[option]] <<- value
    }
})

textcat_profile_db <-
function(x, id, ...)
{
    opts <- fp_options(...)
    fp <- lapply(split(x, rep(id, length.out = length(x))),
                 create_fp, opts)
    .make_textcat_profile_db(fp, names(fp), opts)
}

c.textcat_profile_db <-
function(...)
{
    args <- list(...)    
    ## Ensure common fp options.
    for(nm in fp_option_names) {
        if(length(unique(sapply(args, attr, nm)) > 1L))
            stop(gettextf("Need common '%s'.", nm))
    }
        
    ## What about duplicated names?  Could merge ...
    if(any(duplicated(unlist(lapply(args, names)))))
        stop("Need unique ids.")

    fp <- NextMethod("c")
    
    .make_textcat_profile_db(fp, names(fp), attributes(args[[1L]]))
}

## Could add more methods eventually ...

print.textcat_profile_db <-
function(x, ...)
{
    writeLines(sprintf("A textcat profile db of length %d.", length(x)))
    invisible(x)
}

textcat <-
function(x, p = ECIMCI_profiles, method = "CT")
{
    if(is.null(method))
        method <- textcat_options("method")
    
    ## Use the profile db options for creating the fingerprints.
    x <- lapply(as.character(x), create_fp, attributes(p))

    d <- textcat_xdist(x, p, method)
    ## For now assume that this really does distances.
    pos <- apply(d, 1L,
                 function(d) {
                     pos <- which(d == min(d))
                     if(length(pos) > 1L) NA else pos
                 })
    ifelse(is.na(pos), NA_character_, colnames(d)[pos])
}

### Internal stuff.

fp_option_names <-
    c("n", "split", "tolower", "reduce", "useBytes", "ignore", "size")

fp_options <-
function(...)
{
    opts <- textcat_options()[fp_option_names]
    args <- list(...)
    if(length(args)) { 
        ind <- pmatch(names(args), names(opts), nomatch = 0L)
        opts[ind] <- args[ind > 0L]
    }
    opts
}

create_fp <-
function(x, opts = textcat_options())
{
    marker <- if(opts$reduce) "\1" else "\2"
    fp <- tau::textcnt(as.character(x),
                       n = opts$n, split = opts$split,
                       marker = marker, method = "ngram",
                       useBytes = opts$useBytes, decreasing = TRUE)
    ignore <- opts$ignore
    if(length(ignore))
        fp <- fp[is.na(match(names(fp), ignore))]
    ## Note that the fp length can be smaller than the size.
    size <- opts$size
    if(length(fp) > size)
        fp <- fp[seq_len(size)]
    fp
}

.make_textcat_profile_db <-
function(fps, ids, opts)
{
    val <- fps
    attributes(val) <- opts[fp_option_names]
    names(val) <- ids
    class(val) <- "textcat_profile_db"
    val
}


textcat_xdist <-
function(x, p, method = NULL)
{
    ## Compute distances between collections of profiles.
    if(is.null(method))
        method <- textcat_options("method")
    if(is.character(method))
        method <- textcat_xdist_methods_db[[method[1L]]]
    else if(!is.function(method))
        stop("Invalid 'method'.")
    
    d <- matrix(0, length(x), length(p))
    dimnames(d) <- list(names(x), names(p))    
    if(is.character(x))                 # Be nice ...
        x <- lapply(x, create_fp, attributes(p))
    for(i in seq_along(x))
        for(j in seq_along(p))
            d[i, j] <- method(x[[i]], p[[j]])
    d
}

textcat_xdist_methods_db <- new.env()

## Cavnar-Trenkle out-of-place measure.
textcat_xdist_methods_db$CT <-
function(x, p)
{
    pos <- match(names(x), names(p))
    (sum(abs(pos - seq_along(x)), na.rm = TRUE)
     + length(p) * sum(is.na(pos)))
}

## Some distance measures as mentioned in Singh (2006), "Study Of Some
## Distance Measures For Language And Encoding Identification",
## http://clair.si.umich.edu/clair/anthology/query.cgi?type=Paper&id=W06-1109.
## We expand profiles to a common set of n-grams and, where necessary,
## replace 0 frequencies by (e.g.) 1e-6.

.expand_x_and_p <-
function(x, p, z = 0)
{
    ngrams <- unique(c(names(x), names(p)))
    ind <- match(ngrams, names(x))
    x <- ifelse(is.na(ind), z, x[ind])
    ind <- match(ngrams, names(p))
    p <- ifelse(is.na(ind), z, p[ind])
    list(x = x, p = p)
}

## Cavnar-Trenkle variant: most likely this should be used for CT.
textcat_xdist_methods_db$ranks <-
function(x, p)
{
    e <- .expand_x_and_p(x, p)
    sum(abs(rank(e$x) - rank(e$p)))
}

## Absolute Log Probability Difference.
textcat_xdist_methods_db$ALPD <-
function(x, p)
{
    e <- .expand_x_and_p(x, p, 1e-6)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    sum(abs(log(p) - log(q)))
}

## Kullback-Leibler divergences are a mess, see e.g.
## http://en.wikipedia.org/wiki/Kullback–Leibler_divergence:
## What is commonly known as "K-L divergence" is called "mean
## information for discrimination" in the original reference; the
## symmetric version is called "divergence".
## Let us use the terms I-divergence and J-divergence ...

textcat_xdist_methods_db$KLI <-
function(x, p)
{
    e <- .expand_x_and_p(x, p, 1e-6)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    sum(p * log(p / q))
}    

textcat_xdist_methods_db$KLJ <-
function(x, p)
{
    e <- .expand_x_and_p(x, p, 1e-6)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    sum((p - q) * log(p / q))
}

## Jensen-Shannon divergence, see e.g.
## http://en.wikipedia.org/wiki/Kullback–Leibler_divergence.

textcat_xdist_methods_db$JS <-
function(x, p)
{
    e <- .expand_x_and_p(x, p, 1e-6)
    p <- e$x / sum(e$x)
    q <- e$p / sum(e$p)
    f <- function(t) t * log(t)
    sum(f(p) + f(q)) / 2 - sum(f((p + q) / 2))
}
