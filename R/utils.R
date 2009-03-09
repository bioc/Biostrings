### Some low-level (not exported) helper functions

isNumericOrNAs <- function(x)
{
    is.numeric(x) || (is.atomic(x) && is.vector(x) && all(is.na(x)))
}

normargStart <- function(start)
{
    if (!isSingleNumber(start))
        stop("'start' must be a single integer")
    if (!is.integer(start))
        start <- as.integer(start)
    if (start < 1L)
        stop("'start' must be >= 1")
    start
}

normargNchar <- function(start, nchar, seq_nchar)
{
    if (!isSingleNumberOrNA(nchar))
        stop("'nchar' must be a single integer or NA")
    if (is.na(nchar)) {
        nchar <- seq_nchar - start + 1L
        if (nchar < 0L)
            stop("cannot read a negative number of letters")
        return(nchar)
    }
    if (!is.integer(nchar))
        nchar <- as.integer(nchar)
    if (nchar < 0L)
        stop("cannot read a negative number of letters")
    end <- start + nchar - 1L
    if (end > seq_nchar)
        stop("cannot read beyond the end of 'seq'")
    nchar
}

normargUseNames <- function(use.names)
{
    if (is.null(use.names))
        return(TRUE)
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    use.names
}

### Returns an integer vector.
pow.int <- function(x, y)
{
    if (!is.numeric(x))
        stop("'x' must be a numeric vector")
    if (!is.integer(x))
        x <- as.integer(x)
    ans <- rep.int(1L, length(x))
    for (i in seq_len(y))
        ans <- ans * x
    ans
}

