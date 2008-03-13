### =========================================================================
### The BStringViews class
### -------------------------------------------------------------------------
###
### The BStringViews class is the basic container for storing a set of views
### (start/end locations) on the same XString object, called the "subject"
### string.
###

setClass("BStringViews",
    contains="IRanges",
    representation(
        subject="XString"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setGeneric("subject", function(x) standardGeneric("subject"))
setMethod("subject", "BStringViews", function(x) x@subject)

### The "views" generic and methods are commented for now because they
### conflict with the views() function currently defined in the
### BStringViews-constructors.R file.
#setGeneric("views", function(x) standardGeneric("views"))
#setMethod("views", "BStringViews", function(x) x@ranges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.BStringViews.subject <- function(object)
{
    if (!is(subject(object), "XString"))
        return("the subject must be an XString object")
    NULL
}

.valid.BStringViews.width <- function(object)
{
    if (length(width(object)) != 0 && min(width(object)) < 1L)
        return("null widths are not allowed")
    NULL
}

.valid.BStringViews <- function(object)
{
    #cat("validating BStringViews object...\n")
    c(.valid.BStringViews.subject(object),
      .valid.BStringViews.width(object))
}

setValidity("BStringViews",
    function(object)
    {
        problems <- .valid.BStringViews(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization (not intended to be used directly by the user).
###

setMethod("initialize", "BStringViews",
    function(.Object, subject, start=integer(0), width=integer(0),
                               desc=NULL, check=TRUE)
    {
        .Object <- callNextMethod(.Object, start=start, width=width,
                                           names=desc, check=check)
        slot(.Object, "subject", check=FALSE) <- subject
        if (check) {
            ## I found that using validObject() in "initialize" doesn't work
            ## properly (validation is called too many times and not in an
            ## order that makes sense to me...)
            #validObject(.Object)
            problems <- .valid.BStringViews(.Object)
            if (!is.null(problems)) stop(paste(problems, collapse="\n  "))
        }
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Replacement methods.
###
### BStringViews objects inherit the replacement methods defined for parent
### class IRanges (see IRanges-class.R).
### However, the "width" method needs to be overridden because of the
### additional constraint that applies to BStringViews objects.
###

setReplaceMethod("width", "BStringViews",
    function(x, check=TRUE, value)
    {
        x <- callNextMethod(x, check=TRUE, value)
        if (check) {
            problem <- .valid.BStringViews.width(x)
            if (!is.null(problem)) stop(problem)
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "nchar" method.
###

setMethod("nchar", "BStringViews",
    function(x, type="chars", allowNA=FALSE)
    {
        if (length(x) == 0)
            return(integer(0))
        start0 <- pmax.int(start(x), 1L)
        end0 <- pmin.int(end(x), nchar(subject(x)))
        ans <- end0 - start0 + 1L
        ans[ans < 0L] <- 0L
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### The 2 helper functions below convert a given view on an XString object
### into a character-string.
### Both assume that 'start' <= 'end' (so they don't check it) and
### padd the result with spaces to produce the "margin effect"
### if 'start' or 'end' are out of limits.

### nchar(BStringViews.get_view(x, start, end)) is always end-start+1
BStringViews.get_view <- function(x, start, end)
{
    lx <- length(x)
    if (end < 1 || start > lx)
            return(format("", width=end-start+1))
    Lmargin <- ""
    if (start < 1) {
        Lmargin <- format("", width=1-start)
        start <- 1
    }
    Rmargin <- ""
    if (end > lx) {
        Rmargin <- format("", width=end-lx)
        end <- lx
    }
    paste(Lmargin, XString.read(x, start, end), Rmargin, sep="")
}

### nchar(BStringViews.get_snippet(x, start, end, snippetWidth)) is <= snippetWidth
BStringViews.get_snippet <- function(x, start, end, snippetWidth)
{
    if (snippetWidth < 7)
        snippetWidth <- 7
    width <- end - start + 1
    if (width <= snippetWidth) {
        BStringViews.get_view(x, start, end)
    } else {
        w1 <- (snippetWidth - 2) %/% 2
        w2 <- (snippetWidth - 3) %/% 2
        paste(BStringViews.get_view(x, start, start+w1-1),
              "...",
              BStringViews.get_view(x, end-w2+1, end), sep="")
    }
}

BStringViews.show_vframe_header <- function(iW, startW, endW, widthW)
{
    cat(format("", width=iW+1),
        format("start", width=startW, justify="right"), " ",
        format("end", width=endW, justify="right"), " ",
        format("width", width=widthW, justify="right"), "\n",
        sep="")
}

BStringViews.show_vframe_line <- function(x, i, iW, startW, endW, widthW)
{
    start <- start(x)[i]
    end <- end(x)[i]
    width <- end - start + 1
    snippetWidth <- getOption("width") - 6 - iW - startW - endW - widthW
    cat(format(paste("[", i,"]", sep=""), width=iW, justify="right"), " ",
        format(start, width=startW, justify="right"), " ",
        format(end, width=endW, justify="right"), " ",
        format(width, width=widthW, justify="right"), " ",
        "[", BStringViews.get_snippet(subject(x), start, end, snippetWidth), "]\n",
        sep="")
}

### 'half_nrow' must be >= 1
BStringViews.show_vframe <- function(x, half_nrow=9L)
{
    cat("\nviews:")
    lx <- length(x)
    if (lx == 0)
        cat(" NONE\n")
    else {
        cat("\n")
        iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
        startMax <- max(start(x))
        startW <- max(nchar(startMax), nchar("start"))
        endMax <- max(end(x))
        endW <- max(nchar(endMax), nchar("end"))
        widthMax <- max(width(x))
        widthW <- max(nchar(widthMax), nchar("width"))
        BStringViews.show_vframe_header(iW, startW, endW, widthW)
        if (lx <= 2*half_nrow+1) {
            for (i in seq_len(lx))
                BStringViews.show_vframe_line(x, i, iW, startW, endW, widthW)
        } else {
            for (i in 1:half_nrow)
                BStringViews.show_vframe_line(x, i, iW, startW, endW, widthW)
            cat(format("...", width=iW, justify="right"),
                " ",
                format("...", width=startW, justify="right"),
                " ",
                format("...", width=endW, justify="right"),
                " ",
                format("...", width=widthW, justify="right"),
                " ...\n", sep="")
            for (i in (lx-half_nrow+1L):lx)
                BStringViews.show_vframe_line(x, i, iW, startW, endW, widthW)
        }
    }
}

setMethod("show", "BStringViews",
    function(object)
    {
        subject <- subject(object)
        lsub <- length(subject)
        cat("  Views on a ", lsub, "-letter ", class(subject), " subject", sep="")
        cat("\nsubject:", XString.get_snippet(subject, getOption("width") - 9))
        BStringViews.show_vframe(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Extract the i-th views of a BStringViews object as an XString object.
### Return an XString object of the same subtype as subject(x).
### Example:
###   bs <- BString("ABCD-1234-abcd")
###   bsv <- views(bs, 1:7, 13:7)
###   bsv[[3]]
###   bsv[[0]] # Return bs, same as subject(bsv)
###   views(bs)[[1]] # Returns bs too!
###
### Supported 'i' types: numeric vector of length 1.
setMethod("[[", "BStringViews",
    function(x, i, j, ...)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i) || !is.numeric(i))
            stop("invalid subscript type")
        if (length(i) < 1L)
            stop("attempt to select less than one element")
        if (length(i) > 1L)
            stop("attempt to select more than one element")
        if (i == 0)
            return(subject(x))
        if (i < 1L || i > length(x))
            stop("subscript out of bounds")
        start <- start(x)[i]
        end <- end(x)[i]
        if (start < 1L || end > length(subject(x)))
            stop("view is out of limits")
        XString.substr(subject(x), start, end)
    }
)

setReplaceMethod("[[", "BStringViews",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", sQuote(class(x)), " object")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality.
###

### Assume that 'start1', 'end1', 'start2', 'end2' are single integers
### and that start1 <= end1 and start2 <= end2.
BStringViews.view1_equal_view2 <- function(x1, start1, end1, x2, start2, end2)
{
    one <- as.integer(1)
    w1 <- end1 - start1 + one
    w2 <- end2 - start2 + one
    if (w1 != w2)
        return(FALSE)

    lx1 <- length(x1)
    isBlank1 <- end1 < one || start1 > lx1
    lx2 <- length(x2)
    isBlank2 <- end2 < one || start2 > lx2
    if (isBlank1 && isBlank2)
        return(TRUE)
    if (isBlank1 || isBlank2)
        return(FALSE)

    # Left margin
    LmarginSize1 <- start1 < one
    LmarginSize2 <- start2 < one
    if (LmarginSize1 != LmarginSize2)
        return(FALSE)
    if (LmarginSize1) {
        # Both views have a left margin
        if (start1 != start2)
            return(FALSE)
        start1 <- one
        start2 <- one
    }

    # Right margin
    RmarginSize1 <- end1 > lx1
    RmarginSize2 <- end2 > lx2
    if (RmarginSize1 != RmarginSize2)
        return(FALSE)
    if (RmarginSize1) {
        # Both views have a right margin
        if (end1 - lx1 != end2 - lx2)
            return(FALSE)
        end1 <- lx1
        end2 <- lx2
    }

    # At this point, we can trust that 1 <= start1 <= end1 <= lx1
    # and that 1 <= start2 <= end2 <= lx2 so we can call unsafe
    # function XString.substr() with no fear...
    XString.substr(x1, start1, end1) == XString.substr(x2, start2, end2)
}

### 'x' and 'y' must be BStringViews objects.
### Returns a logical vector of length max(length(x), length(y)).
### Recycle its arguments.
BStringViews.equal <- function(x, y)
{
    lx <- length(x)
    ly <- length(y)
    if (lx < ly) {
        tmp <- x
        x <- y
        y <- tmp
        tmp <- lx
        lx <- ly
        ly <- tmp
    }
    if (ly == 0)
        return(logical(0))
    # Now we are sure that lx >= ly >= 1
    ans <- logical(lx)
    j <- 1
    for (i in seq_len(lx)) {
        ans[i] <- BStringViews.view1_equal_view2(
                      subject(x), start(x)[i], end(x)[i],
                      subject(y), start(y)[j], end(y)[j])
        # Recycle
        if (j < ly) j <- j + 1 else j <- 1
    }
    if (j != 1)
        warning(paste("longer object length",
                      "is not a multiple of shorter object length"))
    ans
}

### These methods are called if at least one side of the "==" (or "!=")
### operator is a BStringViews object. They have precedence over the
### corresponding methods defined for XString objects, i.e. they will
### be called if one side is a BStringViews object and the other side
### is an XString object.
### Typical use:
###   v <- views(DNAString("TAATAATG"), -2:9, 0:11)
###   v == v[4]
###   v == v[1]
###   v2 <- views(DNAString("G"), 1, 3)
###   v == v2
### Also works if one side is an XString object:
###   v == DNAString("ATG")
###   RNAString("AUG") == v
### Whitespace matters:
###   v == "TG"
### But this doesn't work neither ("TG " can't be converted to a DNAString
### object):
###   v == "TG "

setMethod("==", signature(e1="BStringViews", e2="BStringViews"),
    function(e1, e2)
    {
        if (!comparableXStrings(subject(e1), subject(e2))) {
            class1 <- class(subject(e1))
            class2 <- class(subject(e2))
            stop("comparison between BStringViews objects with subjects of ",
                 "class \"", class1, "\" and \"", class2, "\" ",
                 "is not supported")
        }
        BStringViews.equal(e1, e2)
    }
)
setMethod("==", signature(e1="BStringViews", e2="XString"),
    function(e1, e2)
    {
        if (!comparableXStrings(subject(e1), e2)) {
            class1 <- class(subject(e1))
            class2 <- class(e2)
            stop("comparison between a BStringViews object with a subject of ",
                 "class \"", class1, "\" and a \"", class2, "\" instance ",
                 "is not supported")
        }
        BStringViews.equal(e1, BStringViews(e2))
    }
)
setMethod("==", signature(e1="BStringViews", e2="character"),
    function(e1, e2)
    {
        if (!is(subject(e1), "BString"))
            stop("comparison between a BStringViews object with a subject of ",
                 "class \"", class(subject(e1)), "\" and a character vector ",
                 "is not supported")
        if (length(e2) == 0 || any(e2 %in% c("", NA)))
            stop("comparison between a BStringViews object and a character ",
                 "vector of length 0 or with empty strings or NAs ",
                 "is not supported")
        BStringViews.equal(e1, BStringViews(e2, subjectClass="BString"))
    }
)
setMethod("==", signature(e1="XString", e2="BStringViews"),
    function(e1, e2) e2 == e1
)
setMethod("==", signature(e1="character", e2="BStringViews"),
    function(e1, e2) e2 == e1
)

setMethod("!=", signature(e1="BStringViews", e2="BStringViews"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="BStringViews", e2="XString"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="BStringViews", e2="character"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="XString", e2="BStringViews"),
    function(e1, e2) !(e1 == e2)
)
setMethod("!=", signature(e1="character", e2="BStringViews"),
    function(e1, e2) !(e1 == e2)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other standard generic methods.
###

setMethod("as.character", "BStringViews",
    function(x, use.names=TRUE, check.limits=TRUE)
    {
        use.names <- normalize.use.names(use.names)
        if (check.limits)
            ans <- sapply(seq_len(length(x)), function(i) as.character(x[[i]]))
        else
            ans <- sapply(seq_len(length(x)),
                          function(i) BStringViews.get_view(subject(x), start(x)[i], end(x)[i]))
        if (use.names)
            names(ans) <- desc(x)
        ans
    }
)

### Supported modes: "integer" (default) and "character".
setMethod("as.matrix", "BStringViews",
    function(x, mode="integer", use.names=TRUE, check.limits=TRUE)
    {
        if (!is.character(mode) || length(mode) != 1
         || !(mode %in% c("integer", "character")))
            stop("'mode' must be either \"integer\" or \"character\"")
        use.names <- normalize.use.names(use.names)
        if (mode == "integer")
            return(callNextMethod())
        nrow <- length(x)
        if (nrow == 0)
            stop("'x' must contain at least 1 view")
        widths <- width(x)
        ncol <- widths[1]
        if (!all(widths == ncol))
            stop("'x' views are not equal-width")
        y <- as.character(x, use.names=FALSE, check.limits=check.limits)
        y <- unlist(strsplit(y, NULL), recursive=FALSE, use.names=FALSE)
        m <- matrix(y, nrow=nrow, byrow=TRUE)
        if (use.names)
            rownames(m) <- desc(x)
        m
    }
)

### This overrides the "as.list" method for IRanges objects with a
### totally different semantic!
setMethod("as.list", "BStringViews",
    function(x)
    {
        lx <- length(x)
        ans <- vector("list", lx)
        for (i in 1:lx) {
            ans[[i]] <- x[[i]]
        }
        ans
    }
)

setMethod("toString", "BStringViews",
    function(x, ...)
    {
        toString(as.character(x), ...)
    }
)
