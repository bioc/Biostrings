### ==========================================================================
### PairwiseAlignedFixedSubject objects
### --------------------------------------------------------------------------
### An PairwiseAlignedFixedSubject object contains the result of the alignment of
### two XString objects of the same subtype.


setClass("PairwiseAlignedFixedSubject",
    representation(
        pattern="AlignedXStringSet",
        subject="AlignedXStringSet",
        type="character",
        score="numeric",
        substitutionArray="array",
        gapOpening="numeric",
        gapExtension="numeric"
    )
)

setClass("PairwiseAlignedFixedSubjectSummary",
    representation(
        type="character",
        score="numeric",
        nmatch="numeric",
        nmismatch="numeric",
        ninsertion="matrix",
        ndeletion="matrix",
        coverage="XRleInteger",
        mismatchSummary="list"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Initialization.
###

setMethod("initialize", "PairwiseAlignedFixedSubject",
    function(.Object, pattern, subject, type, score, substitutionArray,
             gapOpening, gapExtension, check = TRUE)
    {
        if (!identical(class(unaligned(pattern)), class(unaligned(subject))))
            stop("'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same subtype")
        if (length(type) != 1 || !(type %in% c("global", "local", "overlap", "patternOverlap", "subjectOverlap")))
            stop("'type' must be one of 'global', 'local', 'overlap', 'patternOverlap', or 'subjectOverlap'")
        if (length(pattern) != length(subject))
            stop("'length(pattern)' must equal 'length(subject)'")
        gapOpening <- as.double(- abs(gapOpening))
        if (length(gapOpening) != 1 || is.na(gapOpening))
            stop("'gapOpening' must be a non-positive numeric vector of length 1")
        gapExtension <- as.double(- abs(gapExtension))
        if (length(gapExtension) != 1 || is.na(gapExtension))
            stop("'gapExtension' must be a non-positive numeric vector of length 1")
        slot(.Object, "pattern", check = check) <- pattern
        slot(.Object, "subject", check = check) <- subject
        slot(.Object, "type", check = check) <- type
        slot(.Object, "score", check = check) <- score
        slot(.Object, "substitutionArray", check = check) <- substitutionArray
        slot(.Object, "gapOpening", check = check) <- gapOpening
        slot(.Object, "gapExtension", check = check) <- gapExtension
        .Object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.PairwiseAlignedFixedSubject <- function(object)
{
    message <- character(0)
    if (!identical(class(unaligned(pattern(object))), class(unaligned(subject(object)))))
        message <- c(message, "'unaligned(pattern)' and 'unaligned(subject)' must be XString objects of the same subtype")
    if (length(object@type) != 1 || !(object@type %in% c("global", "local", "overlap", "patternOverlap", "subjectOverlap")))
        message <- c(message, "'type' must be one of 'global', 'local', 'overlap', 'patternOverlap', or 'subjectOverlap'")
    if (length(pattern) != length(subject))
        message <- c(message, "'length(pattern)' must equal 'length(subject)'")
    if (length(object@gapOpening) != 1 || is.na(object@gapOpening))
        message <- c(message, "'gapOpening' must be a non-positive numeric vector of length 1")
    if (length(object@gapExtension) != 1 || is.na(object@gapExtension))
        message <- c(message, "'gapExtension' must be a non-positive numeric vector of length 1")
    if (length(message) == 0)
        message <- NULL
    message
}

setValidity("PairwiseAlignedFixedSubject",
    function(object)
    {
        problems <- .valid.PairwiseAlignedFixedSubject(object)
        if (is.null(problems)) TRUE else problems
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("pattern", "PairwiseAlignedFixedSubject", function(x) x@pattern)
setMethod("subject", "PairwiseAlignedFixedSubject", function(x) x@subject)

setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "PairwiseAlignedFixedSubject", function(x) x@type)
setMethod("type", "PairwiseAlignedFixedSubjectSummary", function(x) x@type)

setGeneric("score", function(x) standardGeneric("score"))
setMethod("score", "PairwiseAlignedFixedSubject", function(x) x@score)
setMethod("score", "PairwiseAlignedFixedSubjectSummary", function(x) x@score)

setMethod("nindel", "PairwiseAlignedFixedSubject",
          function(x) new("InDel", insertion = nindel(subject(x)), deletion = nindel(pattern(x))))
setMethod("nindel", "PairwiseAlignedFixedSubjectSummary",
          function(x) new("InDel", insertion = x@ninsertion, deletion = x@ndeletion))

setMethod("length", "PairwiseAlignedFixedSubject", function(x) length(score(x)))
setMethod("length", "PairwiseAlignedFixedSubjectSummary", function(x) length(score(x)))

setMethod("nchar", "PairwiseAlignedFixedSubject", function(x, type="chars", allowNA=FALSE) nchar(subject(x)))
setMethod("nchar", "PairwiseAlignedFixedSubjectSummary", 
          function(x, type="chars", allowNA=FALSE)
          unname(nmatch(x) + nmismatch(x) + x@ninsertion[,"WidthSum"] + x@ndeletion[,"WidthSum"]))

setMethod("alphabet", "PairwiseAlignedFixedSubject", function(x) alphabet(subject(x)))
setMethod("codec", "PairwiseAlignedFixedSubject", function(x) codec(subject(x)))

setGeneric("pid", signature="x", function(x, type="PID1") standardGeneric("pid"))
setMethod("pid", "PairwiseAlignedFixedSubject",
          function(x, type="PID1") {
              type <- match.arg(type, c("PID1", "PID2", "PID3", "PID4"))
              denom <-
                switch(type,
                       "PID1" = nchar(x),
                       "PID2" = nmatch(x) + nmismatch(x),
                       "PID3" = pmin(nchar(unaligned(pattern(x))), nchar(unaligned(subject(x)))),
                       "PID4" = (nchar(unaligned(pattern(x))) + nchar(unaligned(subject(x)))) / 2)
              100 * nmatch(x)/denom
		  })

setMethod("Views", signature = c(subject = "PairwiseAlignedFixedSubject"),
          function(subject, start=NA, end=NA, names=NULL)
          {
              if (all(is.na(start)))
                  start <- start(subject(subject))
              else if (!is.numeric(start) || length(start) > 1)
                  stop("'start' must be either NA or an integer vector of length 1")
              else
                  start <- as.integer(start) + start(subject(subject))
              if (all(is.na(end)))
                  end <- end(subject(subject))
              else if (!is.numeric(end) || length(end) > 1)
                  stop("'end' must be either NA or an integer vector of length 1")
              else
                  end <- as.integer(end) + start(subject(subject))
              Views(super(unaligned(subject(subject))), start=start, end=end, names=names)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### TODO: Make the "show" method to format the alignment in a SGD fashion
### i.e. split in 60-letter blocks and use the "|" character to highlight
### exact matches.
setMethod("show", "PairwiseAlignedFixedSubject", function(object)
          {
              if (length(object) == 0)
                  cat("Empty Pairwise Alignment\n")
              else {
                  cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                             "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                             "local" = "Local"), " Pairwise Alignment (1 of ", length(object), ")\n", sep = "")
                  if (width(pattern(object))[1] == 0 || width(subject(object))[1] == 0) {
                      patternSpaces <- 0
                      subjectSpaces <- 0
                  } else {
                      patternSpaces <-
                        floor(log10(start(subject(object))[1])) - floor(log10(start(pattern(object))[1]))
		              subjectSpaces <- max(0, - patternSpaces)
		              patternSpaces <- max(0, patternSpaces)
                  }
                  cat(paste(c("pattern: ", rep(" ", patternSpaces)), collapse = ""))
                  show(pattern(object)[1])
                  cat(paste(c("subject: ", rep(" ", subjectSpaces)), collapse = ""))
                  show(subject(object)[1])
                  cat("score:", score(object)[1], "\n")
              }
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "summary" method.
###

setMethod("summary", "PairwiseAlignedFixedSubject", function(object, weight=1L, ...)
          {
              if (!is.numeric(weight) || !(length(weight) %in% c(1, length(object))))
                  stop("'weight' must be an integer vector with length 1 or 'length(object)'")
              if (!is.integer(weight))
                weight <- as.integer(weight)
              if (all(weight == 1))
                  new("PairwiseAlignedFixedSubjectSummary",
                      type = type(object),
                      score = score(object),
                      nmatch = nmatch(object),
                      nmismatch = nmismatch(object),
                      ninsertion = nindel(subject(object)),
                      ndeletion = nindel(pattern(object)),
                      coverage = coverage(object),
                      mismatchSummary = mismatchSummary(object))
              else
                  new("PairwiseAlignedFixedSubjectSummary",
                      type = type(object),
                      score = rep(score(object), weight),
                      nmatch = rep(nmatch(object), weight),
                      nmismatch = rep(nmismatch(object), weight),
                      ninsertion = nindel(subject(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      ndeletion = nindel(pattern(object))[rep(seq_len(length(object)), weight), , drop = FALSE],
                      coverage = coverage(object, weight = weight),
                      mismatchSummary = mismatchSummary(object, weight = weight))
          })

setMethod("show", "PairwiseAlignedFixedSubjectSummary", function(object)
          {
              cat(switch(type(object), "global" = "Global", "overlap" = "Overlap",
                         "patternOverlap" = "Pattern Overlap", "subjectOverlap" = "Subject Overlap",
                         "local" = "Local"), " Pairwise Alignment\n", sep = "")
              cat("Number of Alignments:  ", length(score(object)), "\n", sep = "")
              cat("\nScores:\n")
              print(summary(score(object)))
              cat("\nNumber of matches:\n")
              print(summary(nmatch(object)))
              n <- min(nrow(mismatchSummary(object)[["subject"]]), 10)
              cat(paste("\nTop", n, "Mismatch Counts:\n"))
              print(mismatchSummary(object)[["subject"]][
                order(mismatchSummary(object)[["subject"]][["Count"]],
                      mismatchSummary(object)[["subject"]][["Probability"]],
                      decreasing = TRUE)[seq_len(n)],,drop=FALSE])
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.character" method.
###

setMethod("aligned", "PairwiseAlignedFixedSubject",
          function(x) {
              codecX <- codec(x)
              if (is.null(codecX)) {
                  gapCode <- charToRaw("-")
              } else {
                  letters2codes <- codecX@codes
                  names(letters2codes) <- codecX@letters
                  gapCode <- as.raw(letters2codes[["-"]])
              }
              .Call("PairwiseAlignedFixedSubject_align_aligned", x, gapCode, PACKAGE="Biostrings")
          })

setMethod("as.character", "PairwiseAlignedFixedSubject",
          function(x)
          {
              as.character(aligned(x))
          })

setMethod("toString", "PairwiseAlignedFixedSubject", function(x, ...) toString(as.character(x), ...))

setMethod("as.matrix", "PairwiseAlignedFixedSubject",
          function(x) {
              as.matrix(aligned(x))
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "PairwiseAlignedFixedSubject",
    function(x, i, j, ..., drop)
    {
        if (!missing(j) || length(list(...)) > 0)
            stop("invalid subsetting")
        if (missing(i) || (is.logical(i) && all(i)))
            return(x)
        if (is.logical(i))
            i <- which(i)
        if (!is.numeric(i) || any(is.na(i)))
            stop("invalid subsetting")
        if (any(i < 1) || any(i > length(x)))
            stop("subscript out of bounds")
        new("PairwiseAlignedFixedSubject",
            pattern = x@pattern[i],
            subject = x@subject[i],
            type = x@type,
            score = x@score[i],
            substitutionArray = x@substitutionArray,
            gapOpening = x@gapOpening,
            gapExtension = x@gapExtension)
    }
)

setReplaceMethod("[", "PairwiseAlignedFixedSubject",
    function(x, i, j,..., value)
    {
        stop("attempt to modify the value of a ", class(x), " instance")
    }
)

setMethod("rep", "PairwiseAlignedFixedSubject",
    function(x, times)
		x[rep.int(seq_len(length(x)), times)]
)
