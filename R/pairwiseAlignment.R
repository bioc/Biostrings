### =========================================================================
### The pairwiseAlignment() generic & related functions
### -------------------------------------------------------------------------
###
### The pairwiseAligment() function provides optimal pairwise alignment of
### the following types:
### - Global alignment
### - Local alignment
### - Overlap alignment
###
### -------------------------------------------------------------------------


XString.pairwiseAlignment <-
function(pattern,
         subject,
         patternQuality = 22L,
         subjectQuality = 22L,
         type = "global",
         substitutionMatrix = NULL,
         gapOpening = -10,
         gapExtension = -4,
         scoreOnly = FALSE)
{
  ## Check arguments
  if (class(pattern) != class(subject))
    stop("'pattern' and 'subject' must have the same class")
  type <- match.arg(tolower(type), c("global", "local", "overlap"))
  typeCode <- c("global" = 1L, "local" = 2L, "overlap" = 3L)[[type]]
  gapOpening <- as.double(- abs(gapOpening))
  if (length(gapOpening) != 1 || is.na(gapOpening))
    stop("'gapOpening' must be a non-positive numeric vector of length 1")
  gapExtension <- as.double(- abs(gapExtension))
  if (length(gapExtension) != 1 || is.na(gapExtension))
    stop("'gapExtension' must be a non-positive numeric vector of length 1")
  scoreOnly <- as.logical(scoreOnly)
  if (length(scoreOnly) != 1 || any(is.na(scoreOnly)))
    stop("'scoreOnly' must be a non-missing logical value")

  ## Process string information
  if (is.null(codec(pattern))) {
    uniqueBases <-
      unique(c(unique(charToRaw(as.character(pattern))),
               unique(charToRaw(as.character(subject)))))
    alphabetToCodes <- as.integer(uniqueBases)
    names(alphabetToCodes) <- rawToChar(uniqueBases, multiple = TRUE)
    gapCode <- charToRaw("-")
  } else {
    stringCodec <- codec(pattern)
    alphabetToCodes <- stringCodec@codes
    names(alphabetToCodes) <- stringCodec@letters
    gapCode <- as.raw(alphabetToCodes["-"])
  }

  ## Generate quality-based and constant substitution matrix information
  if (is.null(substitutionMatrix)) {
    if (is.numeric(patternQuality)) {
      if (any(is.na(patternQuality)) || any(patternQuality < 0 || patternQuality > 99))
        stop("integer 'patternQuality' values must be between 0 and 99")
      patternQuality <- rawToChar(as.raw(33L + as.integer(patternQuality)))
    }
    if (!is(patternQuality, "XString"))
      patternQuality <- BString(patternQuality)

    if (is.numeric(subjectQuality)) {
      if (any(is.na(subjectQuality)) || any(subjectQuality < 0 || subjectQuality > 99))
        stop("integer 'subjectQuality' values must be between 0 and 99")
      subjectQuality <- rawToChar(as.raw(33L + as.integer(subjectQuality)))
    }
    if (!is(subjectQuality, "XString"))
      subjectQuality <- BString(subjectQuality)

    nAlphabet <-
      switch(class(pattern),
             DNAString =, RNAString = 4L,
             AAString = 20L,
             length(alphabetToCodes))

    errorProbs <- 10^seq(0, -9.9, by = -0.1)
    errorMatrix <-
      outer(errorProbs, errorProbs,
            function(e1,e2,n) e1 + e2 - (n/(n - 1)) * e1 * e2,
            n = nAlphabet)
    qualityLookupTable <- buildLookupTable(33:(99 + 33), 0:99)
    qualityMatchMatrix <- log2((1 - errorMatrix) * nAlphabet)
    qualityMismatchMatrix <- log2(errorMatrix * (nAlphabet / (nAlphabet - 1)))

    constantLookupTable <- integer(0)
    constantMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  } else {
    patternQuality <- BString("")
    subjectQuality <- BString("")
    qualityLookupTable <- integer(0)
    qualityMatchMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    qualityMismatchMatrix <- matrix(numeric(0), nrow = 0, ncol = 0)
    if (is.character(substitutionMatrix)) {
      if (length(substitutionMatrix) != 1)
        stop("'substitutionMatrix' is a character vector of length != 1")
      tempMatrix <- substitutionMatrix
      substitutionMatrix <- try(getdata(tempMatrix), silent = TRUE)
      if (is(substitutionMatrix, "try-error"))
        stop("unknown scoring matrix \"", tempMatrix, "\"")
    }
    if (!is.matrix(substitutionMatrix) || !is.numeric(substitutionMatrix))
      stop("'substitutionMatrix' must be a numeric matrix")
    if (!identical(rownames(substitutionMatrix), colnames(substitutionMatrix)))
      stop("row and column names differ for matrix 'substitutionMatrix'")
    if (is.null(rownames(substitutionMatrix)))
      stop("matrix 'substitutionMatrix' must have row and column names")
    if (any(duplicated(rownames(substitutionMatrix))))
      stop("matrix 'substitutionMatrix' has duplicated row names")
    availableLetters <-
      intersect(names(alphabetToCodes), rownames(substitutionMatrix))
    constantMatrix <-
      matrix(as.double(substitutionMatrix[availableLetters, availableLetters]),
             nrow = length(availableLetters),
             ncol = length(availableLetters),
             dimnames = list(availableLetters, availableLetters))
    constantLookupTable <-
      buildLookupTable(alphabetToCodes[availableLetters],
                       0:(length(availableLetters) - 1))
  }
  answer <- .Call("align_pairwiseAlignment",
                  pattern,
                  subject,
                  patternQuality,
                  subjectQuality,
                  gapCode,
                  typeCode,
                  scoreOnly,
                  gapOpening,
                  gapExtension,
                  qualityLookupTable,
                  qualityMatchMatrix,
                  qualityMismatchMatrix,
                  dim(qualityMatchMatrix),
                  constantLookupTable,
                  constantMatrix,
                  dim(constantMatrix),
                  PACKAGE="Biostrings")
  if (scoreOnly) {
    output <- answer[["score"]]
  } else {
    align1 <-
      new(class(pattern),
          xdata = answer[["align1"]],
          length = length(answer[["align1"]]))
    align2 <-
      new(class(subject),
          xdata = answer[["align2"]],
          length = length(answer[["align2"]]))
    output <- new("XStringAlign",
                  align1 = align1,
                  align2 = align2,
                  patternQuality = patternQuality,
                  subjectQuality = subjectQuality,
                  type = type,
                  score = answer[["score"]],
                  constantMatrix = constantMatrix,
                  gapOpening = gapOpening,
                  gapExtension = gapExtension)
  }
  return(output)
}


setGeneric("pairwiseAlignment", signature = c("pattern", "subject"),
           function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                    type = "global", substitutionMatrix = NULL,
                    gapOpening = -10, gapExtension = -4,
                    scoreOnly = FALSE)
           standardGeneric("pairwiseAlignment"))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(BString(pattern), BString(subject),
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "character", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(XString(class(subject), pattern), subject,
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "character"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(pattern, XString(class(pattern), subject),
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))

setMethod("pairwiseAlignment",
          signature(pattern = "XString", subject = "XString"),
          function(pattern, subject, patternQuality = 22L, subjectQuality = 22L,
                   type = "global", substitutionMatrix = NULL,
                   gapOpening = -10, gapExtension = -4,
                   scoreOnly = FALSE)
          XString.pairwiseAlignment(pattern, subject,
                                    patternQuality = patternQuality,
                                    subjectQuality = subjectQuality,
                                    type = type,
                                    substitutionMatrix = substitutionMatrix,
                                    gapExtension = gapExtension,
                                    gapOpening = gapOpening,
                                    scoreOnly = scoreOnly))
