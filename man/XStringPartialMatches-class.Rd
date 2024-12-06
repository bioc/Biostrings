\name{XStringPartialMatches-class}
\docType{class}

% Classes:
\alias{class:XStringPartialMatches}
\alias{XStringPartialMatches-class}

% Accessor methods:
\alias{subpatterns}
\alias{subpatterns,XStringPartialMatches-method}

% Standard generic methods:
\alias{show,XStringPartialMatches-method}
\alias{[,XStringPartialMatches-method}

\title{XStringPartialMatches objects}

\description{
  WARNING: XStringPartialMatches objects are deprecated!
}

\section{Accessor methods}{
  In the code snippets below,
  \code{x} is an XStringPartialMatches object.

  \describe{
    \item{\code{subpatterns(x)}:}{
      Not ready yet.
    }
    \item{\code{pattern(x)}:}{
      Not ready yet.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x} is an XStringPartialMatches objects,
  and \code{i} can be a numeric or logical vector.

  \describe{
    \item{\code{x[i]}:}{
      Return a new XStringPartialMatches object made of the selected views.
      \code{i} can be a numeric vector, a logical vector, \code{NULL} or missing.
      The returned object has the same subject as \code{x}.
    }
  }
}

\author{H. Pagès}

\seealso{
  \link{XStringViews-class},
  \link{XString-class},
  \code{\link{letter}}
}

\keyword{methods}
\keyword{classes}
