\name{prime_allocate}
\alias{prime_allocate}
\keyword{internal}
\title{Allocate Precomputed Subset Values}
\usage{
prime_allocate(subset_values, zero_tolerance = sqrt(.Machine$double.eps))
}
\arguments{
  \item{subset_values}{Finite numeric vector ordered by binary subset masks from
  zero through \eqn{2^d-1}. Its length must be \eqn{2^d}.}
  \item{zero_tolerance}{Nonnegative relative threshold used to identify null
  subsets after projection.}
}
\value{
A named numeric allocation. The importance values sum to the projected
full-set value. Attributes contain:
\itemize{
  \item \code{projected_values}: the nonnegative monotone projection;
  \item \code{projection}: iteration count, convergence flag, and maximum
  remaining violation.
}
}
\description{
Projects supplied subset values onto the set of
nonnegative monotone set functions and evaluates the canonical
proportional-value extension by ratio-potential recursion.
}
\details{
The first value is set to zero. The full-set value is preserved by the
projection. Computation and storage scale exponentially with dimension
\eqn{d}.
}
\references{
Herin, M., Il Idrissi, M., Chabridon, V., and Iooss, B. (2024).
Proportional marginal effects for global sensitivity analysis.
\emph{SIAM/ASA Journal on Uncertainty Quantification}, 12(2), 667--692.
\doi{10.1137/22M153032X}.
}
\examples{
values <- c(0, 0.25, 0.60, 1)
allocation <- prime_allocate(values)
allocation
attr(allocation, "projected_values")
}
