\name{prime_screen}
\alias{prime_screen}
\title{Nonparametric Variable Selection}
\usage{
prime_screen(X, y, rescale = TRUE, twin_r = NULL, lambda = 1e-6)
}
\arguments{
  \item{X}{Numeric predictor matrix or data frame with unique column names.}
  \item{y}{Finite numeric response vector with length \code{nrow(X)}.}
  \item{rescale}{Logical; standardize predictors before distance calculations.}
  \item{twin_r}{Optional integer inverse Data-Twinning split ratio. \code{NULL}
  uses every observation as a query location.}
  \item{lambda}{Nonnegative slope-only numerical nugget used in local-linear
  fits.}
}
\value{
An object of class \code{prime_screen} with components:
\itemize{
  \item \code{call}: matched function call;
  \item \code{candidates} and \code{candidate_names}: inclusive forward set;
  \item \code{selected} and \code{selected_names}: backward-cleaned set;
  \item \code{forward}: ridge seed, local path, stopping reason, and trace;
  \item \code{backward}: retained set and deletion-comparison trace;
  \item \code{settings}: neighborhood, query, nugget, and calibration rules.
}
}
\description{
Implements the nonparametric model-free variable selection technique in
Cho and Joseph (2026).
}
\details{
The ridge ordering is searched through
\eqn{\lceil1.5\log\{\max(2,p)\}\rceil} coordinates, including the intercept-only model,
and BIC selects the initial prefix. The local path adds the coordinate with the
largest explained-variance gain when its paired querywise calibrated gain is
positive. Backward reduction sequentially removes zero calibrated deletion
losses and recomputes all projected neighborhoods after every removal.

Use \code{summary(screen)} to see the ridge seed, forward path, retained
coordinates, and backward-reduction counts.
}
\references{
Cho, Y. and Joseph, V. R. (2026). Measuring Variable Importance Using
Proportional Inclusive Marginal Effects. arXiv:2610.17053.
\doi{10.48550/arXiv.2610.17053}.
}
\seealso{\code{\link{prime}}}
\examples{
\donttest{
set.seed(5)
X <- matrix(runif(960, -1, 1), nrow = 120)
colnames(X) <- paste0("X", seq_len(ncol(X)))
y <- X[, 1] - X[, 2] + X[, 5]^2 + rnorm(120, sd = 0.08)

screen <- prime_screen(X, y)
screen$selected_names
summary(screen)
screen$forward$trace
screen$backward$trace
}
}
