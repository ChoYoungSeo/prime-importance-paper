\name{prime-package}
\alias{prime-package}
\docType{package}
\title{PRIME: Proportional Inclusive Marginal Effects}
\description{
The \pkg{prime} package implements the proportional inclusive marginal effects
of Cho and Joseph (2026), which estimates model-free variable importance
from noisy input-output data. It also provides model-free variable selection.
}
\details{
PRIME estimates explained-variance subset values with local-linear
neighborhoods, constructs a BIC-gated ridge-start candidate set, applies
paired querywise forward and backward calibration, and evaluates proportional
marginal effects from noise-adjusted risks on the retained coordinates.

The main interface is \code{prime()}, which runs coordinate selection,
subset-value estimation, and proportional allocation. The supporting interface
is \code{prime_screen()} for coordinate selection and diagnostic traces
without the allocation step. Use \code{summary(fit)} or
\code{summary(screen)} for a concise result summary, and \code{plot(fit)}
to display variable importance.

Exact subset enumeration grows exponentially with the retained dimension. The
ambient dimension affects screening, while the selected dimension determines
the final allocation cost. Data Twinning (Vakayil and Joseph, 2022) can reduce
the query-average cost for large datasets. Fitted objects have concise print,
summary, and plot methods.
}
\references{
Herin, M., Il Idrissi, M., Chabridon, V., and Iooss, B. (2024).
Proportional marginal effects for global sensitivity analysis.
\emph{SIAM/ASA Journal on Uncertainty Quantification}, 12(2), 667--692.
\doi{10.1137/22M153032X}.

Cho, Y. and Joseph, V. R. (2026). Measuring Variable Importance Using
Proportional Inclusive Marginal Effects. arXiv:2610.17053.
\doi{10.48550/arXiv.2610.17053}.

Vakayil, A. and Joseph, V. R. (2022). Data Twinning.
\emph{Statistical Analysis and Data Mining: The ASA Data Science Journal},
15, 598--610. \doi{10.1002/sam.11574}.
}
\section{Acknowledgments}{
This research is supported by a U.S. National Science Foundation grant
DMS-2310637.
}
\seealso{
\code{\link{prime}}, \code{\link{plot.prime_fit}},
\code{\link{prime_screen}}
}
\author{
Youngseo Cho and V. Roshan Joseph
}
