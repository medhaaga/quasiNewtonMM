\name{LBQN}
\alias{LBQN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Limited memory Broyden's quasi-Newton method for accelerating slowly-convergent fixed-point iterations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Partially monotone, acceleration schemes for fast convergence of smooth, monotone, slowly-converging contraction mapping in high dimension where storing the Hessian approximations is difficult. It can be used to accelerate the convergence of a wide variety of fixed-point iterations including the expectation-maximization (EM) algorithms and its variants, majorization-minimization (MM) algorithm.
}
\usage{
LBQN(par, fixptfn, objfn, ..., control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
%%     ~~Describe \code{x} here~~
  A vector of parameters denoting the initial guess for the fixed-point.
  }
  \item{fixptfn}{
%%     ~~Describe \code{x} here~~
  A function $F$ that denotes the fixed-point MM algorithm map. It accepts the     current parameter vector $x1$ of length $p$ and returns the update $x2 = F(x1)$   of same length such that the objective is smaller at $x2$ in majorization        -minimization (and vice versa in minimization-majorization).

  }
  \item{objfn}{
%%     ~~Describe \code{x} here~~
  This is the objective function $L$ that MM algorithm tries to minimize. It takes a vector of parameters as input and returns a scalar that needs to be minimized. It attains its local minimum at the fixed point of the MM algorithm map $F$. Note: in case of maximization problem like maximizing the log-likelihood, this function should be the negative log-likelihood to be minimized.

  }
  \item{control}{
  %%     ~~Describe \code{x} here~~
  A list of control parameters specifing any changes to default values of algorithm control parameters. Full names of control list elements must be specified, otherwise, user-specifications are ignored. See *Details*.
  }
  \item{...}{
  %%     ~~Describe \code{x} here~~
  Arguments passed to fixptfn and objfn.
  }
}
\details{
  These control parameters can be specified by the user, otherwise default values are used.
  \describe{
  \item{m}{
  The number of previous iterates to be stored for LBQN. Suggested values are between 3 and 20. Default is 10.
  }
  \item{objfn.inc}{
  A non-negative scalar that dictates the degree of non-montonicity. Default is 1. Set objfn.inc = 0 to obtain monotone convergence. Setting objfn.inc = Inf gives a non-monotone scheme. In-between values result in partially-monotone convergence.
  }
  \item{tol}{
  A small, positive scalar that determines when iterations should be terminated. Iteration is terminated when abs(x[k] - F(x[k])) <= tol. Default is 1.e-07.
  }
  \item{maxiter}{
  An integer denoting the maximum limit on the number of evaluations of fixptfn, F. Default is 1500.
  }
  \item{obj.stop}{
  A logical variable indicating whether the iterations should be terminated when the decrease in value of objective function is below a defined threshold. Default is FALSE.
  }
  \item{obj.tol}{
  When obj.stop==TRUE. A small, positive scalar that determines when iterations should be terminated. Iteration is terminated when abs(L(x2) - L(x1)) <= obj.tol. Default is 1.e-07.

  }
  \item{verbose}{
  A logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
  }
  \item{intermed}{
  A logical variable denoting whether the intermediate results of iterations       should be returned. If set to TRUE, the function will return a matrix where     each row corresponds to parameters at each iteration, along with the correspondi   ng log-likelihood value. This option is inactive when objfn is not specified.   Default is FALSE.
  }
  }
}
\value{
  A list with the following components:
  \item{par}{
  Parameter, x* that are the fixed-point of F such that x* = F(x*), if convergence is successful.
  }
  \item{value.objfn}{
  The value of the objective function $L$ at termination.
  }
  \item{fpevals}{
  Number of times the fixed-point function fixptfn was evaluated.
  }
  \item{objfevals}{
  Number of times the objective function objfn was evaluated.
  }
  \item{convergence}{
  An integer code indicating type of convergence. 0 indicates successful convergence, whereas 1 denotes failure to converge.
  }
  \item{p.inter}{
  A matrix where each row corresponds to parameters at each iteration, along with the corresponding log-likelihood value. This object is returned only when the control parameter intermed is set to TRUE.
  }

}
\references{
%% ~put references to the literature/web site here ~
Agarwal, M and Xu, J (2020+) Quasi-Newton acceleration of MM algorithm using Broyden's method, \emph{arxiv}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BQN}} which performs the Broyden's quasi-Newton based acceleration scheme for fast convergence of MM algorithm.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


fixptfn <- function(x) x + sin(x)  ## Fixed-point iteration
objfn <- function(x) cos(x)        ## Objective function to be minimized

par <- (pi/2)-.1                   ## Starting point
fp <- LBQN(par, fixptfn = fixptfn, objfn = objfn, control = list())
}
