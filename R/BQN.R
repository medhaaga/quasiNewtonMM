#############################################################
## Broyden's quasi-Newton method
#############################################################

BQN <- function(par, fixptfn, objfn, ... , control=list()) {

  control.default <- list(objfn.inc=1, qn=1, tol=1e-07, obj.tol=1e-07, obj.stop=FALSE, step.max=1000, step.min=1, maxiter=1500, verbose=FALSE, intermed=FALSE)

  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  ctrl <- modifyList(control.default, control)

  maxiter <- ctrl$maxiter
  tol <- ctrl$tol
  objfn.inc <- ctrl$objfn.inc
  verbose <- ctrl$verbose
  intermed <- ctrl$intermed
  obj.stop <- ctrl$obj.stop
  obj.tol <- ctrl$obj.tol
  step.max <- ctrl$step.max
  step.min <- ctrl$step.min
  qn <- ctrl$qn

  if (qn > length(par))
    stop("q cannot be bigger than problem dimension")
  if (verbose) cat("BQN \n")
  if (missing(objfn)) stop("\n objective function is not available \n\n")
  if(missing(par)) stop("Missing starting points")

  P <- length(par)
  iter <- 1
  p.now <- par
  H.old <- -diag(P)
  H.new <- H.old
  U <- matrix(0, nrow = P, ncol = qn)
  W <- matrix(0, nrow = P, ncol = qn)


  lold <- objfn(p.now, ...)
  leval <- 1
  if (verbose) cat("Objective fn: ", lold, "\n")
  feval <- 0
  conv <- TRUE
  p.inter <- c(p.now, lold)

  for (i in 1:qn) {
    p.new <- try(fixptfn(p.now, ...),silent=TRUE)
    U[,i] <- p.new - p.now
    p.now <- p.new
    l.old <- (objfn(p.now, ...))
    if(intermed) p.inter <- rbind(p.inter, c(p.now, lold))
  }

  feval <- feval + qn
  leval <- leval + qn

  if(qn > 1)
    W[,1:(qn-1)] <- U[,2:qn] - U[,1:(qn-1)]

  p.new <- try(fixptfn(p.now, ...),silent=TRUE)
  W[,qn] <- p.new - p.now - U[,qn]
  feval <- feval + 1
  iter <- iter + qn
  old_secant <- 1

  while (iter < maxiter) {


    if(iter %% 100 ==0)  print(iter)
    extrap <- TRUE
    p1 <- try(fixptfn(p.now, ...),silent=TRUE)
    if (inherits(p1, "try-error") | any(is.nan(unlist(p1)))) stop("Error in function evaluation")
    u <- p1 - p.now
    U[,old_secant] <- u
    sr2 <- crossprod(u)
    if (sqrt(sr2) < tol) break

    p2 <- try(fixptfn(p1, ...),silent=TRUE)
    feval <- feval + 2
    if (inherits(p2, "try-error") | any(is.nan(unlist(p2)))) stop("Error in function evaluation")

    q2 <- p2 - p1
    sq2 <- (crossprod(q2))
    if (sqrt(sq2) < tol) break
    v <- q2 - u
    sv2 <- crossprod(v)
    W[,old_secant] <- v
    old_secant <- (old_secant %% qn) + 1

    prod1 <- solve(t(W) %*% W) %*% t(W)
    H.new <- H.old - (H.old %*% W %*% prod1) + U %*% prod1
    H.old <- H.new


    dir <- H.old %*% u
    step <- max(step.min, min(step.max, sqrt(sr2/sv2)))
    p.new <- p.now - step * as.numeric(sqrt(sr2)) * (dir/norm(dir, type='2'))

    if (inherits(p.new, "try-error") | any(is.nan(p.new))) {
      p.new <- p2
      lnew <- try(objfn(p2, ...), silent=TRUE)
      leval <- leval + 1
      extrap <- FALSE
    } else {
      if (is.finite(objfn.inc)) {
        lnew <- try(objfn(p.new, ...), silent=TRUE)
        leval <- leval + 1
      } else lnew <- lold

      if (inherits(lnew, "try-error") | is.nan(lnew) |  (lnew > lold + objfn.inc)) {
        if(verbose) print(paste("Fallback by: ", lnew - lold))
        p.new <- p2
        lnew <- try(objfn(p2, ...), silent=TRUE)
        leval <- leval + 1
        extrap <- FALSE
      }
    }

    if (obj.stop)
      if(abs(lnew - lold) <= obj.tol) break

    p.now <- p.new
    if(iter %% 100 == 0)print(lnew)
    if (!is.nan(lnew)) lold <- lnew
    if (verbose) cat("Objective fn: ", lnew, "  Extrapolation: ", extrap,  "\n")
    if(intermed) p.inter <- rbind(p.inter, c(p.now, lnew))

    iter <- iter+1
  }
  if (feval >= maxiter) conv <- FALSE
  if (is.infinite(objfn.inc)) {
    lold <- objfn(p.now, ...)
    leval <- leval + 1
  }

  rownames(p.inter) <- NULL
  if (!intermed) return(list(par=p.now, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv))
  else return(list(par=p.now, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv, p.inter=p.inter))
}

##########################################################
## Limited memory BQN method
#########################################################


LBQN <- function(par, fixptfn, objfn, ... , control=list()) {

  control.default <- list(m=10, objfn.inc=1, tol=1e-07, obj.tol=1e-07, maxiter=1500, verbose=FALSE, obj.stop=FALSE, intermed=FALSE)

  namc <- names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  ctrl <- modifyList(control.default, control)

  m <- ctrl$m
  maxiter <- ctrl$maxiter
  tol <- ctrl$tol
  objfn.inc <- ctrl$objfn.inc
  verbose <- ctrl$verbose
  intermed <- ctrl$intermed
  obj.tol <- ctrl$obj.tol
  obj.stop <- ctrl$obj.stop

  if (verbose) cat("LBQN \n")
  if (missing(objfn)) stop("\n objective function is not available \n\n")
  if (missing(par)) stop("\n Starting vector not available \n")
  P <- length(par)
  iter <- 1
  p.now <- par
  m.u <- list()
  m.v <- list()
  lold <- objfn(p.now, ...)
  leval <- 1
  if (verbose) cat("Objective fn: ", lold, "\n")
  feval <- 0
  conv <- TRUE
  p.inter <- c(p.now, lold)

  while (iter < maxiter) {

    extrap <- TRUE
    p1 <- try(fixptfn(p.now, ...),silent=TRUE)
    feval <- feval + 1
    if (inherits(p1, "try-error") | any(is.nan(unlist(p1)))) stop("Error in function evaluation")
    u <- p1 - p.now
    sr2 <- crossprod(u)
    if (sqrt(sr2) < tol) break

    p2 <- try(fixptfn(p1, ...),silent=TRUE)
    feval <- feval + 1
    if (inherits(p2, "try-error") | any(is.nan(unlist(p2)))) stop("Error in function evaluation")

    q2 <- p2 - p1
    sq2 <- sqrt(crossprod(q2))
    if (sq2 < tol) break
    v <- q2-u

    m.u[[1]] <- u
    m.v[[1]] <- v


    gamma_t <- as.numeric(crossprod(u, v))/as.numeric(crossprod(v, v))
    H_init <- gamma_t*diag(P)
    q <- u
    alpha <- rep(0,min(m, iter-1))

    if (iter >= 2){
      for (i in 1:min(m, (iter-1))){
        rho <- 1/as.numeric(crossprod(m.v[[i]], m.v[[i]]))
        alpha[i] <- rho*crossprod(m.v[[i]], q)
        q <- q - alpha[i]*m.v[[i]]
      }
    }

    r <- H_init %*% q
    if(iter >= 2){
      for (i in 1:min(m, iter-1)){
        r <- r + alpha[i]*m.u[[i]]
      }
    }

    p.new <- p.now - r


    if (inherits(p.new, "try-error") | any(is.nan(p.new))) {
      p.new <- p2
      lnew <- try(objfn(p2, ...), silent=TRUE)
      leval <- leval + 1
      extrap <- FALSE
    } else {
      if (is.finite(objfn.inc)) {
        lnew <- try(objfn(p.new, ...), silent=TRUE)
        leval <- leval + 1
      } else lnew <- lold

      if (inherits(lnew, "try-error") | is.nan(lnew) | (lnew > lold + objfn.inc)) {
        if (verbose) print(paste("Fallback by:", lnew - lold))
        p.new <- p2
        lnew <- try(objfn(p2, ...), silent=TRUE)
        leval <- leval + 1
        extrap <- FALSE
      }
    }

    if(obj.stop)
      if(abs(lnew - lold) < obj.tol) break

    p.now <- p.new

    if (!is.nan(lnew)) lold <- lnew
    if (verbose) cat("Objective fn: ", lnew, "  Extrapolation: ", extrap, "\n")
    if(intermed) p.inter <- rbind(p.inter, c(p.now, lnew))

    if(iter >=2){
      for (i in min(iter,m):2){
        m.u[[i]] <- m.u[[i-1]]
        m.v[[i]] <- m.v[[i-1]]
      }
    }

    iter <- iter+1
  }
  if (feval >= maxiter) conv <- FALSE
  if (is.infinite(objfn.inc)) {
    lold <- objfn(p.now, ...)
    leval <- leval + 1
  }

  rownames(p.inter) <- NULL
  if (!intermed) return(list(par=p.now, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv))
  else return(list(par=p.now, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv, p.inter=p.inter))
}

