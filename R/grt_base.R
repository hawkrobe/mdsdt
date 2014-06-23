source("two_by_two_fit.R")

#  Revision to fit model using Newton-Raphson iteration.
#    Old version uses fit.grt.old

# Define a grt object.
#  dists: Matrix giving means, standard deviations, and correlation
#           for each stimuls type
#  fit:   Optional information from fitting program
#    obs:      Observed frequencies
#    fitted:   Fitted frequencies
#    estimate: Estimated parameter vector
#    expd2:    Hessian (expected second derivative) matrix at estimate
#    map:      Parameter map
#    loglik:   Log likelihood at estimate
#    code:     Convergence code (as nlm)
#    iter:     Number of iterations required
#  rcuts: Optional cutpoints for the rows
#  ccuts: Optional cutpoints for the columns
grt <- function (dists, fit=NULL, rcuts = 0, ccuts = 0) {
  if (is.null(colnames(dists)))
    colnames(dists) <- c('mu_r','sd_r','mu_c','sd_c','rho')
  structure(list(dists=dists,fit=fit, rowcuts=rcuts, colcuts=ccuts),
            class = 'grt')
}

fit.grt <- function(freq, PS_x = FALSE, PS_y = FALSE, PI = 'none') {
  if (length(freq) == 16) {
    two_by_twomod <- create_two_by_two_mod(PS_x, PS_y, PI);
    return(two_by_twofit.grt(freq, two_by_twomod))
  } else {
    n_by_nmap <- create_n_by_n_mod(PS_x, PS_y, PI);
    print(n_by_nmap);
    return(n_by_n_fit.grt(freq, pmap=n_by_nmap))
  }
}

# Prints a grt object, but not fitting information
print.grt <- function (b) {
  cat('Row cuts:',format(b$rowcuts,digits=3),'\n')
  cat('Col cuts:',format(b$colcuts,digits=3),'\n')
  cat('Distributions:\n')
  print(round(b$dists,3))
  invisible(b)
}

# Summary prints the object and, if available, information about
# standard derrors and the fit.
summary.grt <- function(b) {
  print.grt(b)
  if (!is.null(fit <- b$fit)){
    cat('Standard errors:\n')
    print(round(distribution.se(b),3))
    print(GOF(b,test='X2'))
    print(GOF(b,test='G2'))
    cat('Log likelihood',fit$loglik,'\nConvergence code',fit$code,
        'in',fit$iter,'iterations\n')
  }
  invisible(b)
}


# Plots a grt object.  Nonobvious arguments are
#   bb:      A grt object.
#   level:   Proportion of distribution enclosed in ellipse (default 0.5).
#   connect: Connect points in the order given in a vector.
#   llty, lcol: Line type and color for cutpoint lines.
#   clty, ccol: Line type and color for connecting line.
#   names:   Print distribution names in center of distributions.
# Standard commands for axes, overall label, etc. are also available.  Unless
# otherwise specified, axis labels are taken from the dimension names.
plot.grt <- function(bb,level=.5,xlab=NULL,ylab=NULL,lim.sc=2, # lim.sc added 1.24.14
                        connect=NULL,names=NULL,clty=1,ccol='Black',llty=1,lcol='Black',
                        main=deparse(substitute(bb)),...) {
  require(ellipse)
  xc <- bb$colcuts
  yc <- bb$rowcuts
  dd <- bb$dists
  mx <- dd[,3]; my <- dd[,1]
  sx <- dd[,4]; sy <- dd[,2]
  rho <- dd[,5]
  min.ax <- min(min(mx-lim.sc*sx),min(my-lim.sc*sy))
  max.ax <- max(max(mx+lim.sc*sx),max(my+lim.sc*sy))
  X <- c(min.ax,max.ax)
  Y <- c(min.ax,max.ax)
  if (is.null(xlab)) xlab <- if(is.null(bb$fit)) 'Y' else
    names(dimnames(bb$fit$obs)[2])
  if (is.null(ylab)) ylab <- if(is.null(bb$fit)) 'X' else
    names(dimnames(bb$fit$obs)[1])
  # axes=F, box(which="plot") added 1.24.14
  plot(X,Y,type='n',main=main,xlab="",ylab="",axes=F,...)
  mtext(text=xlab,side=1,line=1)
  mtext(text=ylab,side=2,line=1)
  box(which="plot")
  for (i in 1:length(yc)) abline(h=yc[i],lty=llty,col=lcol)
  for (i in 1:length(xc)) abline(v=xc[i],lty=llty,col=lcol)
  for (i in 1:dim(dd)[1]) {
    v <- matrix(c(sx[i]^2,rep(sx[i]*sy[i]*rho[i],2),sy[i]^2),2)
    lines(ellipse(v,centre=c(mx[i],my[i]),level=level))
  } 
  if (!is.null(connect[1])) 
    lines(mx[c(connect,connect[1])],my[c(connect,connect[1])],
          lty=clty,col=ccol)
  # names changed 1.24.14
  if (!is.null(names)) text(mx,my,names)
}


# Overall fitting function.  Fitting either uses Newton-Raphson iteration
# or one of the method provided by the R function constrOptim 
# xx: The frequency table.  It can be entered in two ways
#     1)  A three-way 'xtabs' table with the stimulis as the third index
#     2)  A data frame contiaing the three indices with condition last,
#         and frequencies as the variable 'x' (see 'form' if not this way)
# pmap:     Parameter-mapping array (default: complete parameteriation)
# form:     A formula to convert a data frame (default x~.)
# p0:       Initial point---calculated if not given
# verbose:  Whether to print results (default FALSE)
# trace:    Print steps in minimization (default FALSE)
# method:   NA for method of scoring (default) or a method used by
#           constrOptim
# maxiter:  Maximum number of iterations
# stepsize: Size of iteration step (usually 1) in method of scoring
# maxch:    Maximum change for Newton Raphson convergence
# ...:      Arguments passed to minimization or likelihood routines
# Returns a grt object
n_by_n_fit.grt <- function (xx, pmap=NA, formula=x~., p0=NA, method=NA,
                        verbose=FALSE, trace=FALSE, maxiter=100, stepsize=1.0, maxch=1.0e-5, ...) {
  if (identical(class(xx)[1],'data.frame')) xx <- xtabs(formula,xx)
  if (!any(class(xx) == 'table'))
    stop('First argument must be data frame or contingency table')
  #  if (!(is.na(method) || (method == 0) || (method == 1)))
  #     stop('Method must be NA, 0, or 1')
  dxx <- dim(xx)
  if (length(dxx) != 3) stop('Table must have dimension 3')
  KK <- dxx[3];
  if (is.na(pmap)[1]) pmap <- matrix(c(rep(0:(KK-1),4),1:KK),4)
  colnames(pmap) <- c('mu','sigma','nu','tau','rho')
  if (is.null(rownames(pmap))) rownames(pmap) <- dimnames(xx)[[3]]
  bsd.valid.map(pmap,KK)
  # Create initial vector if required
  if (is.na(p0[1])) p0 <- bsd.initial(xx,pmap)
  # Construct index arrays
  imap <- bsd.imap(pmap,dxx)
  if (verbose) {
    cat('Parameter mapping vector\n'); print(pmap)
    cat('Initial parameter vector\n'); # print(round(p0,3))
    cat('Row cutpoints', round(p0[imap$xi],4),'\n')
    cat('Col cutpoints', round(p0[imap$eta],4),'\n')
    cat('Parameters by groups\n')
    print(round(bsd.map2array(p0,pmap,imap),4))
  }
  # Do the minimization
  if (is.na(method)){
    if (verbose) cat('Fitting by Newton-Raphson iteration\n')
    found <- FALSE
    pold <- p0
    for (iter in 1:maxiter) {
      #print(pold)
      #print(pmap)
      #print(imap)
      bb <- bsd.llike(pold,xx,pmap,imap,d.level=2,...)
      #print(attr(bb,'ExpD2'))
      #print(attr(bb,'gradient'))
      dlt <- solve(attr(bb,'ExpD2'),attr(bb,'gradient'))
      if (trace){
        cat('Iteration number',iter,'\n')
        cat('Row cutpoints', round(pold[imap$xi],4),'\n')
        cat('Col cutpoints', round(pold[imap$eta],4),'\n')
        print(round(bsd.map2array(pold,pmap,imap),4))
        cat('Value',bb[1],'\n')
      }
      s <- stepsize
      repeat{
        pp <- pold + s*dlt
        if (all(c(diff(pp[imap$xi]),diff(pp[imap$eta]))>0)) break
        s <- s/2
        warning('Reduced stepsize to',s,call.=FALSE)
        if (s < 0.001)
          stop('Stepsize reduction too large: check problem definition')
      }
      if (max(abs(dlt)) < maxch) {
        found <- TRUE
        iterations <- iter
        break
      }
      else iterations <- maxiter
      pold <- pp
    }
    code <- if (found) 0 else 4
  }
  else {
    if (verbose) cat('Fitting using optim\n')
    ixxi <- imap$xi; lxi <- length(ixxi)
    ixeta <- imap$eta; leta <- length(ixeta)
    ixvar <- c(imap$sigma,imap$tau) ; lvar <- length(ixvar)
    ixrho <- imap$rho; lrho <- length(ixrho)
    lp <- length(p0)
    if(length(ixxi)>1 && length(ixeta)>1){
      nconst <- lxi+leta+lvar+2*lrho-2
    }else{
      nconst <- 2*lrho
      b <- 1
    }
    cm <- matrix(0,nrow=nconst,ncol=lp)
    if(length(ixxi)>1){
      for(i in 1:(lxi-1)) cm[i,ixxi[i:(i+1)]] <- c(-1,1)
      b <- lxi-1
    }
    if(length(ixeta)>1){
      for(j in 1:(leta-1)) cm[j+b,ixeta[j:(j+1)]] <- c(-1,1)
      b <- b + leta-1
    }
    if (lvar > 0){
      for (i in 1:lvar) cm[i+b,ixvar[i]] <- 1
      b <- b + lvar
    }
    if (lrho > 0) for (i in 1:lrho){
      if(length(ixxi)>1 && length(ixeta)>1){
        cm[(b+1):(b+2),ixrho[i]] <- c(1,-1)
        b <- b+2
      }else{
        cm[b:(b+1),ixrho[i]] <- c(1,-1)
        b <- b+2
      }
    }
    cv <- c(rep(0,nconst-2*lrho),rep(-1,2*lrho))
    ffit <- constrOptim(p0,bsd.like,bsd.grad,cm,cv,method=method,
                        x=xx, pmap=pmap, imap=imap,...)
    pp <- ffit$par
    code <- ffit$convergence
    iterations <- ffit$counts
    bb <- bsd.llike(pp,xx,pmap,imap,d.level=2)
    found <- code < 4
  }
  if (!found) warning('Minimization routine encountered difficultites',
                      call.=FALSE)
  # Assemble results and print if required
  names(pp) <- names(p0)
  xi <- pp[imap$xi]
  eta <- pp[imap$eta]
  dists <- bsd.map2array(pp,pmap,imap)
  nk <- apply(xx,3,sum)
  muh <- array(dim=dim(xx),dimnames=dimnames(xx))
  for (k in 1:KK) muh[,,k] <- bsd.freq(xi,eta,dists[k,],nk[k])
  fit <- list(obs=xx,fitted=muh,estimate=pp,
              expd2=attr(bb,'ExpD2'),map=pmap,
              loglik=bb[1],code=code,iter=iterations)
  if (verbose) {
    if (found) cat('\nConvergence required',iterations,'iterations\n\n')
    else cat (iterations, 'iterations used without reaching convergence\n\n')
    cat('Parameter estimate vector\n'); # print(round(pp,3))
    cat('Row cutpoints', round(xi,4),'\n')
    cat('Col cutpoints', round(eta,4),'\n')
    print(round(dists,4))
    cat('Log likelihood =',bb[1],'\n')
  }
  print(dists);  
  return(grt(dists,fit=fit,rcuts = xi,ccuts = eta));
}

create_n_by_n_mod <- function(PS_x, PS_y, PI) {
  # Each row is distribution, cols are x_mean, x_std, y_mean, y_std, rho
  map <- matrix(data = 0, nrow = 4, ncol= 5)
  if (PS_x) { 
    map[2,1:2] <- map[4,1:2] <- c(1,1); 
  } else for (i in 1:4) map[i,1:2] <- c(i-1,i-1);
  if (PS_y) {
    map[3,3:4] <- map[4,3:4] <- c(1,1);
  } else for (i in 1:4) map[i,3:4] <- c(i-1,i-1);
  if (PI == 'same_rho') {
    for (i in 1:4) map[i,5] <- 1;
  } else if (PI == 'none') {
    for (i in 1:4) map[i,5] <- i;
  } 
  return(map);
}


# Get parameter map
parameter.map <- function(bb){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  if (is.null(ff <- bb$fit)) null else ff$map
}


# Estimated parameters
coef.grt <- function(bb)
  if (is.null(ff <- bb$fit)) NULL else ff$estimate


# Covariance matrix of parameters
vcov.grt <- function(bb){
  if (is.null(ff <- bb$fit)) return(NULL)
  vcv <- solve(-ff$expd2)
  rownames(vcv) <- colnames(vcv) <- names(ff$estimate)
  vcv
}


# Parameters by stimuli
distribution.parameters <- function(bb){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  bb$dists
}


# Standard errors by stimuli
distribution.se <- function(bb){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  if (is.null(ff <- bb$fit)) return(NULL)
  pmap <- ff$map
  dimx <- dim(ff$obs)
  imap <- bsd.imap(pmap,dimx)
  sds <- bsd.map2array(sqrt(diag(vcov(bb))),pmap,imap,0,0)
  rownames(sds) <- rownames(bb$dists)
  colnames(sds) <- colnames(bb$dists)
  sds
}


# Log likelihood
logLik.grt <- function(bb){
  if (is.null(ff <- bb$fit)) return(NULL)
  sum(lfactorial(apply(ff$obs,3,sum))) - sum(lfactorial(ff$obs)) +
    ff$loglik
}


# Pearson residuals
residuals.grt <- function(bb){
  if (is.null(ff <- bb$fit)) return(NULL)
  xx <- ff$obs
  ex <- ff$fitted
  (xx-ex)/sqrt(ex)
}


# Fitted values
fitted.grt <- function(bb)
  if (is.null(ff <- bb$fit)) NULL else ff$fitted


# Goodness of fit
GOF <- function(bb,teststat='X2',observed=NULL){
  if (!identical(class(bb),'grt'))
    stop('Argument must be object of class "grt"')
  if (is.null(ff <- bb$fit) && is.null(observed))
    stop('Must have fitted model, observed frequencies, or both')
  # added AIC, AICc, BIC 1.27.14
  statlist <- c('X2','G2','AIC','AIC.c','BIC')
  #print(c(statlist,teststat))
  #print(statlist %in% teststat)
  #teststat <- toupper(teststat)
  test <- pmatch(teststat,statlist)
  if (is.na(test)) stop('Test statistic unrecognized')
  teststat <- statlist[test]
  df <- if (is.null(observed)) length(ff$estimate) else 0
  if (is.null(observed)) observed <- ff$obs
  if (!is.null(ff)) ex <- ff$fitted else{
    nk <- apply(observed,3,sum)
    ex <- array(dim=dim(observed))
    for (k in 1:dim(observed)[3])
      ex[,,k] <- bsd.freq(bb$rowcuts,bb$colcuts,bb$dists[k,],nk[k])
  }
  df <- length(observed) - dim(observed)[3] - df
  if (test == 1){
    tstat <- sum((observed-ex)^2/ex)}
  if (test == 2){
    ex <- ex[observed>0]; observed <- observed[observed>0]
    tstat <- 2*sum(observed*log(observed/ex))
  }
  if (test == 3){
    k <- nrow(bb$dists)
    tstat <- 2*bb$fit$loglik + 2*k
  }
  if (test == 4){
    k <- nrow(bb$dists)
    n <- sum(observed)
    tstat <- 2*bb$fit$loglik + (2*k*(k+1))/(n-k-1)
  }
  if (test == 5){
    k <- nrow(bb$dists)
    n <- sum(observed)
    tstat <- 2*bb$fit$loglik + log(n)*k
  }
  if (test < 3){
    structure(tstat,names=teststat,df=df,class='bsdGOF')    
  }else{
    structure(round(tstat,1),names=teststat)
  }
  
}

print.bsdGOF <- function(gof){
  df <- attr(gof,'df')
  cat(names(gof),'(',df,') = ',gof,', p = ',
      round(pchisq(gof,df,lower.tail=F),5), '\n',sep='')
}


# Wald test of a linear hypothesis m p = c 
#   b:   fitted grt model containing estimates of p
#   m:   contrast vector or matrix with contrast vectors as rows
#   c:   a vector of numerical values (default zero)
#   set: set of parameters to test (means, sds, correlations,
#           distribution parameters, cutpoints, or all parameters)
linear.hypothesis <- function(b,m,c=0,set='means'){
  if (is.null(ff <- b$fit)) stop('Must test a fitted model')
  imap <- bsd.imap(ff$map,dim(ff$obs))
  if(is.na(set <- pmatch(set,
                         c('means','sds','correlations','distributions','cutpoints','all'))))
    stop('Test set unrecognized')
  set <- switch(set,
                '1'=c(imap$mu,imap$nu),
                '2'=c(imap$sigma,imap$tau),
                '3'=imap$rho,
                '4'=-c(imap$xi,imap$eta),
                '5'=c(imap$xi,imap$eta),
                '6'=1:(length(ff$estimate)))
  p <- ff$estimate[set]
  varp <- vcov(b)[set,set]
  if (!is.matrix(m)) m <- t(m)
  if (length(p) != dim(m)[2])
    stop('Size of hypothesis matrix must agree with number of parameters')
  df <- dim(m)[1]
  h <- m %*% p
  v <- m %*%varp %*% t(m)
  W <- t(h) %*% solve(v)%*% h
  structure(W,names='W',df=df,class='bsdGOF')
}


# Comparison of models
anova.grt <- function(b1,b2){
  g21 <- GOF(b1,test='G2')
  df1 <- attr(g21,'df')
  g22 <- GOF(b2,test='G2')
  df2 <- attr(g22,'df')
  # p.val added 1.24.14 -NHS
  DG2 <- round(g21-g22,3)
  ddf <- df1-df2
  p.val <- round(pchisq(DG2,ddf,lower.tail=F),4)
  table <- matrix(c(round(g21,3),round(g22,3),df1,df2,'',
                    DG2,'',ddf,'',p.val),2)
  dimnames(table) <- list(c(substitute(b1),substitute(b2)),
                          c('G2', 'df', 'DG2', 'df','p-val'))
  as.table(table)
}


# Test parameters for equality
test.parameters <- function(bb,set='means'){
  if (is.null(ff <- bb$fit)) stop('Must work with fitted model')
  imap <- bsd.imap(ff$map,dim(ff$obs))
  if(is.na(set <- pmatch(set, c('means','sds','correlations','cutpoints'))))
    stop('Test set unrecognized')
  c0 <- if (set==2) 1 else 0
  set <- switch(set,
                '1'=c(imap$mu,imap$nu),
                '2'=c(imap$sigma,imap$tau),
                '3'=imap$rho,
                '4'=c(imap$xi,imap$eta))
  p <- ff$estimate[set]
  vp <- vcov(bb)[set,set]
  ix2 <- lower.tri(vp)
  ix1 <- col(vp)[ix2]
  ix2 <- row(vp)[ix2]
  tv <- c(p-c0,p[ix1]-p[ix2])
  dv <- diag(vp)
  se <- sqrt(c(dv,dv[ix1]+dv[ix2]-2*vp[cbind(ix1,ix2)]))
  z <- tv/se
  structure(cbind(tv,se,z,2*pnorm(abs(z),lower.tail=FALSE)),
            dimnames=list(c(paste(names(p),'-',c0),
                            paste(names(p)[ix1],'-',names(p)[ix2])),
                          c('Est','s.e.','z','p')),
            class=c('bsd.test','matrix'))
}

print.bsd.test <- function(mx,digits=3){
  class(mx) <- 'matrix'
  print(round(mx,digits))
}

# Likelihood calculation routines

# Calculate minus the log-likelihood for a bivariate detection model
# pv:      Parameter vector
# x:       Data as a table
# pmap:    A 5 by KK array pmam mapping tables onto parameter vectors
# imap:    List of parameter indices by type in parameter vector
# d.level: Derivitive level to return
#           0: Likelihood only
#           1: Likelihood and first derivative
#           2: Likelihood, first derivative, and expected second derivatives
# diag:    Diagnostic print code: 0, 1, 2, or 3
# fracw:   Value use to space overlapping criteria
bsd.llike <- function (pv,x,pmap,imap,d.level=0,diag=0,fracw=10) {
  tt <- dim(x); II <- tt[1]; JJ <- tt[2]; KK <- tt[3]
  lpv <- length(pv)
  llike <- 0
  if (d.level > 0){
    gradient <- numeric(lpv)
    grx <- numeric(lpv-II-JJ+2)
    if (d.level == 2){
      ExpD2 <- matrix(0,lpv,lpv)
      d2eta <- 1 + (lpv+1)*(0:(II+JJ-1))
      d2xi <- d2eta[imap$xi];   d2xi1 <- (d2xi-1)[-1]
      d2eta <- d2eta[imap$eta]; d2eta1 <- (d2eta-1)[-1]
      nk <- apply(x,3,sum)
    }
  }
  aa <- bsd.map2array(pv,pmap,imap)
  xi <- pv[imap$xi];
  eta <- pv[imap$eta];
  if (diag > 0){
    cat('Call of bsd.llike\n')
    cat('xi ',xi,'\neta',eta,'\n')
    print(aa)
  }
  # Check the validity of parameters
  if (any(c(diff(xi),diff(eta),pv[imap$sigma],pv[imap$tau])<0)
      || any(abs(pv[imap$rho])>1)) {
    warning('Criteria out of order or invalid distribution',
            call.=FALSE)
    llike <- Inf
    if (d.level > 0) attr(llike,'gradient') <- gradient
    return(llike)
  }
  # Loop over tables
  for (k in 1:KK) {
    bf <- bsd.freq(xi,eta,aa[k,],if(d.level==0) 1 else NULL)
    gamma <- if (d.level==0) bf else bf$pi
    llike <- llike - sum(x[,,k]*log(gamma))
    if (diag > 1) {
      cat('Table level ',k,'\nProbabilities\n')
      if(is.list(bf)){
        print(bf$pi)  
      }
      cat('Likelihood contribution ',-sum(x[,,k]*log(gamma)),
          '\nCumulated log-likelihood',llike,'\n')
    }
    # Calculate gradients if required
    if (d.level > 0){
      xg <- x[,,k]/gamma
      t1 <- rbind(bf$d1,0)
      t2 <- rbind(bf$d2,0)
      if (diag > 1) {print(bf)
                     cat('x_ijk/pi_ijk\n'); print(xg)
      }
      vxi <-  D2(t1)[-II,]/aa[k,2]
      veta <- D2(t2)[-JJ,]/aa[k,4]
      if (diag > 2) {
        cat('vxi\n');  print(vxi)
        cat('veta\n'); print(veta)
      }
      # NHS: added 'length != 1' bit; not sure if it's working right...
      if(length(xi)!=1 && length(eta)!=1){
        gd <- c(rowSums(vxi*(xg[-II,]-xg[-1,])),
                rowSums(veta*t(xg[,-JJ]-xg[,-1])), grx)
      }else{
        gd <- c(sum(vxi*(xg[-II,]-xg[-1,])),
                sum(veta*t(xg[,-JJ]-xg[,-1])), grx)
      }
      ipm <- pmap[k,1]
      ips <- pmap[k,2]
      ipn <- pmap[k,3]
      ipt <- pmap[k,4]
      ipr <- pmap[k,5]
      if (ipm > 0){
        ixa <- imap$mu[ipm]
        vmu <- -D12(t1)/aa[k,2]
        if (diag > 2) {cat('vmu\n'); print(vmu)}
        gd[ixa] <- sum(xg*vmu)
      }
      if (ips > 0) {
        ixb <- imap$sigma[ips]
        vsigma <- -D12(((c(xi,0)-aa[k,1])/aa[k,2]^2)*t1)
        if (diag > 2) {cat('vsigma\n'); print(vsigma)}
        gd[ixb] <- sum(xg*vsigma)
      }
      if (ipn > 0) {
        ixk <- imap$nu[ipn]
        vnu <- t(-D12(t2))/aa[k,4]
        if (diag > 2) {cat('vnu\n'); print(vnu)}
        gd[ixk] <- sum(xg*vnu)
      }
      if (ipt > 0) {
        ixl <- imap$tau[ipt]
        vtau <- t(-D12(((c(eta,0)-aa[k,3])/aa[k,4]^2)*t2))
        if (diag > 2) {cat('vtau\n'); print(vtau)}
        gd[ixl] <- sum(xg*vtau)
      }
      if (ipr > 0) { 
        ixr <- imap$rho[ipr]
        vrho <- D12(rbind(cbind(bf$phi,0),0))
        if (diag > 2) {cat('vrho\n'); print(vrho)}
        gd[ixr] <- sum(xg*vrho)
      }
      gradient <- gradient + gd
      if (diag > 0) {
        cat('First derivative contribution\n')
        names(gd) <- names(pv)
        print(round(gd,4))
      }
      # Calculate expected second derivatives if requested
      # if(length(xi)!=1) and if(length(eta)!=1) clauses added 1.25.14 -NHS
      if (d.level == 2){
        oog <- 1/gamma
        if(length(eta)!=1){
          veta <- t(veta)          
        }
        hs <- matrix(0,lpv,lpv)
        
        if(length(xi)!=1){
          hs[d2xi] <- rowSums(vxi^2*(oog[-II,]+oog[-1,]))
          tt <- vxi[-(II-1),]*vxi[-1,]*oog[2:(II-1),]
        }else{
          hs[d2xi] <- sum(vxi^2*(oog[-II,]+oog[-1,]))
          tt <- vxi[-(II-1)]*vxi[-1]*oog[2:(II-1),]
        }
        hs[d2xi1] <- -if (is.matrix(tt)) rowSums(tt) else sum(tt)
        if(length(eta)!=1){
          hs[d2eta] <- colSums(veta^2*(oog[,-JJ]+oog[,-1]))
          tt <- veta[,-(JJ-1)]*veta[,-1]*oog[,2:(JJ-1)]          
        }else{
          hs[d2eta] <- sum(veta^2*(oog[,-JJ]+oog[,-1]))
          tt <- veta[-(JJ-1)]*veta[-1]*oog[,2:(JJ-1)]
        }
        hs[d2eta1] <- -if (is.matrix(tt)) colSums(tt) else sum(tt)
        
        if(length(xi)!=1 && length(eta)!=1){
          hs[imap$xi,imap$eta] <- 
            vxi[,-JJ]*(veta[-II,]*oog[-II,-JJ] - veta[-1,]*oog[-1,-JJ]) -
            vxi[,-1]*(veta[-II,]*oog[-II,-1] - veta[-1,]*oog[-1,-1])          
        }else{
          hs[imap$xi,imap$eta] <- 
            vxi[-JJ]*(veta[-II]*oog[-II,-JJ] - veta[-1]*oog[-1,-JJ]) -
            vxi[-1]*(veta[-II]*oog[-II,-1] - veta[-1]*oog[-1,-1])
        }
        if (ipm > 0){
          tt <- vmu*oog
          hs[ixa,ixa] <- sum(vmu*tt)
          if(length(xi)!=1){
            hs[imap$xi,ixa] <- rowSums(vxi*(tt[-II,] - tt[-1,]))
          }else{
            hs[imap$xi,ixa] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixa] <- colSums(veta*(tt[,-JJ] - tt[,-1]))            
          }else{
            hs[imap$eta,ixa] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
          if (ips > 0) hs[ixa,ixb] <- sum(tt*vsigma)
          if (ipn > 0) hs[ixa,ixk] <- sum(tt*vnu)
          if (ipt > 0) hs[ixa,ixl] <- sum(tt*vtau)
          if (ipr > 0) hs[ixa,ixr] <- sum(tt*vrho)
        }
        if (ips > 0){
          tt <- vsigma*oog
          hs[ixb,ixb] <- sum(tt*vsigma)
          hs[imap$xi,ixb] <- rowSums(vxi*(tt[-II,] - tt[-1,]))
          hs[imap$eta,ixb] <- colSums(veta*(tt[,-JJ] - tt[,-1]))
          if (ipn > 0) hs[ixb,ixk] <- sum(tt*vnu)
          if (ipt > 0) hs[ixb,ixl] <- sum(tt*vtau)
          if (ipr > 0) hs[ixb,ixr] <- sum(tt*vrho)
        }
        if (ipn > 0){
          tt <- vnu*oog
          hs[ixk,ixk] <- sum(tt*vnu)
          if(length(xi)!=1){
            hs[imap$xi,ixk] <- rowSums(vxi*(tt[-II,] - tt[-1,]))            
          }else{
            hs[imap$xi,ixk] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixk] <- colSums(veta*(tt[,-JJ] - tt[,-1])) 
          }else{
            hs[imap$eta,ixk] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
          if (ipt > 0) hs[ixk,ixl] <- sum(tt*vtau)
          if (ipr > 0) hs[ixk,ixr] <- sum(tt*vrho)
        }
        if (ipt > 0){
          tt <- vtau*oog
          hs[ixl,ixl] <- sum(tt*vtau)
          if(length(xi)!=1){
            hs[imap$xi,ixl] <- rowSums(vxi*(tt[-II,] - tt[-1,]))            
          }else{
            hs[imap$xi,ixl] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixl] <- colSums(veta*(tt[,-JJ] - tt[,-1])) 
          }else{
            hs[imap$eta,ixl] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
          if (ipr > 0) hs[ixl,ixr] <- sum(tt*vrho)
        }
        if (ipr > 0){
          tt <- vrho*oog
          hs[ixr,ixr] <- sum(tt*vrho)
          if(length(xi)!=1){
            hs[imap$xi,ixr] <- rowSums(vxi*(tt[-II,] - tt[-1,]))            
          }else{
            hs[imap$xi,ixr] <- sum(vxi*(tt[-II,] - tt[-1,]))
          }
          if(length(eta)!=1){
            hs[imap$eta,ixr] <- colSums(veta*(tt[,-JJ] - tt[,-1]))            
          }else{
            hs[imap$eta,ixr] <- sum(veta*(tt[,-JJ] - tt[,-1]))
          }
        }
        ExpD2 <- ExpD2 - nk[k]*hs
        if (diag > 2) {
          cat('Expected second derivative contribution\n')
          print(round(-nk[k]*hs,3))
        }
      }
    }
  }
  if (diag > 0) cat('Log-likelihood:',llike,'\n')
  if (d.level > 0){
    attr(llike,'gradient') <- -gradient
    if (d.level == 2) {
      hs <- ExpD2
      diag(hs) <- 0
      attr(llike,'ExpD2') <- ExpD2+t(hs)
    }
  }
  llike
}


# Differencing of first and second subscript of matrix
D1 <- function(x){
  dx <- dim(x)[1]
  rbind(x[1,],x[-1,]-x[-dx,])
}

D2 <- function(x){
  dx <- dim(x)[2]
  cbind(x[,1],x[,-1]-x[,-dx])
}

D12 <- function(x){
  r <- dim(x)[1]; c <- dim(x)[2]
  x <- rbind(x[1,],x[-1,]-x[-r,])
  cbind(x[,1],x[,-1]-x[,-c])
}

bsd.like <- function(p,...) bsd.llike(p,d.level=0,...)
bsd.grad <- function(p,...) attr(bsd.llike(p,d.level=1,...),'gradient')

# Calculate the cell frequencies for a bivariate SDT model.
# xi and eta: Row and column criteria
# m:          A vector of the distributional parameters
#               (mu_r, sigma_r, mu_c, sigma_c, rho).
# n:          Sample size or NULL
# When a sample size n is given, the function returns expected frequencies;
# when it is NULL, the function returns a list containing the probabilities pi,
# the densities phi, and the row and column derivative terms (the latter
# as its transpose).
bsd.freq <- function (xi,eta,m,n=NULL) {
  require(mvtnorm)
  fracw <- 10
  nrow <- length(xi) +1
  ncol <- length(eta)+ 1
  Xi  <- c(-Inf, (xi-m[1])/m[2], Inf)
  Eta <- c(-Inf, (eta-m[3])/m[4], Inf)
  rho <- m[5]
  pii <- matrix(nrow=nrow, ncol=ncol)
  cx <- matrix(c(1,rho,rho,1),2)
  for (i in 1:nrow) for (j in 1:ncol) 
    pii[i,j] <- pmvnorm(c(Xi[i],Eta[j]), c(Xi[i+1],Eta[j+1]),corr=cx)
  if (is.null(n)){
    Xis <- Xi[2:nrow]
    Etas <- Eta[2:ncol]
    phi <- matrix(0, nrow=nrow-1, ncol=ncol-1)
    for (i in 1:(nrow-1)) for (j in 1:(ncol-1))
      phi[i,j] <- dmvnorm(c(Xis[i],Etas[j]),sigma=cx)
    list(pi = pii, phi = phi,
         d1 = dnorm(Xis)*pnorm(outer(-rho*Xis,c(Etas,Inf),'+')/sqrt(1-rho^2)), 
         d2 = dnorm(Etas)*pnorm(outer(-rho*Etas,c(Xis,Inf),'+')/sqrt(1-rho^2)))
  }
  else pii*n
}


# Checks that a parameter map has the correct form
bsd.valid.map <- function(map,K){
  dm <- dim(map)
  if (!is.matrix(map)) stop('Map must be a matrix')
  if (dm[1] != K) stop('Map must have same number of rows as conditions')
  if (dm[2] != 5) stop('Map must have 5 columns')
  for (i in 1:4){
    u <- unique(map[,i])
    if (min(u) != 0) stop(paste('Values in map column',i,'must start at 0'))
    if (max(u) != length(u)-1) 
      stop(paste('Map column',i,'must be dense series'))
  }
  u <- unique(map[,5])
  nu <- min(u); xu <- max(u)
  if (((nu==0) && (xu!=length(u)-1)) || ((nu==1) & (xu!=length(u))))
    stop('Map column 5 must start at 0 or 1 and be dense series')
  TRUE
}


# Changes notation from cutpoints to differences of cutpoints and back
#bsd.todiff <- function(p,ixxi,ixeta){
#  c(p[ixxi[1]],diff(p[ixxi]),p[ixeta[1]],diff(p[ixeta]),
#     p[(max(ixeta)+1):length(p)])
#  }
#bsd.tocuts <- function(p,ixxi,ixeta){
#  c(cumsum(p[ixxi]),cumsum(p[ixeta]),p[(max(ixeta)+1):length(p)])
#  }


# Constructs the parameter indices
# pmap is parameter map and dimx is dimension of data
bsd.imap <- function(pmap,dimx){
  II <- dimx[1]; JJ <- dimx[2]
  ixxi <- 1:(II-1);  ib <- II+JJ-1; ixeta <- II:(ib-1);
  inp <- max(pmap[,1]); ixmu  <- ib:(ib+inp-1);     ib <- ib+inp
  # if-else clause added 1.20.14 -NHS
  if(max(pmap[,2])>0){
    inp <- max(pmap[,2]); ixsigma <- ib:(ib+inp-1); ib <- ib+inp    
  }else{
    ixsigma <- NULL#inp <- 1; ixsigma <- ib; ib <- ib+inp
  }
  inp <- max(pmap[,3]); ixnu  <- ib:(ib+inp-1);     ib <- ib+inp
  # if-else clause added 1.20.14 -NHS
  if(max(pmap[,4])>0){
    inp <- max(pmap[,4]); ixtau <- ib:(ib+inp-1); ib <- ib+inp    
  }else{
    ixtau <- NULL#inp <- 1; ixtau <- ib; ib <- ib+inp
  }
  # if-else clause added 1.20.14 -NHS
  if(max(pmap[,5])>0){
    inp <- max(pmap[,5]); ixrho <- ib:(ib+inp-1)    
  }else{
    ixrho <- NULL
  }
  list(xi=ixxi,eta=ixeta,mu=ixmu,sigma=ixsigma,nu=ixnu,
       tau=ixtau,rho=ixrho)
}


# Takes a parameter vector 'p' and a map of parameters are returns a
# table of parameter values.
# 'm0' and 'x0' are mean and variance values for parameter with map 0.
bsd.map2array <- function(p,pmap,imap,m0=0,s0=1){
  KK <- dim(pmap)[1]
  aa <- matrix(m0,KK,5)
  rownames(aa) <- rownames(pmap)
  colnames(aa) <- colnames(pmap)
  for (k in 1:KK) {
    if (pmap[k,1] != 0) aa[k,1] <- p[imap$mu[pmap[k,1]]]
    aa[k,2] <- if (pmap[k,2] != 0) p[imap$sigma[pmap[k,2]]] else s0
    if (pmap[k,3] != 0) aa[k,3] <- p[imap$nu[pmap[k,3]]]
    aa[k,4] <- if (pmap[k,4] != 0) p[imap$tau[pmap[k,4]]] else s0
    if (pmap[k,5] != 0) aa[k,5] <- p[imap$rho[pmap[k,5]]]
  }
  aa}


# Construct an initial vector
# The value delta is added to all frequencies to avoid problems with zeros
bsd.initial <- function(xx,pmap,delta=0.5){
  require(polycor)
  pnames <- c('mu','sigma','nu','tau','rho')
  dxx <- dim(xx)
  II <- dxx[1]; JJ <- dxx[2];  KK <- dxx[3];
  ixx <- 1:(II-1)
  ixy <- 1:(JJ-1)
  xx <- xx + delta
  ptx <- prop.table(margin.table(xx,c(3,1)),1)
  pty <- prop.table(margin.table(xx,c(3,2)),1)
  xi <- rep(0,II); eta <- rep(0,JJ)
  ni <- nj <- 0
  for (k in 1:KK){
    ptx[k,] <- qnorm(cumsum(ptx[k,]))
    if (pmap[k,1]+pmap[k,2]==0) {xi <- xi+ptx[k,]; ni <- ni+1}
    pty[k,] <- qnorm(cumsum(pty[k,]))
    if (pmap[k,3]+pmap[k,4]==0) {eta <- eta+pty[k,]; nj <- nj+1}
  }
  # NHS 1.20.14
  # added 'as.matrix' terms to maintain KK rows in ptx, pty
  # code still doesn't work with 2x2 data
  # si matrix (below) ends up with NAs in cols 2 and 4
  ptx <- as.matrix(ptx[,ixx]); pty <- as.matrix(pty[,ixy])
  xi <- (xi/ni)[ixx]; eta <- (eta/nj)[ixy]
  pv <- c(xi,eta)
  nv <- c(paste('xi',ixx,sep=''),paste('eta',ixy,sep=''))
  np <- apply(pmap,2,max)
  si <- matrix(0,KK,4)
  for (k in 1:KK){
    si[k,1:2] <- coef(lm(ptx[k,]~xi)) 
    si[k,3:4] <- coef(lm(pty[k,]~eta))
  }
  if(II>2){
    si[,1] <- -si[,1]/si[,2]
    si[,2] <- 1/si[,2]    
  }else{
    si[,1] <- -si[,1]
    si[,2] <- 1
  }
  if(JJ>2){
    si[,3] <- -si[,3]/si[,4]
    si[,4] <- 1/si[,4]    
  }else{
    si[,3] <- -si[,3]
    si[,4] <- 1
  }
  for (i in 1:4) if (np[i] > 0) {
    pv <- c(pv,tapply(si[,i],pmap[,i],mean)[-1])
    nv <- c(nv, paste(pnames[i],1:np[i],sep=''))
  }
  if (np[5]>0){
    r <- rep(0,KK)
    for (k in 1:KK) r[k] <- polychor(xx[,,k])
    rr <- tapply(r,pmap[,5],mean)
    if (min(pmap[,5]) == 0) rr <- rr[-1]
    pv <- c(pv,rr)
    nv <- c(nv, paste('rho',1:np[5],sep=''))
  }
  names(pv) <- nv
  pv
}


# Some test functions

test.bsd.llike <- function(pv,xx,pmap,n=1,d.level=0,diag=d.level+1){
  dxx <- dim(xx)
  imap <- bsd.imap(pmap,dxx)
  cat('Parameter mapping vector\n');   print(pmap)
  cat('Parameters by groups\n')
  print(bsd.map2array(pv,pmap,imap))
  bsd.llike(pv,xx,pmap,imap,d.level=d.level,diag=diag)
}

test.bsd.freq <- function(xi=c(-.6,.15,.65), eta=c(-.5,.25,.75),
                          m=c(0,1,.2,1.2,-.3),n=NULL) {
  print('Calling bsd.freq')
  bsd.freq(xi,eta,m,n)
}