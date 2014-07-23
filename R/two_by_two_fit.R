#' Tests sampling independence for each stimulus response distribution
#'
#' @param x four-by-four confusion matrix 
#' @return data frame containing z-scores and p-values for all four tests
#' @details If p value is sufficiently low, we're justified in rejecting the null hypothesis of sampling independence. 
#' @examples
#' data(thomasA)
#' siTest(thomasA)
#' @export
siTest <- function(x) {
  if(!checkConfusionMatrix(x)) {
    return(FALSE)
  }
  
  stimulus <- c("(1,1)", "(1,2)", "(2,1)", "(2,2)")
  statistic <- rep(NA,4)
  for ( i in 1:4 ) { 
    x1 <- matrix(x[i,], 2,2, byrow=T)
    ex1 <- c(apply(x1,1,sum)*apply(x1,2,sum),
             rev(apply(x1,1,sum)) * apply(x1,2,sum))
    ex1 <- matrix(ex1[c(1,4,3,2)],2,2,byrow=TRUE) / sum(x1)
    statistic[i] <- sum( (x1 - ex1)^2/ ex1 )
  }
  
  return(data.frame(stimulus=stimulus, statistic=statistic, 
                    p.value=1-pchisq(statistic, 1)) )
}

#' Tests marginal response invariance at both levels on each dimension
#'
#' @param x four-by-four confusion matrix 
#' @return data frame containing z-scores and p-values for all four tests
#' @details If the p value for either level of the x dimension is significant, 
#' we are justified in rejecting the null hypothesis of perceptual separability on the x dimension. 
#' Similarly for the y dimension.
#' @examples
#' data(thomasA)
#' mriTest(thomasA)
#' @export
mriTest <- function(x) {
  if(!checkConfusionMatrix(x)) {
    return(FALSE)
  }
  
  stimulus <- c("(1,-)", "(2,-)", "(-,1)", "(-,2)")
  statistic <- rep(NA,4)
  
  for ( A in 1:2 ) {
    rw <- 2*(A-1)+1
    #rA.sAB1 <- sum( x[rw,  rw:(rw+1)] )
    #rA.sAB2 <- sum( x[rw+1,rw:(rw+1)] )
    rA.sAB1 <- sum( x[rw,  1:2] )
    rA.sAB2 <- sum( x[rw+1,1:2] )
    nAB1 <- sum(x[rw,])
    nAB2 <- sum(x[rw+1,])
    
    p.s <- (rA.sAB1 + rA.sAB2)/(nAB1 + nAB2)
    statistic[A] <- ((rA.sAB1/nAB1 - rA.sAB2/nAB2)/
                       sqrt(p.s*(1-p.s)*(1/nAB1+1/nAB2)) )
  }
  
  for ( B in 1:2 ) {
    rw <- 2*(A-1)+1
    rB.sA1B <- sum( x[B,c(1,3)] )
    rB.sA2B <- sum(x[B+2,c(1,3)] )
    nA1B <- sum(x[B,])
    nA2B <- sum(x[B+2,])
    
    p.s <- (rB.sA1B + rB.sA2B)/(nA1B + nA2B)
    statistic[B+2] <- ((rB.sA2B/nA2B - rB.sA1B/nA1B)/
                         sqrt(p.s*(1-p.s)*(1/nA1B+1/nA2B)) )
  }
  return(data.frame(stimulus=stimulus, statistic=statistic, 
                    p.value= 2*(pmin(1-pnorm(statistic),pnorm(statistic))) ))
}


#' @export
two_by_two_fit.grt <- function(freq, PS_x = FALSE, PS_y = FALSE, PI = 'none') {
  if(!checkConfusionMatrix(freq)) return(FALSE); # Make sure confusion matrix valid
  delta <- 1/10000; # Tolerance
  freq = freq + 1;   # protection against zeros, which cause the algorithm to explode
  w = 1;  # initialize weight for adjusting step size/preventing oscillation

  # initialize predicted probability matrix
  prob <- matrix(data=0, nrow = 4, ncol = 4)

  # use observed relative frequencies as first estimates of probs
  for (i in 1:4) {
    prob[i,] = freq[i,]/sum(freq[i,]);
  }

  # Get good initial estimates
  initial = initial_point(prob, PS_x, PS_y, PI);
  xpar=initial$xpar; ypar=initial$ypar; rpar=initial$rpar; 
  ps_old=initial$ps_old; rows=initial$rows; npar = length(xpar)+length(ypar)+length(rpar);
  d = matrix(data = 0, nrow = npar, ncol = 1); # Store gradient at estimate
  v <- array(0, dim = c(4,4,3)); # Store estimate variances
  E = matrix(data = 0, nrow = npar, ncol = npar); # Store information matrix
  
  # calculate prob estimates and co-variances for first step
  temp <- estimate_prob_and_var(xpar,ypar,rpar,ps_old);
  prob = temp$prob;
  v = temp$v;
  
  # initial computation of log-likelihood gradient d & information matrix E
  for (i in 1:npar) {
    g <- if(sum(i==xpar)) 1 else (if(sum(i==ypar)) 2 else 3);
    ir <- rows[i,];
    d[i] = sum(sum((freq[ir,]/prob[ir,])
                   *v[ir,, g]));   
    for (j in 1:npar) {
      h <- if(sum(j==xpar)) 1 else (if(sum(j==ypar)) 2 else 3);
      jr <- rows[j,];
      ijr <- ir == 1 & jr == 1;
      if (sum(ijr) > 0) {
        sum_f <- if(is.null(dim(ijr))) sum else rowSums; # rowSums doesn't reduce to sum when d=1...
        K = sum_f(freq[ijr,]) / sum_f(prob[ijr,]);
        L = (sum_f(v[ijr,,g] * v[ijr,,h] / prob[ijr,]) -
             sum_f(v[ijr,,g])* sum_f(v[ijr,,h]) / sum_f(prob[ijr,]));
        E[i,j] = -t(K)*L; 
      }
    }
  }
  Ei = solve(E, diag(npar));  
  ps_new = ps_old - Ei %*% d;
  # iterate!
  
  it = 1;
  df = abs( ps_new - ps_old ) / abs(ps_new);
  dfp_new = t(df) %*% df;
  dfp_old = dfp_new;
  while (dfp_new > delta) {    
    ps_old = ps_new;
    
    # calculate prob estimates and co-variances
    temp <- estimate_prob_and_var(xpar,ypar,rpar,ps_old);
    prob = temp$prob;
    v = temp$v;
    
    # log-likelihood gradient, information matrix
    for (i in 1:npar) {
      g <- if(sum(i==xpar)) 1 else (if(sum(i==ypar)) 2 else 3);
      ir <- rows[i,];
      d[i] = sum(sum((freq[ir,]/prob[ir,])
                     *v[ir,, g]));   
      for (j in 1:npar) {
        h <- if(sum(j==xpar)) 1 else (if(sum(j==ypar)) 2 else 3);
        jr <- rows[j,];
        ijr <- ir == 1 & jr == 1;
        if (sum(ijr) > 0) {
          sum_f <- if(is.null(dim(ijr))) sum else rowSums; # rowSums doesn't reduce to sum when d=1...
          K = sum_f(freq[ijr,]) / sum_f(prob[ijr,]);
          L = (sum_f(v[ijr,,g] * v[ijr,,h] / prob[ijr,]) -
                 sum_f(v[ijr,,g])* sum_f(v[ijr,,h]) / sum_f(prob[ijr,]));
          E[i,j] = -t(K)*L; 
        }
      }
    }
    Ei = solve(E, diag(npar));
    
    if (dfp_new > dfp_old) { #w halving procedure
      w = .5*w;
    }
    
    ps_new = ps_old - w * Ei %*% d;
    df = abs(ps_new-ps_old) / abs(ps_new);
    dfp_old = dfp_new;
    dfp_new = t(df) %*% df;
    it = it + 1;
  }
  
  parameters <- make_parameter_mat(xpar, ypar, rpar, ps_new);
  temp <- estimate_prob_and_var(xpar,ypar,rpar,ps_new);
  prob = temp$prob; 

  # calculate various fit statistics and put them in output structure
  loglike <- sum(sum(freq * log(prob)));
  info_mat <- -solve(E, diag(npar));  
  nll = -loglike;
  aic = 2*npar - 2*loglike;
  bic = npar*log(sum(freq)) - 2*loglike;
  icomp = -loglike + (npar/2)*log(tr(info_mat)/npar)- .5*log(det(info_mat));
  fit <- list(obs=freq2xtabs(freq),fitted=freq2xtabs(prob), estimate=ps_new,
            expd2=E, map=create_n_by_n_mod(PS_x, PS_y, PI, from_2x2 = TRUE), iter=it, 
            loglik=nll);#, aic = aic, bic = bic, icomp = icomp)
  return(grt(parameters, fit, 0, 0))  
}

estimate_prob_and_var <- function(xpar,ypar,rpar,ps_old){
  prob <- matrix(data=0, nrow = 4, ncol = 4)
  v <- array(0, dim = c(4,4,3)); 
  for (i in 1:4){
    # Bookkeeping
    x_i <- if(length(xpar)==2) ceiling(i/2) else i;
    alpha = ps_old[xpar[x_i]];
    y_i <- if(length(ypar)==2) ((i-1) %% 2) + 1 else i;
    kappa = ps_old[ypar[y_i]];
    if (is.null(rpar)) {
      rho = 0;
    } else if (length(rpar) == 1) {
      rho = ps_old[rpar];
    } else {
      rho = ps_old[rpar[i]];
    }
    prob[i,] = prcalc(c(alpha, kappa), matrix(data = c(1, rho, rho, 1), nrow = 2, ncol = 2));
    v[i,,] = vcalc(alpha,kappa,rho);
  }
  return(list(prob=prob,v=v));
}

make_parameter_mat <- function(xpar, ypar, rpar, ps_new){
  if (length(xpar) == 2) {
    mu = c(ps_new[1], ps_new[1], ps_new[2], ps_new[2]); offset = 2;
  } else {
    mu = c(ps_new[1],ps_new[2],ps_new[3],ps_new[4]); offset = 4;
  }
  if (length(ypar) == 2) {
    nu = c(ps_new[offset+1], ps_new[offset+2], ps_new[offset+1], ps_new[offset+2]);
    offset = offset + 2;
  } else {
    nu = c(ps_new[offset+1],ps_new[offset+2], ps_new[offset+3], ps_new[offset+4]);
    offset = offset + 4;
  }
  if (is.null(rpar)) {
    rho = rep(0,4);
  } else if (length(rpar) == 1) {
    rho = rep(ps_new[offset+1],4);
  } else {
    rho = c(ps_new[offset+1], ps_new[offset+2], ps_new[offset+3], ps_new[offset+4]);
  }
  sigma = rep(1.0,4);
  tau = rep(1.0,4);
  return(cbind(mu,sigma,nu,tau,rho));
}
# initialize various scalars and arrays 

# parameters in order:
# mu_x_** mu_y_** rho_**
# where ** = aa, ab, ba, bb
initial_point <- function(prob, PS_x, PS_y, PI) {
  nx =0; ny = 0; nr = 0;xpar = NULL;ypar=NULL;rpar=NULL;
  # Figure out how many params we need
  if (PS_x) {
    xpar=1:2; nx=2; 
    rows=matrix(data=rbind(c(1,1,0,0),c(0,0,1,1)),nrow=2,ncol = 4);
  } else {
    xpar=1:4; nx=4;
    rows=matrix(data=diag(4),nrow=4,ncol=4);
  }
  if (PS_y) {
    ypar = nx + 1:2; ny=2;
    rows=rbind(rows,rbind(c(1,0,1,0),c(0,1,0,1)));
  } else {
    ypar = nx + 1:4; ny=4;
    rows=rbind(rows,diag(4));
  }
  if (PI == 'same_rho') {
    rpar = nx+ny+1; nr = 1;
    rows=rbind(rows,c(1,1,1,1));
  } 
  if (PI == 'none') {
    rpar = nx+ny+1:4; nr=4;
    rows=rbind(rows, diag(4));
  }  
  npar = nx + ny + nr;
  rows = matrix(data=as.logical(rows),ncol=4,nrow=npar);
  param_estimate = matrix(data = 0, nrow= npar, ncol = 1); # For param estimates
  # initial estimates: y means
  if (PS_x) {
    for (i in c(1,3)) { 
      param_estimate[xpar[ceiling(i/2)]] = -qnorm(.5*(prob[i,1]   + prob[i,2]) 
                                                + .5*(prob[i+1,1] + prob[i+1,2]));}
  } else {
    for (i in 1:4) {
      param_estimate[xpar[i]] = -qnorm(prob[i,1] + prob[i,2]);}
  }
  # initial estimates: y means
  if (PS_y) {
    for (i in c(1,2)) {
      param_estimate[ypar[i]] = -qnorm(.5*(prob[i,1]   + prob[i,3])
                                     + .5*(prob[i+2,1] + prob[i+2,3])); }
  } else {
    for (i in 1:4) {
      param_estimate[ypar[i]] = -qnorm(prob[i,1] + prob[i,3]);}
  }  
  # initial estimates: correlation  
  if (PI=='same_rho') {
     r = cos(pi/(1+sqrt((.25*sum(prob[,4])*.25*sum(prob[,1]))
                       /(.25*sum(prob[,3])*.25*sum(prob[,2])))));
    if (r <= -1) { r = -.95; }
    else if (r >= 1) { r = .95; }
    param_estimate[rpar] = r;
  } else if (PI == 'none') {
    for (i in 1:4) {
      r = cos(pi/(1+sqrt((prob[i,4]*prob[i,1])/(prob[i,3]*prob[i,2]))));
      if (r <= -1) { r = -.95; }
      else if (r >= 1) { r = .95; }
      param_estimate[rpar[i]] = r;
    }
  }
  return(list(xpar=xpar, ypar=ypar, rpar=rpar, ps_old=param_estimate, rows=rows))
}  
  
create_two_by_two_mod <- function(PS_x, PS_y ,PI) {
  mod <- matrix(data = 0, nrow = 1, ncol = 7);
  if (PS_x) mod[1] = 1;
  if (PS_y) mod[2] = 1;
  if (PI == 'all') {
    mod[3:6] = rep(1,times=4);
  } else { 
    if (PI == 'same_rho') mod[7] = 1;
  }
  return(mod);
}

# In 2x2 case, it's typical to use a 4x4 frequency matrix w/ each row being a stim
# and each col being the freqency of responding "aa", "ab", "ba", "bb", respectively, 
# to that stim. Wickens' code for nxn case requires data in xtabs format.
freq2xtabs <- function(freq) {
  xdim = dim(freq)[1]; 
  ydim = dim(freq)[2];
  d = as.data.frame(matrix(rep(x=0,times=xdim*ydim*4), nrow = xdim*ydim, ncol = 4));
  stimuli = c('aa','ab','ba','bb')
  names(d) <- c("Stim", "L1", "L2", "x");
  for (i in 1:4) {
    for (j in 1:4) {
      d[4*(i-1) + j,] = c(stimuli[i], floor((j+1) / 2), ((j-1) %% 2) + 1, freq[i,j]);
    }
  }
  d$Stim <- ordered(d$Stim,levels=stimuli);
  d$L1  <- ordered(d$L1);
  d$L2 <- ordered(d$L2);
  d$x <- as.numeric(d$x);
  return(xtabs(x~L1+L2+Stim, d));
}

# calculate predicted response probabilities for the given stimulus
# b/c we assume decisional sep, each response is a quadrant of Cartesian plane
prcalc <- function(mean, cov) {
  pr <- matrix(data = 0, nrow = 1, ncol = 4);
  pr[1,1] = sadmvn(lower = c(-Inf, -Inf), upper = c(0, 0), mean, cov);
  pr[1,2] = sadmvn(lower = c(-Inf, 0), upper = c(0, +Inf), mean, cov);
  pr[1,3] = sadmvn(lower = c(0, -Inf), upper = c(+Inf, 0), mean, cov);
  pr[1,4] = sadmvn(lower = c(0, 0), upper = c(+Inf, +Inf), mean, cov);
  return(pr);
}

# calculate v-matrix elements
vcalc <- function(ap,kp,rh) {        
  ve <- matrix(data = 0, ncol = 3, nrow = 4)
  d_mx <- matrix(data = 0, ncol = 2, nrow = 2)
  d_my <- matrix(data = 0, ncol = 2, nrow = 2)

  d_mx_arg = (rh*ap-kp)/sqrt(1-rh^2);
  d_mx[1,1] = -dnorm(-ap)*pnorm( d_mx_arg );
  d_mx[1,2] = -dnorm(-ap);
  d_mx[2,1] = 0;
  d_mx[2,2] = 0;
  
  ve[1,1] =  d_mx[1,1];
  ve[2,1] =  d_mx[1,2] - d_mx[1,1];
  ve[3,1] =  d_mx[2,1] - d_mx[1,1];
  ve[4,1] =  d_mx[2,2] - d_mx[2,1] - d_mx[1,2] + d_mx[1,1];
  
  d_my_arg = (rh*kp-ap)/sqrt(1-rh^2);
  d_my[1,1] = -dnorm(-kp)*pnorm( d_my_arg );
  d_my[1,2] = 0;
  d_my[2,1] = -dnorm(-kp);
  d_my[2,2] = 0;
  
  ve[1,2] = d_my[1,1];
  ve[2,2] = d_my[1,2] - d_my[1,1];
  ve[3,2] = d_my[2,1] - d_my[1,1];
  ve[4,2] = d_my[2,2] - d_my[1,2] - d_my[2,1] + d_my[1,1];
  
  S_aa = matrix(c(1, rh, rh, 1), ncol = 2, nrow = 2);
  d_rho = dmnorm(c(-ap, -kp),mean = c(0, 0), varcov = S_aa);
  
  ve[1,3] = d_rho;
  ve[2,3] = -d_rho;
  ve[3,3] = -d_rho;
  ve[4,3] = d_rho;

  return(ve);
}

get_fit_params <- function(grt_obj) {
  d = distribution.parameters(bb=grt_obj);
  return(list(aa=d[1,c(1,3,5)], ab = d[2,c(1,3,5)], 
              ba = d[3, c(1,3,5)], bb = d[4,c(1,3,5)]));
}

#' @export
two_by_two_plot.grt <- function(fit_params, xlab1, ylab1) {
  bin_width= .25;
  bivar_obj = bivar_norm(fit_params, bin_width);
  dist = bivar_obj$dist;
  x = bivar_obj$x;
  y = bivar_obj$y;
  xlims = c(x[1], x[length(x)]);
  ylims = c(y[1], y[length(y)]);
  level = c(.1, .1); # Determines which contour to plot
  
  # Plot Gaussian contours 
  plot.new();
  old_mar <- par()$mar;
  par(fig=c(.3,1,.3,1), mar=c(2,.5,1,1));
  for (i in 1:4) {
    cond = fit_params[[i]];
    par(new = TRUE);
    contour(x = x, y = y, z = dist[i,,], 
            levels = level, xlim = xlims, ylim = ylims,
            drawlabels = FALSE, axes = FALSE, asp=1);
    points(cond[1], cond[2], pch = '+');
    title(xlab = 'x', ylab = 'y');
  }
  
  # Plot decision bounds
  abline(h=0);
  abline(v=0);
  
  # compute marginals
  margx = margy = list(aa=NULL,ab=NULL,ba=NULL, bb=NULL);
  for (i in 1:4) {
    cond = fit_params[[i]];
    xvar = bivar_obj$covars[i,1,1];
    yvar = bivar_obj$covars[i,2,2];
    margx[[i]] = (1 / sqrt(2*pi*xvar)) * exp(-.5*((x-cond[1])/sqrt(xvar))^2);
    margy[[i]] = (1 / sqrt(2*pi*yvar)) * exp(-.5*((y-cond[2])/sqrt(yvar))^2);
  }
  
  # Plot X marginals
  par(fig=c(.3,1,.05,.4), mar = c(3.8,.5,1,1), pty='m', xaxt = 'n', yaxt = 'n', new=TRUE);
  plot(margx$aa,type='l', xlab = xlab1, ylab = NULL);
  lines(margx$ab,type='l',lty=2);
  lines(margx$ba,type='l',lty=2);
  lines(margx$bb,type='l');
  
  # Plot Y marginals
  par(fig=c(.05,.3,.3,1), mar = c(2, 3.75, 1, 0), pty='m', xaxt = 'n', yaxt = 'n', new=TRUE);
  plot(margy$aa,bivar_obj$y,type='l', xlab = NULL, ylab = ylab1);
  lines(margy$ab,bivar_obj$y, type='l',lty=2);
  lines(margy$ba,bivar_obj$y,type='l',lty=2);
  lines(margy$bb,bivar_obj$y,type='l');
  
  # Reset graphical par for later
  par(mar = old_mar, fig = c(0,1,0,1));
}

# Returns the density of bivariate normal distribution over a grid of x,y values
bivar_norm <- function(fit_params, bin_width) {
  # Compute covariance matrices
  covars <- array(0, dim = c(4,2,2));
  for (i in 1:4) {
    cond = fit_params[[i]];
    covars[i,,] <- matrix(data=c(1, cond[3], cond[3], 1), ncol = 2, nrow = 2);
  }
  # Set the range that will be graphed
  densfact = 3; # how many SDs?
  xrange <- array(0, dim = c(4,2));
  yrange <- array(0, dim = c(4,2));
  for (i in 1:4) {
    cond = fit_params[[i]];
    x_sd = densfact*sqrt(covars[i,1,1]);
    y_sd = densfact*sqrt(covars[i,2,2]);
    xrange[i,] = c(cond[1] - x_sd, cond[1] + x_sd);
    yrange[i,] = c(cond[2] - y_sd, cond[2] + y_sd);
  }  
  x = seq(from = min(xrange), to = max(xrange), by = bin_width);
  y = seq(from = min(yrange), to = max(yrange), by = bin_width);
  lenx = length(x);
  leny = length(y);
  dist = array(0, dim = c(4, lenx, leny));  
  for (i in 1:4) {
    cond = fit_params[[i]];
    for (y_ind in seq(from=leny, to=1, by = -1)) {
      for (x_ind in seq(from=1, to=lenx)) {
        dist[i, x_ind, y_ind] = dmnorm(c(x[x_ind], y[y_ind]), 
                                       mean = c(cond[1], cond[2]), 
                                       varcov = covars[i,,]);
      }
    }
  }
  return(list(dist = dist, 
              x = x,
              y = y,
              covars = covars));
}

checkConfusionMatrix <- function(x) {
  dimx <- dim(x)[1]
  if( dimx != dim(x)[2]){ 
    cat("Confusion matrix must have an equal number of rows and columns!\n")
    return(FALSE)
  }
  
  if(max(x)<=1 & min(x) >=0) {
    if(all( apply(x, 1, sum) == rep(1,dimx))) {
      return(TRUE)
    } else {
      cat("The rows of confusion probability matrix must sum to one!\n")
      return(FALSE)
    }
  } else {
    return(TRUE)}
}

  