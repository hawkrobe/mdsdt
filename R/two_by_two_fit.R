require()

#' Fit a bivariate Gaussian model for each response distribution using Newton-Raphson gradient descent
#'
#' @param freq 4x4 confusion matrix containing counts. Assumes row/col order: aa, ab, ba, bb.
#' @return List containing the following elements:
#'  \item{parameters}{The fitted means and correlations for all 4 response distributions}
#'  \item{prob}{The predicted response probabilities according to the fit}
#'  \item{info_mat}{The information matrix (inverse of Fisher Information matrix) for all parameters}
#'  \item{nll}{The negative log likelihood}
#'  \item{aic}{Akaike information criterion, for model comparison}
#'  \item{bic}{Bayesian information criterion, for model comparison}
#'  \item{icomp}{Bozdogan's information complexity (ICOMP), for model comparison}
#' @examples gaussian_fit(observerB)
#' @export
two_by_twofit.grt <- function(freq, PS_x = FALSE, PS_y = FALSE, PI = 'none') {
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
  prob=temp$prob; 

  # calculate various fit statistics and put them in output structure
  loglike <- sum(sum(freq * log(prob)));
  info_mat <- -solve(E, diag(npar));  
  nll = -loglike;
  aic = 2*npar - 2*loglike;
  bic = npar*log(sum(freq)) - 2*loglike;
  icomp = -loglike + (npar/2)*log(tr(info_mat)/npar)- .5*log(det(info_mat));
  fit <- list(obs=freq,fitted=prob, estimate=ps_new,
            expd2=E, map=create_n_by_n_mod(PS_x, PS_y, PI), iter=it, 
            nll=-loglike, aic = aic, bic = bic, icomp = icomp)
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
    rho = rep(1,4);
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
  print(rows);
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

#' Checks whether the input is a valid confusion matrix
#'
#' @param x: four-by-four confusion matrix of counts of probabilities
#' @return Boolean indicating whether input is valid
#' @export
#' @examples
#' checkConfusionMatrix(matrix(c(3,1,1,4,3,3),3)) 
#' checkConfusionMatrix(matrix(c(.5,.2,.2,.2,
#'                               .2,.5,.1,.2,
#'                               .1,.2,.5,.1,
#'                               .2,.2,.3,.5), 4))
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
  } else if(min(x) >= 0) {
    if(all( apply(x, 1, sum) == sum(x[1,]))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {return(FALSE)}
}

  