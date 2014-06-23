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
two_by_twofit.grt <- function(freq, mod) {
  if(!checkConfusionMatrix(freq)) return(FALSE); # Make sure confusion matrix valid
  delta <- 1/10000; # Tolerance
  freq = freq + 1;   # protection against zeros, which cause the algorithm to explode

  # initialize weight for adjusting step size/preventing oscillation
  w = 1;

  # This would index which model we're fitting if there was more than one row in mods
  m = 1;

  # initialize predicted probability matrix
  nrow = dim(freq)[1];
  ncol = dim(freq)[2]; 
  prob <- matrix(data=0, nrow = nrow, ncol = ncol)

  # use observed relative frequencies as first estimates of probs
  for (i in 1:nrow) {
    prob[i,] = freq[i,]/sum(freq[i,]);
  }

  # initialize various scalars and arrays (assume no constraints)
  
  # parameters in order:
  # mu_x_** mu_y_** rho_**
  # where ** = aa, ab, ba, bb
  rows = matrix(data=0, nrow = 12, ncol = 4);
  nx = 4;
  ny = 4;
  nr = 4;
  xpar = 1:4;
  ypar = nx + 1:4;
  rpar = nx + ny + 1:4;
  rows[xpar,] = diag(4);
  rows[ypar,] = diag(4);
  rows[rpar,] = diag(4);
  rows = as.logical(rows);
  npar = nx + ny + nr;
  d = matrix(data = 0, nrow = npar, ncol = 1);
  ps_old = matrix(data = 0, nrow= npar, ncol = 1);
  v <- array(0, dim = c(4,4,3))
  E = matrix(data = 0, nrow = npar, ncol = npar);
  
  for (i in 1:4) {
    # initial estimates: means
    ps_old[xpar[i]] = -qnorm(prob[i,1] + prob[i,2]);
    ps_old[ypar[i]] = -qnorm(prob[i,1] + prob[i,3]);
    # initial estimates: correlation
    r = cos(pi/(1+sqrt((prob[i,4]*prob[i,1])/(prob[i,3]*prob[i,2]))));
    if (r <= -1) {
      r = -.95;
    } else if (r >= 1) {
      r = .95;
    }
    ps_old[rpar[i]] = r;
  }
  
  # calculate log-like gradient and information matrix for first step
  for (i in 1:4){
    alpha = ps_old[xpar[i]]; # x mean
    kappa = ps_old[ypar[i]]; # y mean
    rho = ps_old[rpar[i]];
    prob[i,] = prcalc(c(alpha, kappa), matrix(data = c(1, rho, rho, 1), nrow = 2, ncol = 2));
    v[i,,] = vcalc(alpha,kappa,rho);
  }  
    
  # log-likelihood gradient, information matrix
  for (i in 1:npar) {
    condi = ((i-1) %% 4) + 1;
    d[i] = sum(sum((freq[condi,]/prob[condi,])
                   *v[condi,, ceiling(i/4)]));   
    for (j in 1:npar) {
      condj = ((j-1) %% 4) + 1;
      if (condi == condj) {
        K = sum( freq[condi,]);#/sum(prob[condi,]);
        L = sum( v[condi,,ceiling(i/4)]   * v[condi,,ceiling(j/4)] / prob[condi,]) -
            sum( v[condi,,ceiling(i/4)]) * sum(v[condi,, ceiling(j/4)]);#./sum(prob(condi,),2);
        E[i,j] = -K*L; 
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
    # calculate log-like gradient and information matrix for first step
    for (i in 1:4) {
      alpha = ps_old[xpar[i]]; # x mean
      kappa = ps_old[ypar[i]]; # y mean
      rho = ps_old[rpar[i]];
      prob[i,] = prcalc(c(alpha, kappa), matrix(data = c(1, rho, rho, 1), nrow = 2, ncol = 2));
      v[i,,] = vcalc(alpha,kappa,rho);
    }  
   
    # log-likelihood gradient, information matrix
    for (i in 1:npar) {
      condi = ((i-1) %% 4) + 1;
      d[i] = sum(sum((freq[condi,]/prob[condi,])
                     *v[condi,, ceiling(i/4)]));   
      for (j in 1:npar) {
        condj = ((j-1) %% 4) + 1;
        if (condi == condj) {
          K = sum( freq[condi,]);#/sum(prob[condi,]);
          L = sum( v[condi,,ceiling(i/4)]   * v[condi,,ceiling(j/4)] / prob[condi,]) -
            sum( v[condi,,ceiling(i/4)]) * sum(v[condi,, ceiling(j/4)]);#./sum(prob(condi,),2);
          E[i,j] = -K*L; 
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
    print(it)
    print(dfp_new)
  }
  
  # calculate fit statistics and put them in output structure
  parameters <- data.frame(aa = c(ps_new[1],ps_new[5],ps_new[9]),
                           ab = c(ps_new[2],ps_new[6],ps_new[10]),
                           ba = c(ps_new[3],ps_new[7],ps_new[11]),
                           bb = c(ps_new[4],ps_new[8],ps_new[12]));
  
  for (i in 1:4) {
    alpha = ps_new[xpar[i]]; # x mean
    kappa = ps_new[ypar[i]]; # y mean
    rho = ps_new[rpar[i]];
    prob[i,] = prcalc(c(alpha, kappa), matrix(data = c(1, rho, rho, 1), nrow = 2, ncol = 2));
  }  
  
  info_mat <- -solve(E, diag(npar));
  loglike = sum(sum(freq * log(prob)));
  nll = -loglike;
  aic = 2*npar - 2*loglike;
  bic = npar*log(sum(freq)) - 2*loglike;
  icomp = -loglike + (npar/2)*log(tr(info_mat)/npar)- .5*log(det(info_mat));
  return(list(parameters = parameters, prob = prob, info_mat = info_mat, 
              nll = nll, aic = aic, bic = bic, icomp = icomp))
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
  