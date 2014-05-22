#' Uses Newton-Raphson gradient descent to fit a Gaussian model 
#'
#' @param freq 4x4 confusion matrix (containing counts). Assumes row/col order: aa, ab, ba, bb.
#' @return data frame containing z-scores and p-values for all four tests
#' @export
#' @examples
gaussian_fit <- function(freq) {
  if(!checkConfusionMatrix(freq)) {
    return(FALSE)
  }

  # Tolerance -- stops when new parameter values are less than delta away from old
  delta <- 1/10000; 

  # In principle, we could add constraints to our model, but 
  # we've only implemented the most general, unconstrained model so far
  mods <- matrix(data = 0, nrow = 1, ncol = 7);
  
  # protection against zeros, which cause the algorithm to explode
  freq = freq + 1;

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
  for (i in 1:4) {
    alpha = ps_old[xpar[i]]; # x mean
    kappa = ps_old[ypar[i]]; # y mean
    rho = ps_old[rpar[i]];
    prob[i,] = prcalc(c(alpha, kappa), matrix(data = c(1, rho, rho, 1), nrow = 2, ncol = 2));
    v[i,,] = vcalc(alpha,kappa,rho);
  }  
    
  # log-likelihood gradient, information matrix
  for (i in 1:npar) {
    cond = ((i-1) %% 4) + 1;
    d[i] = sum(sum((freq[cond,]/prob[cond,])
                   *v[cond,, ceiling(i/4)]));    
    K = sum( freq[cond,]);#./sum( prob(cond,:) , 2 );
    L = sum( v[cond,,ceiling(i/4)]   * v[cond,,ceiling(i/4)] / prob[cond,] , 2 ) -
        sum( v[cond,,ceiling(i/4)], 2) * sum(v[cond,, ceiling(i/4)],2);#./sum(prob(cond,),2);
    E[i,i] = -K*L; 
  }
  
  print(E)
  Ei = solve(E, diag(npar));
  
  ps_new = ps_old - Ei %*% d;
  
  # iterate!
  
  it = 1;
  
  df = abs( ps_new - ps_old ) / abs(ps_new);
  
  dfp_new = t(df) %*% df;
  
  dfp_old = dfp_new;
  
  print(c(it))
  
  while (dfp_new > delta) {
    
    ps_old = ps_new;
    
    # stimulus aa
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
      cond = ((i - 1) %% 4) + 1;
      d[i] = sum(sum((freq[cond,]/prob[cond,])
                     *v[cond,, ceiling(i/4)]));    
      K = sum( freq[cond,]);#./sum( prob(cond,:) , 2 );
      L = sum( v[cond,,ceiling(i/4)]   * v[cond,,ceiling(i/4)] / prob[cond,] , 2 ) -
        sum( v[cond,,ceiling(i/4)], 2) * sum(v[cond,, ceiling(i/4)],2);#./sum(prob(cond,),2);
      E[i,i] = -K*L; 
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
  
  print(ps_new)
  # put fitted parameters into the output structure
  
#   outfit.p(1,m) = ps_new(1);
#   outfit.p(2,m) = ps_new(2);
#   outfit.p(3,m) = ps_new(3);
#   outfit.p(4,m) = ps_new(4);
#   
#   outfit.p(5,m) = ps_new(nx+1);
#   outfit.p(6,m) = ps_new(nx+2);
#   outfit.p(7,m) = ps_new(nx+3);
#   outfit.p(8,m) = ps_new(nx+4);
#   
#   outfit.p(9,m) = ps_new(nx+ny+1);
#   
#   offset = 0+mods(m,3)==0;
#   outfit.p(10,m) = ps_new(nx+ny+offset+1);
#   
#   offset = sum(mods(m,3:4)==0);
#   outfit.p(11,m) = ps_new(nx+ny+offset+1);
#   
#   offset = sum(mods(m,3:5)==0);
#   outfit.p(12,m) = ps_new(nx+ny+offset+1);
#   
#   # put information matrix in the output structure
#   outfit.M(:,:,m) = zeros(12);
#   infoM = (-E)\eye(npar);
#   outfit.M(1:npar,1:npar,m) = infoM;
#   
#   # for calculating log-likelihood, AIC, BIC
#   
#   alpha = ps_new(xpar(1)); # x mean, stim aa
#   kappa = ps_new(ypar(1)); # y mean, stim aa
#   rho = ps_new(rpar(1));
#   
#   prob(1,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
#   
#   alpha = ps_new(xpar(2));
#   kappa = ps_new(ypar(2)); # y mean, stim ab
#   offset = 0+[mods(m,3)==0];
#   rho = ps_new(rpar(offset+1));
#   
#   prob(2,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
#   
#   alpha = ps_new(xpar(3)); # x mean, stim ba
#   kappa = ps_new(ypar(3)); # y mean, stim ba
#   offset = sum(mods(m,3:4)==0);
#   rho = ps_new(rpar(offset+1));
#   
#   prob(3,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
#   
#   alpha = ps_new(xpar(4)); # x mean, stim bb
#   kappa = ps_new(ypar(4)); # y mean, stim bb
#   offset = sum(mods(m,3:5)==0);
#   rho = ps_new(rpar(offset+1));
#   
#   prob(4,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
#   
#   # put fitted probabilities in output structure
#   outfit.pr(:,:,m) = prob;
#   
#   # calculate fit statistics and put them in output structure
#   loglike = sum(sum(freq.*log(prob)));
#   outfit.nll(m) = -loglike;
#   outfit.aic(m) = 2*npar - 2*loglike;
#   outfit.bic(m) = npar*log(sum(sum(freq))) - 2*loglike;
#   outfit.icomp(m) = -loglike + (npar/2)*log(trace(infoM)/npar)-.5*log(det(infoM));
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
  