#' Uses Newton-Raphson gradient descent to fit a Gaussian model 
#'
#' @param freq 4x4 confusion matrix (containing counts). Assumes row/col order: aa, ab, ba, bb.
#' @return data frame containing z-scores and p-values for all four tests
#' @export
#' @examples
gaussian_fit <- function(freq) {
  if(!checkConfusionMatrix(x)) {
    return(FALSE)
  }

  # Step size. Adjust as needed
  delta <- 1/10000; 

  mods <- matrix(data = 0, nrow = 1, ncol = 7);
  
  # protection against zeros, which cause the algorithm to explode
  freq = freq + 1;

  # parameters in order:
  # mu_x_** mu_y_** rho_**
  # where ** = aa, ab, ba, bb

  # initialization

  # notes:
  # should initialize correlation parameters in more informative way
  # need to figure right out criterion for continued iterations

  # initialize weight for adjusting step size/preventing oscillation
  w = 1;

  # Part of our effort to reduce generality of model
  m = 1;

  # initialize predicted probability matrix
  nrow = dim(freq)[1];
  ncol = dim(freq)[2]; 
  prob <- matrix(data=0, nrow = nrow, ncol = ncol)

  # use observed relative frequencies as first estimates of probs
  for (i in 1:nrow) {
    prob[i,] = freq[i,:]/sum(freq[i,:]);
  }

  # initialize various scalars and arrays
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
  E = matrix(data = 0, nrow = npar, ncol = 1);
  xc_aa = 0;
  xc_ab = 0;
  xc_ba = 0;
  xc_bb = 0;
  yc_aa = 0;
  yc_ab = 0;
  yc_ba = 0;
  yc_bb = 0;
  
  # count the # of parameters and specify matrix
  # to pick out relevant cells in data matrix for sums used later on
    
  # stim aa means
  
  xc_aa = qnorm(prob[1,1] + prob[1,2]);
  yc_aa = qnorm(prob[1,1] + prob[1,3]);
  
  ps_old[xpar[1]] = -xc_aa;
  ps_old[ypar[1]] = -yc_aa;
  
  
  # stim ab means
  
  xc_ab = qnorm(prob[2,1] + prob[2,2]);
  yc_ab = qnorm(prob[2,1] + prob[2,3]);
  
  ps_old[xpar[2]] = -xc_ab;
  
  ps_old[ypar(2)] = -yc_ab;
  
  
  # stim ba means
  
  
  xc_ba = qnorm(prob(3,1) + prob(3,2));
  ps_old(xpar(3)) = -xc_ba;
  
  yc_ba = qnorm(prob(3,1) + prob(3,3));
  ps_old(ypar(3)) = -yc_ba;
  
  # stim bb means
  
  xc_bb = qnorm(prob(4,1) + prob(4,2));
  ps_old(xpar(4)) = -xc_bb;
  
  yc_bb = qnorm(prob(4,1) + prob(4,3));
  ps_old(ypar(4)) = -yc_bb;
  
  # correlation parameters
  
  # separate correlation parameters in separate distributions
  
  # aa
  a = prob(1,4);
  b = prob(1,3);
  c = prob(1,2);
  e = prob(1,1);
  r = cos(pi/(1+sqrt((a*e)/(b*c))));
  if (r <= -1) {
    r = -.95;
  } else if (r >= 1) {
    r = .95;
  }
  rho_aa = r;
  
  # ab
  a = prob(2,4);
  b = prob(2,3);
  c = prob(2,2);
  e = prob(2,1);
  r = cos(pi/(1+sqrt((a*e)/(b*c))));
  if (r <= -1) {
    r = -.95;
  } else if (r >= 1) {
    r = .95;
  }
  rho_ab = r;
  
  # ba
  a = prob(3,4);
  b = prob(3,3);
  c = prob(3,2);
  e = prob(3,1);
  r = cos(pi/(1+sqrt((a*e)/(b*c))));
  if (r <= -1) {
    r = -.95;
  } else if (r >= 1) {
    r = .95;
  }
  rho_ba = r;
  
  # bb
  
  a = prob(4,4);
  b = prob(4,3);
  c = prob(4,2);
  e = prob(4,1);
  r = cos(pi/(1+sqrt((a*e)/(b*c))));
  if (r <= -1) {
    r = -.95;
  } else if (r >= 1) {
    r = .95;
  }
  rho_bb = r;
  
  ps_old(rpar) = [rho_aa; rho_ab; rho_ba; rho_bb];
  
  # calculate log-like gradient and information matrix for first step
  
  # stimulus aa
  
  alpha = ps_old(xpar(1)); # x mean, stim aa
  kappa = ps_old(ypar(1)); # y mean, stim aa
  rho = ps_old(rpar(1));
  
  prob(1,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
  
  v(1,:,:) = vcalc(alpha,kappa,rho);
  
  # stimulus ab
  
  alpha = ps_old(xpar(2));
  kappa = ps_old(ypar(2)); # y mean, stim ab
  
  offset = 0+[mods(m,3)==0];
  rho = ps_old(rpar(offset+1));
  
  
  prob(2,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
  
  v(2,:,:) = vcalc(alpha,kappa,rho);
  
  # stimulus ba
  
  alpha = ps_old(xpar(3)); # x mean, stim ba
  kappa = ps_old(ypar(3)); # y mean, stim ba
  offset = sum(mods(m,3:4)==0);
  rho = ps_old(rpar(offset+1));
  
  prob(3,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
  
  v(3,:,:) = vcalc(alpha,kappa,rho);
  
  # stimulus bb
  
  alpha = ps_old(xpar(4)); # x mean, stim bb
  kappa = ps_old(ypar(4)); # y mean, stim bb
  offset = sum(mods(m,3:5)==0);
  rho = ps_old(rpar(offset+1));
  
  prob(4,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
  
  v(4,:,:) = vcalc(alpha,kappa,rho);
  
  # log-likelihood gradient, information matrix
  for (i in 1:npar) {
    if (sum(i==xpar)) {
      g = 1;
    } else if (sum(i==ypar)) {
      g = 2;
    } else if (sum(i==rpar)) {
      g = 3;
    }
    ir = rows(i,:);
    d(i) = sum(sum( (freq(ir,:)./prob(ir,:)).*v(ir,:,g) ));
    for (j in 1:npar) {
      if (sum(j==xpar)) {
        h = 1;
      } else if (sum(j==ypar)) {
        h = 2;
      } else if (sum(j==rpar)) {
        h = 3;
      }
      jr = rows(j,:);
      
      ijr = ir==1 & jr==1;
      
      if (sum(ijr) > 0) {
        K = sum( freq(ijr,:) , 2 )./sum( prob(ijr,:) , 2 );
        L = sum( v(ijr,:,g).*v(ijr,:,h)./prob(ijr,:) , 2 ) - ...
        sum(v(ijr,:,g),2).*sum(v(ijr,:,h),2)./sum(prob(ijr,:),2);
        E(i,j) = -K'*L; 
      }
      } 
      }
        
        Ei = E\eye(npar);
        
        ps_new = ps_old - Ei*d;
        
        # iterate!
        
        it = 1;
        
        df = abs( ps_new-ps_old )./abs(ps_new);
        
        dfp_new = df'*df;
        
        dfp_old = dfp_new;
        
        disp([m it])
        
        while (dfp_new > delta) {
          
          ps_old = ps_new;
          
          # stimulus aa
          
          alpha = ps_old(xpar(1)); # x mean, stim aa
          kappa = ps_old(ypar(1)); # y mean, stim aa
          rho = ps_old(rpar(1));
          
          prob(1,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
          
          v(1,:,:) = vcalc(alpha,kappa,rho);
          
          # stimulus ab
          
          alpha = ps_old(xpar(2));
          kappa = ps_old(ypar(2)); # y mean, stim ab
          offset = 0+[mods(m,3)==0];
          rho = ps_old(rpar(offset+1));
          
          prob(2,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
          
          v(2,:,:) = vcalc(alpha,kappa,rho);
          
          # stimulus ba
          
          alpha = ps_old(xpar(3)); # x mean, stim ba
          kappa = ps_old(ypar(3)); # y mean, stim ba
          offset = sum(mods(m,3:4)==0);
          rho = ps_old(rpar(offset+1));
          
          prob(3,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
          
          v(3,:,:) = vcalc(alpha,kappa,rho);
          
          # stimulus bb
          
          alpha = ps_old(xpar(4)); # x mean, stim bb
          kappa = ps_old(ypar(4)); # y mean, stim bb
          offset = sum(mods(m,3:5)==0);
          rho = ps_old(rpar(offset+1));
          
          prob(4,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
          
          v(4,:,:) = vcalc(alpha,kappa,rho);
          
          # log-likelihood gradient, information matrix
          for (i in 1:npar) {
            if (sum(i==xpar)) {
              g = 1;
            } else if (sum(i==ypar)) {
              g = 2;
            } else if (sum(i==rpar)) {
              g = 3;
            }
            ir = rows(i,:);
            d(i) = sum(sum( v(ir,:,g).*freq(rows(i,:),:)./prob(rows(i,:),:) ));
            for (j in 1:npar) {
              if (sum(j==xpar)) {
                h = 1;
              } else if (sum(j==ypar)) {
                h = 2;
              } else if (sum(j==rpar)) {
                h = 3;
              }
              jr = rows(j,:);
              if (ir+0)*(jr+0)' > 0) {
              ijr = ir==1 & jr==1;
              K = -1*(sum( freq(ijr,:) , 2 )./sum( prob(ijr,:) , 2 ));
              L = sum( v(ijr,:,g).*v(ijr,:,h)./prob(ijr,:) , 2 ) - ...
              sum(v(ijr,:,g),2).*sum(v(ijr,:,h),2)./sum(prob(ijr,:),2);
              E(i,j) = K'*L;
            }
            }
            }
        
        
        Ei = E\eye(npar);
        
        if (dfp_new > dfp_old) { #w halving procedure
                                w = .5*w;
        }
        
        ps_new = ps_old - w*Ei*d;
        df = abs(ps_new-ps_old)./abs(ps_new);
        dfp_old = dfp_new;
        dfp_new = df'*df;
        it = it + 1;
        disp([m it])
            }
        
        # put fitted parameters into the output structure
        
        outfit.p(1,m) = ps_new(1);
        outfit.p(2,m) = ps_new(2);
        outfit.p(3,m) = ps_new(3);
        outfit.p(4,m) = ps_new(4);
        
        outfit.p(5,m) = ps_new(nx+1);
        outfit.p(6,m) = ps_new(nx+2);
        outfit.p(7,m) = ps_new(nx+3);
        outfit.p(8,m) = ps_new(nx+4);
        
        outfit.p(9,m) = ps_new(nx+ny+1);
        
        offset = 0+mods(m,3)==0;
        outfit.p(10,m) = ps_new(nx+ny+offset+1);
        
        offset = sum(mods(m,3:4)==0);
        outfit.p(11,m) = ps_new(nx+ny+offset+1);
        
        offset = sum(mods(m,3:5)==0);
        outfit.p(12,m) = ps_new(nx+ny+offset+1);
        
        # put information matrix in the output structure
        outfit.M(:,:,m) = zeros(12);
        infoM = (-E)\eye(npar);
        outfit.M(1:npar,1:npar,m) = infoM;
        
        # for calculating log-likelihood, AIC, BIC
        
        alpha = ps_new(xpar(1)); # x mean, stim aa
        kappa = ps_new(ypar(1)); # y mean, stim aa
        rho = ps_new(rpar(1));
        
        prob(1,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
        
        alpha = ps_new(xpar(2));
        kappa = ps_new(ypar(2)); # y mean, stim ab
        offset = 0+[mods(m,3)==0];
        rho = ps_new(rpar(offset+1));
        
        prob(2,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
        
        alpha = ps_new(xpar(3)); # x mean, stim ba
        kappa = ps_new(ypar(3)); # y mean, stim ba
        offset = sum(mods(m,3:4)==0);
        rho = ps_new(rpar(offset+1));
        
        prob(3,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
        
        alpha = ps_new(xpar(4)); # x mean, stim bb
        kappa = ps_new(ypar(4)); # y mean, stim bb
        offset = sum(mods(m,3:5)==0);
        rho = ps_new(rpar(offset+1));
        
        prob(4,:) = prcalc([alpha kappa],eye(2) + [0 rho; rho 0]);
        
        # put fitted probabilities in output structure
        outfit.pr(:,:,m) = prob;
        
        # calculate fit statistics and put them in output structure
        loglike = sum(sum(freq.*log(prob)));
        outfit.nll(m) = -loglike;
        outfit.aic(m) = 2*npar - 2*loglike;
        outfit.bic(m) = npar*log(sum(sum(freq))) - 2*loglike;
        outfit.icomp(m) = -loglike + (npar/2)*log(trace(infoM)/npar)-.5*log(det(infoM));
        
        
        # calculate predicted probabilities
        function pr = prcalc(mv,SM);
        
        pr(1,1) = mvncdf([-inf -inf],[0 0],mv,SM);
        pr(1,2) = mvncdf([-inf 0],[0 inf],mv,SM);
        pr(1,3) = mvncdf([0 -inf],[inf 0],mv,SM);
        pr(1,4) = mvncdf([0 0],[inf inf],mv,SM);
        
        # calculate v-matrix elements
        function ve = vcalc(ap,kp,rh)
        
        d_mx_arg = (rh*ap-kp)/sqrt(1-rh^2);
        d_mx(1,1) = -normpdf(-ap)*normcdf( d_mx_arg );
        d_mx(1,2) = -normpdf(-ap);
        d_mx(2,1) = 0;
        d_mx(2,2) = 0;
        
        ve(1,1,1) =  d_mx(1,1);
        ve(1,2,1) =  d_mx(1,2) - d_mx(1,1);
        ve(1,3,1) =  d_mx(2,1) - d_mx(1,1);
        ve(1,4,1) =  d_mx(2,2) - d_mx(2,1) - d_mx(1,2) + d_mx(1,1);
        
        d_my_arg = (rh*kp-ap)/sqrt(1-rh^2);
        d_my(1,1) = -normpdf(-kp)*normcdf( d_my_arg );
        d_my(1,2) = 0;
        d_my(2,1) = -normpdf(-kp);
        d_my(2,2) = 0;
        
        ve(1,1,2) =  d_my(1,1);
        ve(1,2,2) =  d_my(1,2) - d_my(1,1);
        ve(1,3,2) =  d_my(2,1) - d_my(1,1);
        ve(1,4,2) =  d_my(2,2) - d_my(1,2) - d_my(2,1) + d_my(1,1);
        
        S_aa = [1 rh; rh 1];
        d_rho = mvnpdf([-ap -kp],[0 0],S_aa);
        
        ve(1,1,3) =  d_rho;
        ve(1,2,3) = -d_rho;
        ve(1,3,3) = -d_rho;
        ve(1,4,3) =  d_rho;
        