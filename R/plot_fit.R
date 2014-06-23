#' Plot the Gaussian contours and marginals of a given parameter fit.
#'
#' @param A parameter list, as returned by gaussian_fit() 
#' @export
plot.grt <- function(fit_params) {
  bin_width=.25;
  bivar_obj = bivar_norm(fit_params, bin_width);
  dist = bivar_obj$dist;
  x = bivar_obj$x;
  y = bivar_obj$y;
  xlims = c(x[1], x[length(x)]);
  ylims = c(y[1], y[length(y)]);
  level = c(.1, .1); # Determines which contour to plot
  
  # Plot Gaussian contours 
  plot.new();
  par(fig=c(.3,1,.3,1), mar=c(2,.5,1,.5));
  for (i in 1:4) {
    cond = fit_params[[i]];
    par(new = TRUE);
    contour(x = x, y = y, z = dist[i,,], 
            levels = level, xlim = xlims, ylim = ylims,
            drawlabels = FALSE, axes = FALSE, asp=1);
    points(cond[1], cond[2], pch = '+');
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
  par(fig=c(.3,1,.1,.4), pty='m', xaxt = 'n', yaxt = 'n', new=TRUE);
  plot(margx$aa,type='l');
  lines(margx$ab,type='l',lty=2);
  lines(margx$ba,type='l',lty=2);
  lines(margx$bb,type='l');
  
  # Plot Y marginals
  par(fig=c(.1,.3,.3,1), pty='m', xaxt = 'n', yaxt = 'n', new=TRUE);
  plot(margy$aa,bivar_obj$y,type='l');
  lines(margy$ab,bivar_obj$y, type='l',lty=2);
  lines(margy$ba,bivar_obj$y,type='l',lty=2);
  lines(margy$bb,bivar_obj$y,type='l');
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