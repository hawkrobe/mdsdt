plot_fit <- function(fit_params) {
  lab.x = 'x';
  lab.y = 'y';
  bin_width = .25;
  bounds = TRUE;
  marginals = TRUE;
  plot_fit = TRUE;
  
  bivar_obj = bivar_norm(fit_params, bin_width);
  xlims = with(bivar_obj, c(x[1], x[length(x)]));
  ylims = with(bivar_obj, c(y[1], y[length(y)]));
  level = c(.1, .1); # Determines which contour to plot
  minax = -2.5;
  maxax = 4.5;
  plot.new();
  par(fig=c(0,1,0,1), pty='s', new=TRUE);
  # Add axes
  if (marginals || plot_fit) {
    for (i in 1:4) {
      cond = fit_params[[i]];
      print(i);
      par(fig=c(0,1,0,1), pty='s', new=TRUE);
      contour(x = bivar_obj$x, y = bivar_obj$y, z = bivar_obj$dist[i,,], 
              levels = level, xlim = xlims, ylim = ylims,
              drawlabels = FALSE, axes = FALSE);
      points(cond[1], cond[2], pch = '+');
    }
    abline(h=0);
    abline(v=0);
  }
}

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