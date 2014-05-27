plot_fit <- function(fit_params) {
  lab.x = 'x';
  lab.y = 'y';
  bin_width = .05;
  bounds = TRUE;
  marginals = TRUE;
  plot_fit = TRUE;
  

}

bivar_norm <- function(fit_params, bin_width) {
  v <- array(0, dim = c(4,2,2));
  for (i in 1:4) {
    cond = names(fit_params)[i];
    v[i,,] <- matrix(data=c(1,fit_params$cond[3], fit_params$cond[3], 1), ncol = 2, nrow = 2);
  }
  densfact = 3; # How many SDs to compute?
  xrange <- array(0, dim = c(4,2));
  yrange <- array(0, dim = c(4,2));
  for (i in 1:4) {
    cond = fit_params[[i]];
    x_sd = densfact*sqrt(v[i,1,1]);
    y_sd = densfact*sqrt(v[i,2,2]);
    xrange[i,] = c(cond[1] - x_sd, cond[1] + x_sd);
    yrange[i,] = c(cond[2] - y_sd, cond[2] + y_sd);
  }
  x = seq(from = min(xrange), to = max(xrange), by = bin_width);
  y = seq(from = min(yrange), to = max(yrange), by = bin_width);
  lenx = length(x);
  leny = length(y);
  return(xrange);
}