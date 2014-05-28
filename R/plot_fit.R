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
    cond = fit_params[[i]];
    v[i,,] <- matrix(data=c(1, cond[3], cond[3], 1), ncol = 2, nrow = 2);
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
  dist = array(0, dim = c(4, leny, lenx));
  for (i in 1:4) {
    rdcov = sqrt(det(v[i,,]));
    incov = solve(v[i,,]);
    cond = fit_params[[i]];
    for (y_ind in seq(from=leny, to=1, by = -1)) {
      for (x_ind in seq(from=1, to=lenx)) {
        xy_val = c(x[x_ind], y[y_ind]);
        val_mat = xy_val - c(cond[1], cond[2]); # Normalize
        exp_term = t(val_mat) %*% incov %*% val_mat;
        dist[i,y_ind,x_ind] = (1/(2*pi*rdcov)) * exp(-0.5 * exp_term);
      }
    }
  }
  return(list(dist = dist, 
              x = x,
              y = y,
              v = v));
}