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
#' @export
#' @examples
#' mri_test(observerA)
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
    statistic[B+2] <- ((rB.sA1B/nA1B - rB.sA2B/nA2B)/
                      sqrt(p.s*(1-p.s)*(1/nA1B+1/nA2B)) )
  }
  return(data.frame(stimulus=stimulus, statistic=statistic, 
      p.value= 2*(pmin(1-pnorm(statistic),pnorm(statistic))) ))
}
