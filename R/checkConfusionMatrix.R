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
