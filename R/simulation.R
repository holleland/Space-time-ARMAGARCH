#' Structuring neighbourhood array
#'
#' @param m Spatial dimension vector of length 1 or 2 (integer).
#' @param sp Max spatial order (integer).
#' @param type neighbourhood type (rook or queen)
#' @param torus Logical, indicator for circular space.
#' @export
create.neighbourhood.array <- function(m = c(5, 5), sp = 2, type = "queen", torus = TRUE){
  # Checking input:
  if(length(m) == 1) m <- c(m, 1)
  if(!is.numeric(m)) stop("m must be numeric.")
  if(!is.numeric(sp)) stop("sp must be numeric.")
  if(!(type %in% c("rook","queen"))) stop("Type must be either rook or queen")
  if(!is.logical(torus)) stop("Torus must be logical.")

  # Neighbourmatrix:
  W <- spdep::cell2nb(m[1], m[2], type = type, torus = torus)
  # Neighbour array:
  Warr <- array(0, dim = c(prod(m),prod(m), sp+1))
  Warr[,,1] <- diag(prod(m))
  if(sp == 0) return(Warr)
  W <- spdep::nblag(W, max(sp,2))
  for(i in 1:(sp))
    Warr[,,i+1] <- spdep::nb2mat(W[[i]])
  return(Warr)
}


#' Simulation of STARMAGARCH
#'
#'
#' @param parameters List of parameters.
#' @param n Temporal dimension (integer).
#' @param m Spatial dimension vector of length 1 or 2 (integer).
#' @param W Neighbourhood array (optional)
#' @param burnin Temporal burnin.
#' @param type neighbourhood type (rook or queen)
#' @param torus Logical, indicator for circular space.
#' @export
#'
simSTARMAGARCH <- function(parameters, n=100, m= c(5, 5), W = NULL, burnin = 0,
                           type = "rook", torus = TRUE){
  # Controlling that input makes sense:
  if(length(m)==1) m <- c(m,1)
  if(!(type %in% c("rook","queen"))) stop("Type must be either rook or queen")
  if(length(n)>1) n <- n[1]
  if(!is.numeric(n)) stop("n must be numeric.")
  if(!is.numeric(m)) stop("m must be numeric.")
  if(!is.numeric(burnin)) stop("burnin must be numeric.")
  if(!is.logical(torus)) stop("Torus must be logical.")
  if(class(parameters)!="list") stop("parameters must be a list.")
  if(!all(c("mu", "phi", "theta","omega", "alpha", "beta") %in% names(parameters))){
    k <- which(!(c("mu", "phi", "theta","omega", "alpha", "beta") %in% names(parameters)))
    stop(paste("The following elements is missing from parameters: ",
               paste(c("mu", "phi", "theta","omega", "alpha", "beta")[k],collapse=", "), sep=""))
  }


  # Setting up neighbourhood array:
  if(is.null(W)){
    sp <- max(sapply(c(2,3,5,6), function(i) nrow(parameters[[i]])))
    W <- create.neighbourhood.array(m = m, sp = sp, type=type, torus=torus)
  }
  # ____________________________________________________________________
  #                            Simulation:
  # ____________________________________________________________________
  x <-  y <- matrix(0, ncol = n+burnin, nrow = prod(m[1:2]))
  sig <- x + parameters$omega
  # STGARCH:
  for(t in 2:ncol(y)){
    sig[,t]<-parameters$omega +

      rowSums(sapply(1:min(ncol(parameters$alpha),t-1), function(j) {
        if(nrow(parameters$alpha) ==1 ) {
          W[,,1] %*% x[,t-j]^2 * parameters$alpha[1, j]
        }else{
          apply(W[,,1:nrow(parameters$alpha)], 3, function(w) w %*% x[,t-j]^2) %*% parameters$alpha[, j]}
        }
        ))+
      rowSums(sapply(1:min(ncol(parameters$beta),t-1), function(j) {
        if(nrow(parameters$beta) ==1 ){
          W[,,1] %*% sig[,t-j] * parameters$beta[1, j]
        } else {
          apply(W[,,1:nrow(parameters$beta)], 3, function(w) w %*% sig[,t-j]) %*% parameters$beta[, j]}
        }))
    x[,t]<-sqrt(sig[,t])*rnorm(prod(m))
  }
  # STARMA:
  for(t in 2:ncol(y)){
    y[,t]<-parameters$mu +
      rowSums(sapply(1:min(ncol(parameters$phi),t-1), function(j) {
        if(nrow(parameters$phi) ==1 ){ W[,,1] %*% (y[,t-j]-parameters$mu)*parameters$phi[1, j]
        }else{ apply(W[,,1:nrow(parameters$phi)], 3, function(w)
          w %*% (y[,t-j]-parameters$mu))%*% parameters$phi[, j]}
        }))+
      rowSums(sapply(1:min(ncol(parameters$theta),t-1), function(j){
        if(nrow(parameters$theta) ==1 ) {
          W[,,1] %*% x[,t-j] * parameters$theta[1, j]
        } else{apply(W[,,1:nrow(parameters$theta)], 3, function(w) w %*% x[,t-j]) %*% parameters$theta[, j]
        }
        }))+x[,t]
  }
  #Return data with burnin removed:
  return(y[,-(1:burnin)])
}
