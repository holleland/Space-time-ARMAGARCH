#' Function for estimating circular ARMAGARCH
#'
#' Wrapper function for creating a likelihood object from the TMB
#' package.
#'
#' @param data Data matrix (spatial row index, temporal col index)
#' @param W Neighbourhood 3D-array
#' @param init Initial value of sigma
#' @param parameters List of initial parameters
#' @param map List of same size as parameters, determining which parameters should be fixed.
#' @param silent Logical, indicating wheter optimization print-out should be printet.
#' @example inst/CreateLikelihood_example.R
#' @return List containing automatic differentiated likelihood function.
#'
#' @export
#' @useDynLib STARMAGARCH
#'
CreateLikelihood <- function(data, W=NULL, init=apply(data, 1, stats::var), parameters = NULL, map=NULL, silent=TRUE){
  # Controlling input:
  if(length(dim(W))<3) stop("W must be array of dimension 3.")
  if(nrow(data)!=dim(W)[1]) stop("data must have the same number of rows as W.")
  if(dim(W)[1]!=dim(W)[2]) stop("Each neighbour matrix in the array W must be a square matrix.")
  if(length(init) != nrow(data)) stop("Init must be of same length as data.")
  if(class(data)!="matrix") stop("data must be a matrix.")
  if(class(W)!="array") stop("W must be an array.")
  if(class(parameters)!="list") stop("parameters must be a list.")
  if(!all(c("mu", "phi", "theta","omega", "alpha", "beta") %in% names(parameters))){
    k <- which(!(c("mu", "phi", "theta","omega", "alpha", "beta") %in% names(parameters)))
    stop(paste("The following elements is missing from parameters: ",
               paste(c("mu", "phi", "theta","omega", "alpha", "beta")[k],collapse=", "), sep=""))
  }
  if(class(W)!="array") stop("W must be an array.")

  # Createing AD function (the QML function)
    return(TMB::MakeADFun(data = list(y=data, W=W,init = init),
                   parameters = parameters[c("mu", "phi", "theta","omega", "alpha", "beta")],
                   DLL = "STARMAGARCH", map = map, silent = silent))
}

#' Optimize likelihood using nlminb
#'
#' @param f object produced by CreateLikelihood
#' @param data Data matrix (spatial row index, temporal col index)
#' @param upper upper limits for parameter estimates
#' @param lower lower limits for parameters esimates
#' @param print logical, indicator for printing summary
#' @param simple  logical switch for turning on only returning parameter estimates and sd.
#' @export
#' @return \code{starmagarch} object.
fitSTARMAGARCH <- function(f, data=NULL,  print = TRUE, lower = NULL, upper = NULL, simple=FALSE) {
  # Setting limits if not specified by user:
  if(is.null(lower))
    lower <- c(rep(-1e4, length(which(names(f$par) == "mu"))),
      rep(-5,length(which(names(f$par)%in% c("phi","theta")))),
      rep(-1e-8, length(which(names(f$par) =="omega"))),
      rep(1e-8, length(which(names(f$par) %in% c("alpha","beta")))))
  if(is.null(upper))
    upper <- c(rep(1e4, length(which(names(f$par) == "mu"))),
              rep(5,length(which(names(f$par)%in% c("phi","theta")))),
              rep(1e4, length(which(names(f$par) == "omega"))),
              rep(1,length(which(names(f$par) %in% c("alpha","beta")))))
  #if(length(lower)!=length(upper)) stop("upper and lower limits must be of same length.")
  #if(length(lower)!=length(f$par)) stop("lower be of same length as the number of parameters.")
  #if(length(upper)!=length(f$par)) stop("upper be of same length as the number of parameters.")

  # Optimization:
  fit <- stats::nlminb(f$par,f$fn,f$gr, f$he,
                lower = lower,
                upper = upper)
  if(simple){
    return(c(fit$par, sqrt(diag(solve(f$he(fit$par))))))
  }
  # Problemer her hvis man ikke har med mu!!!
  matcoef <- data.frame(Estimates = fit$par,
                    SD = sqrt(diag(solve(f$he(fit$par)))))
  matcoef$Zscore <- matcoef$Estimates/matcoef$SD
  matcoef$Pvalue <- stats::pnorm(abs(matcoef$Zscore), lower.tail=FALSE)
  # Row names problem:
  rownames(matcoef)<-correct.names(names(f$par))
  obj <- list()
  class(obj) <- "starmagarch"
  obj$coefficients <- fit$par
  names(obj$coefficients) <-correct.names(names(f$par))
  obj$matcoef <- matcoef
  obj$hessian <- f$he()
  obj$observations <- data
  obj$fitted.values <- f$report()$yhat
  obj$sigma <- sqrt(f$report()$sigma)
  obj$garch <- f$report()$x
  obj$dim <- dim(obj$fitted.values)
  obj$optimization <- fit
  obj$aic <- fit$objective+2*length(fit$par)
  obj$bic <- fit$objective+log(prod(obj$dim))*length(fit$par)
  if(print) print(matcoef)
  return(obj)
}

#' Information criterions
#'
#' Akaike's and Bayesian information criterion of a fitted starmagarch model.
#'
#' @name aic
#' @param object \code{starmagarch} object
#' @param ... optionally more fitted model objects.
#' @return \code{AIC}: AIC of fitted model.
#' @export
AIC.starmagarch <- function(object,...) object$aic

#' @rdname aic
#' @export
BIC <- function(object, ...) UseMethod("BIC")
#' Bayesian information criterion
#'
#' @rdname aic
#' @return \code{BIC}: BIC of fitted model.
#' @export
BIC.starmagarch <- function(object, ...) object$bic


#' Methods for \code{starmagarch} objects
#'
#' Collection of generic functions for \code{starmagarch} objects.
#'
#' @param object \code{starmagarch} object
#' @param x \code{starmagarch} object
#' @param ... optionally more fitted model objects.
#' @name genfunctions
#'
NULL



#' Extract model coefficients
#'
#'
#' @rdname genfunctions
#' @return \code{coef}: Coefficients of fitted model.
#' @export
coef.starmagarch <- function(object,...) object$coefficients

#' Print starmagarch
#'
#' @rdname genfunctions
#' @return Printout
#' @export
print.starmagarch <- function(x,...){
  summary(x)
}
#' @rdname genfunctions
#' @export
sigma<- function(object,...) UseMethod("sigma")
#' Extract fitted sigma process
#'
#'
#' @rdname genfunctions
#' @return \code{sigma}: Fitted sigma process, \eqn{\{ \sigma_t(u) \}}.
#' @export
sigma.starmagarch <- function(object, ...) object$sigma

fittedgarch <- function(x) UseMethod("fittedgarch")
#' Extract fitted sigma process
#'
#' @rdname genfunctions
#' @return \code{fittedgarch}: Extract ARMA-residuals (the garch process): \eqn{\widehat y_t(u)}
#' @export
fittedgarch.starmagarch <- function(object,...) object$garch

#' Extract fitted sigma process
#'
#' @rdname genfunctions
#' @return \code{fitted}: Fitted values of the ARMA process.
#' @export
fitted.starmagarch <- function(object,...) object$fitted.values

#' Extract model residuals
#'
#'
#'
#' @rdname genfunctions
#' @return \code{residuals}:  Extract standardized residuals of fitted model:
#' \eqn{z_t(u) = x_t(u) / \sigma_t(u)}
#' @export
residuals.starmagarch <- function(object,...) fittedgarch(object)/sigma(object)

#' Plot \code{starmgarch}
#'
#' Plot and compare the fitted values to the original data.
#'
#'
#' @name plotting
#' @param x Class \code{starmagarch}
#' @param ... additional arguments.
#' @return ggplot object
#' @export
plot.starmagarch <- function(x, ...){
  # Creating a long format:
  tmp <- reshape2::melt(fitted(x))
  tmp2 <- reshape2::melt(x$observations)
  tmp$type <- "Fitted values"
  tmp2$type <- "Observations"
  tmp <- rbind(tmp,tmp2)
  #Plotting with ggplot2:
  mytheme <- ggplot2::theme(
    #legend.title = element_blank(),
    axis.line = ggplot2::element_line(),
    strip.placement = "outside",
    strip.text = ggplot2::element_text(size = 11),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank())
  p <- ggplot2::ggplot(data=transform(tmp,
                                      type = factor(type, levels = c("Observations","Fitted values"))),
                       ggplot2::aes(tmp$Var2,tmp$Var1))+ggplot2::geom_raster(ggplot2::aes(fill=tmp$value))+
    ggplot2::xlab("Time")+ggplot2::ylab("Space")+ggplot2::theme_bw()+mytheme+
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = x$coefficients["mu"])+
    ggplot2::facet_wrap(~type, ncol = 1)+
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0))+
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 1, barheight = 10, title="Value"))
  return(p)
}

#' generic function
#'
#' @param x Object
#' @param ... additional arguments
#' @export
plot_garch <- function(x,...) UseMethod("plot_garch")

#' Plot \code{starmgarch}
#'
#' Plot and compare the fitted values to the original data.
#'
#'
#' @rdname plotting
#' @return ggplot object
#' @export
plot_garch.starmagarch <- function(x,...){
  # Creating a long format:
  tmp <- reshape2::melt(fittedgarch(x))
  tmp2 <- reshape2::melt(sigma(x))
  tmp3 <- reshape2::melt(residuals(x))
  tmp$type <- "GARCH process"
  tmp2$type <- "Fitted sigma"
  tmp3$type <- "Standardized residuals"
  tmp <- rbind(tmp,tmp2,tmp3)
  #Plotting with ggplot2:
  mytheme <- ggplot2::theme(
    legend.title = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(),
    strip.placement = "outside",
    strip.text = ggplot2::element_text(size = 11),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank())
  p1 <- ggplot2::ggplot(data=tmp[tmp$type=="GARCH process",],
                       ggplot2::aes(Var2,Var1))+ggplot2::geom_raster(ggplot2::aes(fill=value))+
    ggplot2::xlab("Time")+ggplot2::ylab("Space")+
    ggplot2::theme_bw()+mytheme+
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0.0)+
    ggplot2::facet_wrap(~type, ncol = 1)+
    ggplot2::theme(axis.line.x = ggplot2::element_blank())+
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0))
  p2 <- ggplot2::ggplot(data=tmp[tmp$type=="Fitted sigma",],
                        ggplot2::aes(Var2,Var1))+ggplot2::geom_raster(ggplot2::aes(fill=value))+
    ggplot2::xlab("Time")+ggplot2::ylab("Space")+
    ggplot2::theme_bw()+mytheme+
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = mean(x$sigma))+
    ggplot2::facet_wrap(~type, ncol = 1)+
    ggplot2::theme(axis.line.x = ggplot2::element_blank())+
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0))
  p3 <- ggplot2::ggplot(data=tmp[tmp$type=="Standardized residuals",],
                        ggplot2::aes(Var2,Var1))+ggplot2::geom_raster(ggplot2::aes(fill=value))+
    ggplot2::xlab("Time")+ggplot2::ylab("Space")+
    ggplot2::theme_bw()+mytheme+
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0.0)+
    ggplot2::facet_wrap(~type, ncol = 1)+
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_continuous(expand = c(0,0))

  ggfigure <- ggpubr::ggarrange(
    p1+ggplot2::ylab("")+ggplot2::xlab("")+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.line.x=ggplot2::element_blank(),
                     plot.margin=ggplot2::unit(c(0.05,0,0,0),"cm" ),
                     axis.ticks.x=ggplot2::element_blank()),
    p2+ggplot2::ylab("")+ggplot2::xlab("")+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.line.x=ggplot2::element_blank(),
                     plot.margin=ggplot2::unit(c(0.05,0,0,0),"cm" ),
                     axis.ticks.x=ggplot2::element_blank()),
    p3+ggplot2::ylab("")+ggplot2::xlab("")+
      ggplot2::theme(plot.margin=ggplot2::unit(c(0.05,0,0,0),"cm" )), ncol=1,nrow=3,align="v",
                     common.legend = FALSE, legend = "right")
  ggpubr::annotate_figure(ggfigure,
                  bottom = ggpubr::text_grob("Time", color = "black", size = 11,vjust=-1.5),
                  left = ggpubr::text_grob("Space", color = "black", rot = 90,size=11,vjust=1.7,
                                   hjust=-.0))


}


#' Function for generating a valid list of parameters from a vector
#'
#'
#' @param par Parameter vector (numeric)
#' @param ar Integer of length 2 (spatial, temporal): Spatio-temporal order of AR part
#' @param ma Integer of length 2 (spatial, temporal): Spatio-temporal order of MA part
#' @param arch Integer of length 2 (spatial, temporal): Spatio-temporal order of ARCH part
#' @param garch Integer of length 2 (spatial, temporal): Spatio-temporal order of GARCH part
#' @export
parvector2list <- function(par = NULL, ar = c(2,1), ma = c(2,1),
                           arch = c(2,1), garch = c(2,1)){
  if(is.null(par)) par <- rep(0, times = 2+sum(prod(ar), prod(ma),prod(arch),prod(garch)))
  if(length(par) != 2+sum(prod(ar), prod(ma),prod(arch),prod(garch)))
    stop("Length of parameter vector does not match the order of the model.")
  if(length(ar)!=2) stop("ar must be of length 2.")
  if(length(ma)!=2) stop("ma must be of length 2.")
  if(length(arch)!=2) stop("arch must be of length 2.")
  if(length(garch)!=2) stop("garch must be of length 2.")
  parameter <- list()
  parameter$mu <- par[1]
  parameter$phi <- matrix(par[1+1:prod(ar)], nrow=ar[1], ncol = ar[2])
  parameter$theta <- matrix(par[1+prod(ar)+1:prod(ma)],
                            nrow=ma[1], ncol = ma[2])
  parameter$omega <- par[2+prod(ar)+prod(ma)]
  parameter$alpha <- matrix(par[2+prod(ar)+prod(ma)+1:prod(arch)],
                            nrow=arch[1], ncol = arch[2])
  parameter$beta <- matrix(par[2+prod(ar)+prod(ma)+prod(arch)+1:prod(garch)],
                            nrow=garch[1], ncol = garch[2])
  return(parameter)

}



#' Generate a valid map from parameter list
#'
#'
#' Generate a valid map list from parameter list. The map should be used to
#' specify fixed parameters and which parameters should be set equal. A fixed
#' parameter is set by setting the corresponding map element to \code{NA} while
#' two elements in the \code{parameters} list which share the same factor in the
#' map list are set to be equal.
#'
#'
#' @param parameters Parameter list
#' @return List of factors with same structure as \code{parameters}.
#' @export
#'
parameterlist2maptemplate <- function(parameters){
  map <- list()
  map$mu <- as.factor(1)
  map$phi <- matrix(as.factor(1+1:prod(dim(parameters$phi))),
                    nrow=nrow(parameters$phi),
                    ncol = ncol(parameters$phi))
  map$theta <- matrix(as.factor(1+prod(dim(parameters$phi))+1:prod(dim(parameters$theta))),
                      nrow=nrow(parameters$theta),
                      ncol = ncol(parameters$theta))
  map$omega <- as.factor(2+prod(dim(parameters$phi))+prod(dim(parameters$theta)))
  map$alpha <- matrix(as.factor(2+prod(dim(parameters$phi))+prod(dim(parameters$theta))+
                               1:prod(dim(parameters$alpha))),
                            nrow=nrow(parameters$alpha),
                            ncol = ncol(parameters$alpha))
  map$beta <- matrix(as.factor(2+prod(dim(parameters$phi))+prod(dim(parameters$theta)) +
                              prod(dim(parameters$alpha))+1:prod(dim(parameters$beta))),
                     nrow=nrow(parameters$beta),
                     ncol = ncol(parameters$beta))
  return(map)
}

#' Function for fixing parameter names
#'
#' Correct the parameter names.
#'
#' @param names vector of names which are non-unique
#' @export
#' @return Vector with correct names
correct.names <- function(names){
  if(all(table(names)==1)) return(names)
  k <- table(names)
  k <- k[which(k>1)]
  for(i in names(k)){
    names[which(names==i)] <- paste(i, 1:k[i], sep="")
  }
  names
}

#' Summary of \code{starmagarch} object.
#'
#'
#' @rdname genfunctions
#' @return \code{summary}: Summary statistics of fitted model.
#' @export
summary.starmagarch <- function(object,...) {
  print(object$matcoef)
  cat("\n\nStandardized residual standard error: ", round(sd(object$garch/object$sigma),3),".\n",
      "AIC: ", AIC.starmagarch(object),"\t BIC: ", BIC.starmagarch(object),"\n")
}
