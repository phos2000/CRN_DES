library(pracma)

### use point estimate and 95% CI to get the distribution
### basically1 it is (qdist() - lower)^2 - (qdist() - upper)^2

randomNorm = function(lower, upper, n){
  # norm_obj = function(parameters){
  #   mean = parameters[1]
  #   sd = parameters[2]
  #   return((qnorm(0.025, mean, sd) - lower)^2 + (qnorm(0.975, mean, sd) - upper)^2)
  # }
  # 
  # result = optim(c(0,1), norm_obj)
  
  # Transfinite scaling using logit
  Tf <- function(x, mu, sigma)
    pnorm(x,mu,sigma,log=TRUE) - 
    pnorm(x,mu,sigma,log=TRUE, lower.tail = FALSE)
  
  # Function to minimize (note log(p1) is okay as this is constant).
  norm <- function(x1, p1, x2, p2, mu, sigma)
    (Tf(x1, mu, sigma)-(log(p1)-log(1-p1)) )^2 +
    (Tf(x2, mu, sigma)-(log(p2)-log(1-p2)) )^2
  
  # Bundle it all together into a single function.
  fn <- function(x) norm(lower, 0.025, upper, 0.975, x[1], x[2])
  result = 
    optim(c(0.5, 0.1),
        fn,
        gr = function(x) pracma::grad(fn, x),
        lower=c(-Inf, 1e-4),
        method = "L-BFGS-B",
        control=list(factr=1e-10, maxit=100))
  
  mean = result$par[1]
  sd = result$par[2]
  
  return(rnorm(n, mean, sd))
}

randomBeta = function(lower, upper, n){
  # beta_obj <- function(parameters) {
  #   alpha <- parameters[1]
  #   beta <- parameters[2]
  #   return((qbeta(0.025, alpha, beta) - lower)^2 + (qbeta(0.975, alpha, beta) - upper)^2)
  #          #+ 2*(alpha / (alpha + beta) - mean)^2
  # }
  # 
  # # Optimization
  # result <- optim(c(10, 10), beta_obj)
  
  # Transfinite scaling using logit
  Tf <- function(x, alpha, beta)
    pbeta(x,alpha,beta,log=TRUE) - 
    pbeta(x,alpha,beta,log=TRUE, lower.tail = FALSE)
  
  # Function to minimize (note log(p1) is okay as this is constant).
  norm <- function(x1, p1, x2, p2, alpha, beta)
    (Tf(x1, alpha, beta)-(log(p1)-log(1-p1)) )^2 +
    (Tf(x2, alpha, beta)-(log(p2)-log(1-p2)) )^2
  
  # Bundle it all together into a single function.
  fn <- function(x) norm(lower, 0.025, upper, 0.975, x[1], x[2])
  result = 
    optim(c(10, 1),
          fn,
          gr = function(x) pracma::grad(fn, x),
          lower=c(-Inf, 1e-4),
          method = "L-BFGS-B",
          control=list(factr=1e-10, maxit=100))
  
  alpha <- result$par[1]
  beta <- result$par[2]
  
  return(rbeta(n, alpha, beta))
}

randomGamma = function(lower, upper, n){
  # Known values: Mean value of the distribution, 2.5% quantile, 97.5% quantile
  # Objective function to minimize
  # gamma_obj <- function(parameters) {
  #   shape <- parameters[1]
  #   rate <- parameters[2]
  #   return((qgamma(0.025, shape, rate) - lower)^2 + (qgamma(0.975, shape, rate) - upper)^2)
  #          #+ 2*(shape / rate - mean)^2
  # }
  # 
  # # Optimization
  # result <- optim(c(1, 1), gamma_obj)
  
  # Transfinite scaling using logit
  Tf <- function(x, shape, rate)
    pgamma(x,shape,rate,log=TRUE) - 
    pgamma(x,shape,rate,log=TRUE, lower.tail = FALSE)
  
  # Function to minimize (note log(p1) is okay as this is constant).
  norm <- function(x1, p1, x2, p2, shape, rate)
    (Tf(x1, shape, rate)-(log(p1)-log(1-p1)) )^2 +
    (Tf(x2, shape, rate)-(log(p2)-log(1-p2)) )^2
  
  # Bundle it all together into a single function.
  fn <- function(x) norm(lower, 0.025, upper, 0.975, x[1], x[2])
  result = 
    optim(c(1, 1),
          fn,
          gr = function(x) pracma::grad(fn, x),
          lower=c(-Inf, 1e-4),
          method = "L-BFGS-B",
          control=list(factr=1e-10, maxit=100))
  
  shape <- result$par[1]
  rate <- result$par[2]

  return(rgamma(n,shape = shape,rate = rate))
}

randomLnorm = function(lower, upper, n){
  # Known values: Mean value of the distribution, 2.5% quantile, 97.5% quantile
  # Objective function to minimize
  # lognorm_obj <- function(parameters) {
  #   meanlog <- parameters[1]
  #   sdlog <- parameters[2]
  #   return((qlnorm(0.025, meanlog, sdlog) - lower)^2 + (qlnorm(0.975, meanlog, sdlog) - upper)^2)
  #            #+ 2*(exp(meanlog + (sdlog^2)/2) - mean)^2
  # }
  # 
  # # Optimization
  # result <- optim(c(1, 1), lognorm_obj)
  
  # Transfinite scaling using logit
  Tf <- function(x, meanlog, sdlog)
    plnorm(x,meanlog,sdlog,log=TRUE) - 
    plnorm(x,meanlog,sdlog,log=TRUE, lower.tail = FALSE)
  
  # Function to minimize (note log(p1) is okay as this is constant).
  norm <- function(x1, p1, x2, p2, meanlog, sdlog)
    (Tf(x1, meanlog, sdlog)-(log(p1)-log(1-p1)) )^2 +
    (Tf(x2, meanlog, sdlog)-(log(p2)-log(1-p2)) )^2
  
  # Bundle it all together into a single function.
  fn <- function(x) norm(lower, 0.025, upper, 0.975, x[1], x[2])
  result = 
    optim(c(1, 1),
          fn,
          gr = function(x) pracma::grad(fn, x),
          lower=c(-Inf, 1e-4),
          method = "L-BFGS-B",
          control=list(factr=1e-10, maxit=100))
  
  meanlog <- result$par[1]
  sdlog <- result$par[2]

  return(rlnorm(n, meanlog, sdlog))
}

### random drawing

randomDraw = function(param, params, n){
  # param: the parameter to be drawed
  # params: the source of all parameters
  # n: how many times drawing
  params_to_use <- params[grepl(param,names(params))]
  if (params_to_use[[paste0(param, "_dist")]] == "beta") {
    return(randomBeta(params_to_use[[paste0(param, "_low")]], params_to_use[[paste0(param, "_upp")]], n))
  } else if (params_to_use[[paste0(param, "_dist")]] == "gamma") {
    return(randomGamma(params_to_use[[paste0(param, "_low")]], params_to_use[[paste0(param, "_upp")]], n))
  } else if (params_to_use[[paste0(param, "_dist")]] == "lognormal"){
    return(randomLnorm(params_to_use[[paste0(param, "_low")]], params_to_use[[paste0(param, "_upp")]], n))
  } else if (params_to_use[[paste0(param, "_dist")]] == "normal") {
    return(randomNorm(params_to_use[[paste0(param, "_low")]], params_to_use[[paste0(param, "_upp")]], n))
  } else {stop("No such distribution!")}
}