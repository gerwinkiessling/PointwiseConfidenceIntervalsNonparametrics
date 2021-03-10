##### On Pointwise Confidence Intervals in Nonparametric Regression#####
##### Gerwin Janis Kiessling #####

### This file includes the functions.


#Installing some packages
install.packages("sfsmisc")
library(sfsmisc)

install.packages("mvtnorm")
library(mvtnorm)

install.packages("expm")
library(expm)

install.packages("ggplot2")
library(ggplot2)

install.packages("pracma")
library(pracma)

install.packages('latex2exp')
library(latex2exp)

install.packages("readxl")
library(readxl)

install.packages("matrixStats")
library(matrixStats)

install.packages("foreign")
library(foreign)

epa_kernel<- function(v){
  
  # The Epanechnikov Kernel
  
  #Arguments: v: point of evaluation
  #Returns: K: value of kernel function
  
  K<-3/4*(1-v^2)*(-1<=v)*(v<=1)
  return(K)
}

density_beta<-function(x, p1, p2, ncp){
  
  #Pdf of the beta distribution
  
  #Arguments: - x: point in the support
  #           - p1: parameter 1
  #           - p2: parameter 2
  #           - ncp: non-centrality parameter
  #Returns:   - f: value of the pdf
  
  f=dbeta(x, p1, p2, ncp=ncp)
  return(f)
}

data_generating_processes<-function(n, label, p1, p2, e, r, mu, var){
  
  #Data generating processes
  
  #Arguments: - n: sample size
  #           - label: string for number of simulation
  #           - p1: parameter 1 of beta distribution
  #           - p2: parameter 2 of beta distribution
  #           - e: lower end of uniform distribution
  #           - r: upper end of uniform distribution
  #Returns:   - x: vector of x data
  #           - g_x: vector of true conditional mean function value
  #           - Y: vector outcome variable data
  
  
  if (label=='sim1'){
    x<-rbeta(n, p1, p2)
    g_x <- cos(3*pi*x)
    Y<-g_x+runif(n, min=e, max=r)
  }
  
  if (label=="sim2"){
  x<-rbeta(n, p1, p2)
  g_x <- 0.5*sin(4*pi*x)+dnorm(x, mean=mu, sd=sqrt(var))
  Y<-g_x+runif(n, min=e, max=r)
  }
  
  if (label=="sim4con"){
    x<-rbeta(n, p1, p2)
    g_x <- x*sin(4/(x+0.01))
    Y<-g_x+(exp(x)+0.1)*runif(n, min=e, max=r)
  }
  
  if (label=="sim4lin"){
    x<-rbeta(n, p1, p2)
    g_x <- x*sin(4/(x+0.01))
    Y<-g_x+sqrt((exp(x)+0.1))*runif(n, min=e, max=r)
  }
  
  if (label=="sim5con"){
    x<-rbeta(n, p1, p2)
    g_x <- -sin(2*pi+3*pi*x)
    Y<-g_x+sqrt(x^4+0.5)*runif(n, min=e, max=r)
  }
  
  if (label=="sim5lin"){
    x<-rbeta(n, p1, p2)
    g_x <- -sin(2*pi+3*pi*x)
    Y<-g_x+sqrt(x^4+0.5)*runif(n, min=e, max=r)
  }
  
  return(list(x=x, Y=Y, g_x=g_x))
}

g_x_calc<-function(x, label, mu, var){
  
  #Calculates true conditional mean function
  
  #Arguments: - x: regressor data
  #           - label: label for simulation
  #           - mu: mean of normal distribution for sim2
  #           - var: variance of normal distribution for sim2
  #Returns:   - g_x: true conditional mean function value
  
  if(label=="sim1"){
    g_x=cos(3*pi*x)
  }
  if(label=="sim2"){
    g_x=0.5*sin(4*pi*x)+dnorm(x, mean=mu, sd=sqrt(var))
  }
  
  if (label=="sim4con"){
    g_x <- x*sin(4/(x+0.01))
  }
  
  if (label=="sim4lin"){
    g_x <- x*sin(4/(x+0.01))
  }
  
  if (label=="sim5con"){
    g_x <- -sin(2*pi+3*pi*x)
  }
  
  if (label=="sim5lin"){
    g_x <- -sin(2*pi+3*pi*x)
  }
  return(g_x)
}

calculate_bias<-function(label, x, p1, p2, bandwidth, mu, var){
  
  #Calculates the bias
  
  #Arguments:  - label: string for simulation number
  #            - x: regressor data point
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - bandwidth: bandwidth
  #Returns:    - bias: bias at point x
  #            - derivative beta pdf: derivative of beta distribution at point x
  
    derivative_beta_pdf<-(1/(gamma(p1)*gamma(p2)/gamma(p1+p2)))*((p1-1)*x^(p1-2)*(1-x)^(p2-1)-x^(p1-1)*(1-x)^(p2-2)*(p2-1))
    if(label=="sim1"){
    derivative_g_x <- -3*pi*sin(3*pi*x)
    second_derivative_g_x <- -9*pi^2*cos(3*pi*x)
    bias <- bandwidth^2*(3/15)*(1/2)*(2*derivative_beta_pdf*derivative_g_x/dbeta(x, p1, p2)+second_derivative_g_x)
    }
    
    if(label=="sim2"){
      second_derivative_g_x <- -8*pi^2*sin(4*pi*x)+((x-mu)^2-var)/(sqrt(2*pi)*var^(5/2))*exp(-(x-mu)^2/(2*var))
      bias <- bandwidth^2*(3/15)*(1/2)*second_derivative_g_x
    }
    
    if(label=="sim4con"){
      
      derivative_g_x <- sin(4/(x+0.01))-4*x*cos(4/(x+0.01))/(x+0.01)^2
      second_derivative_g_x <- -80000*(20000*x*sin(4/(x+0.01))+(100*x+1)*cos(4/(x+0.01)))/(100*x+1)^4
      bias <- bandwidth^2*(3/15)*(1/2)*(2*derivative_beta_pdf*derivative_g_x/dbeta(x, p1, p2)+second_derivative_g_x)
    }
    
    if(label=="sim4lin"){
      second_derivative_g_x <- -80000*(20000*x*sin(4/(x+0.01))+(100*x+1)*cos(4/(x+0.01)))/(100*x+1)^4
      bias <- bandwidth^2*(3/15)*(1/2)*second_derivative_g_x
    }
    
    if(label=="sim5con"){
      derivative_g_x <- -3^2*pi^2*cos(2*pi+3*pi*x)#-1/(sqrt(2*pi*var^3))*(x-mu)*exp(-(x-mu)^2)/(2*var)
      second_derivative_g_x <- 3^3*pi^3*sin(1.25*pi+2.5*pi*x)
      bias <- bandwidth^2*(3/15)*(1/2)*(2*derivative_beta_pdf*derivative_g_x/dbeta(x, p1, p2)+second_derivative_g_x)
    }
    
    if(label=="sim5lin"){
      second_derivative_g_x <- -2.5^3*pi^3*sin(1.25*pi+2.5*pi*x)
      bias <- bandwidth^2*(3/15)*(1/2)*second_derivative_g_x
    }
    
    return(list(bias=bias, derivative_beta_pdf=derivative_beta_pdf))
}


cv_local_constant<-function(kernel_function, x, Y, bw_seq){
  # Performs least-squares cross-validation for local constant regression
  
  # Arguments: - kernel_function: the kernel function
  #            - x: vector of x data
  #            - Y: vector of Y data
  #            - bw_seq: bandwidth sequence
  # Returns:   - h_loc_opt: bandwidth which minimizes the cv criterion
  
  cv_criterion_loc<-c()
  
  for (m in 1:length(bw_seq)){
    #iterating over the bandwidth grid
    h<-bw_seq[m]
    #denominator and numerator of weights
    denominator = rowSums( kernel_function( (t(replicate(length(x), x))-x)/h ))
    numerator <- kernel_function((t(replicate(length(x), x))-x)/h)
    # ii entry of smoothing matrix
    L_ii<-kernel_function(0)/denominator
    # nadaraya watson estimator at xi
    nw_i<-rowSums(numerator*t(replicate(length(x),Y) ) )/denominator
  
    #taking only observations away from the boundaries
    Y_subset<-Y[h<x&x<(1-h)]
    L_ii_subset<-L_ii[h<x&x<(1-h)]
    nw_i_subset<-nw_i[h<x&x<(1-h)]
    
    #formula for cross validation from Wasserman (2006, p.70) 
    cv_criterion_loc[m] <- mean(((Y_subset-nw_i_subset)/(1-L_ii_subset))^2, na.rm=TRUE)
  }
  index<-which.min(cv_criterion_loc)
  h_loc_opt=bw_seq[index]
  return(h_loc_opt)
}

local_constant<-function(kernel_function, bandwidth, x_vector, evaluation_vector){
  #Calculates numerator and denominator of weights for the local constant estimator
  
  # Arguments: - kernel_function: the kernel function
  #            - bandwidth: bandwidth
  #            - x_vector: vector of x data points
  #            - evaluation_vector: grid of points where estimator is evaluated
  # Returns:   - numerator: matrix with elements used for numerators of weights
  #            - denominator: vector with denominator of weights
  
  numerator <- kernel_function((t(replicate(length(evaluation_vector), x_vector))-evaluation_vector)/bandwidth)
  denominator <- rowSums(kernel_function((t(replicate(length(evaluation_vector), x_vector))-evaluation_vector)/bandwidth))
  
  return(list(numerator=numerator, denominator=denominator))
}

local_constant_estimator<-function(kernel_function, label, n, p1, p2, mu, var, x_seq, e, r, undersmooth_seq){
  #Perform local constant kernel estimation
  
  # Arguments: - kernel_function: kernel funtion
  #            - label: string for simulation number
  #            - n: sample size
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - x_seq: points of evaluation
  #            - e: lower end of uniform distribution
  #            - r: upper end of uniform distribution
  #            - undersmooth_seq: the sequence of undersmoothing factors
  # Returns:   - x: x data
  #            - Y: y data
  #            - gh_x: estimated function values at grid points
  #            - gh_x_i: estimated function values at x_i points
  #            - smoothing list: list of smoothing matrices for the grid points
  #            - trace_L: trace of smoothing matrix accounting for boundary problems
  #            - vtilde: trace of smoothing matrix transposed*smoothing matrix accounting for boundary problems
  #            - bandwidths: sequence of bandwidths used for estimation
  
  dgp<-data_generating_processes(n=n, label=label, p1=p1, p2=p2, e=e, r=r, mu=mu, var=var)
  x<-dgp$x
  Y<-dgp$Y
  #thumb bandwidth with epanechnikov kernel
  h_baseline<- sd(x)*n^(-1/5)#thumbBw(x, Y, 0, EpaK)
  #reduced bandwidth sequence
  red_h_seq <- seq(from=h_baseline*0.5, to=h_baseline*1.5, length.out=7)
  h_cv<-cv_local_constant(kernel_function=kernel_function, x=x, Y=Y, bw_seq=red_h_seq)
  
  #initializing objects
  gh_x <- matrix(NA, ncol=length(undersmooth_seq), nrow=length(x_seq))
  gh_x_i<-matrix(NA, ncol=length(undersmooth_seq), nrow=n)
  smoothing_list<-list()
  smoothing_list_i<-list()
  trace_L<-c()
  vtilde<-c()
  bandwidths<-c()
  

  for (f in 1:length(undersmooth_seq)){
    smoothing_matrix <- matrix(NA, nrow=length(x), ncol=length(x))
    smoothing_matrix_grid <- matrix(NA, nrow=length(x_seq), ncol=length(x))
    
    factor_undersmooth<-undersmooth_seq[f]
    h<-n^(factor_undersmooth)*h_cv
    
    local_constant_results<-local_constant(kernel=kernel_function, bandwidth=h, x_vector=x, evaluation_vector = x)
    
    denominator <-local_constant_results$denominator
    numerator <- local_constant_results$numerator
    
    for (i in 1:length(x)){
      for (j in 1:length(x)){
        smoothing_matrix[i,j]<-numerator[i,j]/denominator[i]
      }
    }
    booleans<- ((h<x)&(x<(1-h)))
    # I do not want to evaluate the weights 
    # at the boundary, but I do want to include the effect of boundary 
    # observations on weights of observations
    # in the interior of the support as these
    # were used to estimate g. I checked in simulations 
    # that the variance is consistently estimated.
    smoothing_matrix_trace_LT_L<-smoothing_matrix[booleans,]
    smoothing_matrix_reduced<-smoothing_matrix[booleans,booleans]
    gh_x_i[,f]=rowSums(numerator*t(replicate(length(x),Y)))/denominator
    
    
    local_constant_results_grid<-local_constant(kernel=kernel_function,bandwidth=h, x_vector=x, evaluation_vector = x_seq)
    denominator_grid<-local_constant_results_grid$denominator
    numerator_grid<-local_constant_results_grid$numerator
    
    for (i in 1:length(x_seq)){
      for (j in 1:length(x)){
        smoothing_matrix_grid[i,j]<-numerator_grid[i,j]/denominator_grid[i]
      }
    }
    
    gh_x[,f]=rowSums(numerator_grid*t(replicate(length(x_seq),Y)))/denominator_grid
    
    trace_L[f]<-sum(diag(smoothing_matrix_reduced))
    vtilde[f]<-sum(diag(crossprod(smoothing_matrix_trace_LT_L, smoothing_matrix_trace_LT_L)))
    smoothing_list[[f]]<-smoothing_matrix_grid
    smoothing_list_i[[f]]<-smoothing_matrix
    
    bandwidths[f]<-h
  }
  return(list(x=x, Y=Y, gh_x=gh_x, gh_x_i=gh_x_i, smoothing_list=smoothing_list, 
              smoothing_list_i=smoothing_list_i, trace_L=trace_L, vtilde=vtilde, bandwidths=bandwidths))
}

error_variance_estimation<-function(x, Y, gh_x_i, smoothing_matrix, smoothing_matrix_grid, kernel_function){
  
  #calculating standard error of g using the local constant estimator
  
  #arguments: - x: x values
  #           - Y: y values
  #           - gh_x_i: estimated values of regression function at xi
  #           - smoothing_matrix_grid: smoothing matrix at xi values
  #           - kernel_function: the kernel function
  #returns:   - s_x: standard error of g hat
  
  loo_squared_errors<- (((Y-gh_x_i)/(1-diag(smoothing_matrix)))^2)
  h_baseline<- sd(x)*length(x)^(-1/5)
  red_h_seq<-seq(from=h_baseline*0.5, to=h_baseline*1.5, length.out=5)
  h<-cv_local_constant(kernel_function=kernel_function, x=x, Y=loo_squared_errors,
                    bw_seq=red_h_seq)
  
  local_constant_results<-local_constant(kernel=kernel_function, bandwidth=h, x_vector=x, evaluation_vector = x)
  
  denominator <-local_constant_results$denominator
  numerator <- local_constant_results$numerator
  
  sigma_x_i<-rowSums(numerator*t(replicate(length(x),loo_squared_errors)))/denominator
  v_x<-as.vector(sigma_x_i%*%t(smoothing_matrix_grid^2))
  s_x <- sqrt(v_x)
  
  return(s_x)
}


pointwise_intervals_local_constant<-function(label, kernel_function, n, p1, p2, 
                                             mu, var, x_seq, e, r, alpha, undersmooth_seq, 
                                             variance_estimator){
  
  #calculates pointwise confidence intervals for local constant regression
  
  ##arguments: - label: string for the simulation number
  #            - kernel_function: kernel_function
  #            - n: sample size
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - x_seq: sequence of x values 
  #            - e: lower end of uniform distribution
  #            - r: upper end of uniform distribution
  #            - undersmooth_seq: sequence of undersmoothing factors
  #            - variance_estimator: string that specifies whether the robust
  #               or non-robust estimator is used
  
  ##returns:   - coverage_matrix: matrix with booleans for coverage 
  #            - mse_matrix: matrix with pointwise MSEs
  
  g_x<-g_x_calc(x=x_seq, label=label, mu=mu, var=var)
  
  results<-local_constant_estimator(kernel_function=kernel_function, label=label,  
                                    n=n, p1=p1, p2=p2, mu=mu, var=var, x_seq=x_seq, e=e, r=r, 
                                    undersmooth_seq=undersmooth_seq)
  
  coverage_matrix<-matrix(NA, nrow=length(x_seq), ncol=length(undersmooth_seq)*2)
  mse_matrix<-matrix(NA, nrow=length(x_seq), ncol=length(undersmooth_seq))
  
  for (f in 1:length(undersmooth_seq)){
    results_df<-data.frame('gh_x'=results$gh_x[,f], 'g_x'=g_x)
    trace_L<- results$trace_L[f]
    
    l_norm<-sqrt(rowSums(results$smoothing_list[[f]]^2, na.rm = TRUE))
    
    Y_subset<-results$Y[((results$bandwidths[f]<results$x)&(results$x<(1-results$bandwidths[f])))]
    gh_x_i_subset<-results$gh_x_i[((results$bandwidths[f]<results$x)&(results$x<(1-results$bandwidths[f]))),f]
    
    #calculating standard errors
    if (variance_estimator=="non-robust"){
    standard_error_not_robust<-sqrt(sum((Y_subset-gh_x_i_subset)^2)/(length(Y_subset)-2*trace_L+results$vtilde[f]))
    se<-standard_error_not_robust*l_norm
    #standard_deviation_true<-sqrt((1/12)*(r-e)^2)
    
    #s_true_variance<- sqrt((1/12)*(r-e)^2)*l_norm
    }
    
    if (variance_estimator=="robust"){
      se<-error_variance_estimation(x=results$x, Y=results$Y, gh_x_i=results$gh_x_i[,f],
                                smoothing_matrix=results$smoothing_list_i[[f]], 
                                smoothing_matrix_grid=results$smoothing_list[[f]], 
                                kernel_function=kernel_function)
      #standard_deviation_true<-sqrt((1/12)*(r-e)^2)
      
      #s_true_variance<- sqrt((1/12)*(r-e)^2)*l_norm
    }
    
    bias_correction<-calculate_bias(label=label, x=x_seq, p1=p1, p2=p2, bandwidth=results$bandwidths[f], 
                                    mu=mu, var=var)$bias
    E_gh<-g_x+bias_correction
    
    quantile<-qnorm(1-alpha/2)
    
    CI_lower<-results_df$gh_x-quantile*se
    CI_upper<-results_df$gh_x+quantile*se
    
    #getting coverage and MSE matrices
    
    for (i in 1:length(x_seq)){
    coverage_matrix[i,f] <- (CI_lower[i]<g_x[i]&g_x[i]<CI_upper[i])
    coverage_matrix[i, length(undersmooth_seq)+f] <- (CI_lower[i]<E_gh[i]&E_gh[i]<CI_upper[i])
    }
    
    for (i in 1:length(x_seq)){
      mse_matrix[i,f] <- (results_df$gh_x[i]-results_df$g_x[i])^2
    }

  }
  return(list(coverage_matrix=coverage_matrix,
              mse_matrix=mse_matrix))
}



calc_S1_S2_func<-function(kernel_function, x, xeval, h){
  #Calculating terms for local linear regression 
  #using Wasserman's (2007, p.77) S1, S2 representation. 
  
  # Arguments:- kernel_function: vector of data points
  #           - x: vector of data points
  #           - xeval: fixed x
  #           - h: bandwidth
  # Returns:  - b: numerator (and denominator) values for of effective smoother
  #           - S1: S1 part of effective smoother
  #           - S2: S2 part of effective smoother
  #           - l: effective smoother
  
  
  S2_ind<-c()
  for (j in 1:length(x)){
    S2_ind[j]<-kernel_function((x[j]-xeval)/h)*(x[j]-xeval)^2
  }
  S2<-sum(S2_ind)
  
  S1_ind<-c()
  for (j in 1:length(x)){
    S1_ind[j]<-kernel_function((x[j]-xeval)/h)*(x[j]-xeval)
  }
  S1<-sum(S1_ind)
  
  b<-c()
  for (i in 1:length(x)){
    b[i] <- kernel_function((x[i]-xeval)/h)*(S2-(x[i]-xeval)*S1)
  }
  
  l<-c()
  for (i in 1:length(x)){
    l[i]<-b[i]/sum(b)
  }
  l_norm<-sqrt(sum(l^2))
  
  return(list(b=b, l_norm=l_norm, S1=S1, S2=S2))
}


cv_local_linear<-function(kernel_function, x, Y, h){
  #performs least-squares cross-validation for the local linear estimator
  
  #Arguments:  - kernel_function: kernel function 
  #            - x: vector with x values
  #            - Y: vector with Y values
  #            - h: vector with bandwidths
  
  ##returns:   - h_lin_opt: optimal bandwidth
  
  cv_criterion<-c()
  for (m in 1:length(h)){
    L_ii<-c()
    gh_x<-c()
    for (k in 1:length(x)){
      xeval=x[k]
      # b is a vector with values for the numerator of the effective smoother
      b<-calc_S1_S2_func(kernel_function=kernel_function, x=x, xeval=xeval, h=h[m])$b
      l<-c()
      for (i in 1:length(x)){
        l[i]<-b[i]/sum(b)
        if(i==k){
          #making sure I use leave-one-out errors
          L_ii[k]=l[i]
        }
      }
      gh_x[k]<-l%*%Y
    }
    cv_criterion[m] <- mean(((Y-gh_x)/(1-L_ii))^2)
  }
  index<-which.min(cv_criterion)
  h_lin_opt=h[index]
  return(h_lin_opt=h_lin_opt)
}

local_linear_estimator<-function(kernel_function, label, n, p1, p2, mu, var, x_seq, e, r, undersmooth_seq){
  #implements the local linear estimator for a given data generating process
 
  # Arguments: - kernel_function: kernel funtion
  #            - label: string for simulation number
  #            - n: sample size
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - mu: normal density mean for sim 2
  #            - var: normal density variance for sim2
  #            - x_seq: points of evaluation
  #            - e: lower end of uniform distribution
  #            - r: upper end of uniform distribution
  #            - undersmooth_seq: the sequence of undersmoothing factors
  # Returns:   - x: x data
  #            - Y: y data
  #            - gh_x: estimated function values at grid points
  #            - gh_x_i: estimated function values at x_i points
  #            - smoothing list: list of smoothing matrices for the grid points
  #            - trace_L: trace of smoothing matrix
  #            - vtilde: trace of smoothing matrix transposed*smoothing matrix
  #            - bandwidths: sequence of bandwidths used for estimation
  
  dgp<-data_generating_processes(n=n, label=label, p1=p1, p2=p2, e=e, r=r, mu=mu, var=var)
  x<-dgp$x
  Y<-dgp$Y
  
  h_baseline<- sd(x)*length(x)^(-1/5)
  red_h_seq <- seq(from=h_baseline*0.5, to=h_baseline*1.5, length.out=7)
  gh_x <- matrix(NA, ncol=length(undersmooth_seq), nrow=length(x_seq))
  gh_x_i<-matrix(NA, ncol=length(undersmooth_seq), nrow=n)
  smoothing_list<-list()
  smoothing_list_i<-list()
  trace_L<-c()
  vtilde<-c()
  bandwidths<-c()
  
  h_cv<-cv_local_linear(kernel_function=kernel_function, x=x, Y=Y, h=red_h_seq)
  #iterating over the undersmoothing vector
  for (f in 1:length(undersmooth_seq)){
    smoothing_matrix <- matrix(NA, nrow=length(x), ncol=length(x))
    smoothing_matrix_grid <- matrix(NA, nrow=length(x_seq), ncol=length(x))
    
    factor_undersmooth<-undersmooth_seq[f]
    h<-n^(factor_undersmooth)*h_cv
    #implementing the estimator
    for (k in 1:length(x_seq)){
      xeval=x_seq[k]
      b<-calc_S1_S2_func(kernel_function=kernel_function, x=x, xeval=xeval, h=h)$b
      l<-c()
      for (i in 1:length(x)){
        l[i]<-b[i]/sum(b)
      }
      smoothing_matrix_grid[k,] <- l
      gh_x[k,f]<-l%*%Y
    }
    
    for (k in 1:length(x)){
      xeval=x[k]
      b<-calc_S1_S2_func(kernel_function=kernel_function, x=x, xeval=xeval, h=h)$b
      l<-c()
      for (i in 1:length(x)){
        l[i]<-b[i]/sum(b)
      }
      smoothing_matrix[k,] <- l
      gh_x_i[k, f]<-l%*%Y
    }
    
    trace_L[f]<-sum(diag(smoothing_matrix))
    vtilde[f]<-sum(diag(crossprod(smoothing_matrix, smoothing_matrix)))
    smoothing_list[[f]]<-smoothing_matrix_grid
    smoothing_list_i[[f]]<-smoothing_matrix
    bandwidths[f]<-h
  }
  return(list(x=x, Y=Y, gh_x=gh_x, gh_x_i=gh_x_i, smoothing_list=smoothing_list,
              smoothing_list_i=smoothing_list_i, trace_L=trace_L, vtilde=vtilde, bandwidths=bandwidths))
}

pointwise_intervals_local_linear<-function(label, kernel_function, n, p1, p2, mu,
                                           var, x_seq, e, r, alpha, undersmooth_seq, 
                                           variance_estimator){
  
  #calculates pointwise confidence intervals for local linear regression
  
  ##arguments: - label: string for the simulation number
  #            - kernel_function: kernel_function
  #            - n: sample size
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - mu: mean of normal density
  #            - var: variance of normal density
  #            - x_seq: sequence of x values 
  #            - e: lower end of uniform distribution
  #            - r: upper end of uniform distribution
  #            - undersmooth_seq: sequence of undersmoothing factors
  #            - variance_estimator: string to denote the type of variance
  #               estimation
  
  ##returns:   - coverage_matrix: matrix with booleans for coverage 
  #            - mse_matrix: matrix with pointwise MSEs
  
  
  g_x<-g_x_calc(x=x_seq, label=label, mu=mu, var=var)
  
  results<-local_linear_estimator(kernel_function=kernel_function, label=label,  
                                    n=n, p1=p1, p2=p2, mu=mu, var=var, x_seq=x_seq, e=e, r=r, 
                                    undersmooth_seq=undersmooth_seq)
  
  coverage_matrix<-matrix(NA, nrow=length(x_seq), ncol=length(undersmooth_seq)*2)
  mse_matrix<-matrix(NA, nrow=length(x_seq), ncol=length(undersmooth_seq))
  
  for (f in 1:length(undersmooth_seq)){
    results_df<-data.frame('gh_x'=results$gh_x[,f], 'g_x'=g_x)
    trace_L<- results$trace_L[f]
    
    l_norm<-sqrt(rowSums(results$smoothing_list[[f]]^2, na.rm = TRUE))
    
    
    #calculating standard errors
    if (variance_estimator=="non-robust"){
      standard_error_not_robust<-sqrt(sum((results$Y-results$gh_x_i[,f])^2)/(n-2*trace_L+results$vtilde[f]))
      se<-standard_error_not_robust*l_norm
     
    }
    
    if (variance_estimator=="robust"){
      se<-error_variance_estimation(x=results$x, Y=results$Y, gh_x_i=results$gh_x_i[,f],
                                    smoothing_matrix=results$smoothing_list_i[[f]], 
                                    smoothing_matrix_grid=results$smoothing_list[[f]], 
                                    kernel_function=kernel_function)
    
    }
    
    bias_correction<-calculate_bias(label=label, x=x_seq, p1=p1, p2=p2, bandwidth=results$bandwidths[f], 
                                    mu=mu, var=var)$bias
    E_gh<-g_x+bias_correction
    
    quantile<-qnorm(1-alpha/2)
    
    CI_lower<-results_df$gh_x-quantile*se
    CI_upper<-results_df$gh_x+quantile*se
    
    #getting coverage and MSE matrices
    
    for (i in 1:length(x_seq)){
      coverage_matrix[i,f] <- (CI_lower[i]<g_x[i]&g_x[i]<CI_upper[i])
      coverage_matrix[i, length(undersmooth_seq)+f] <- (CI_lower[i]<E_gh[i]&E_gh[i]<CI_upper[i])
    }
    
    for (i in 1:length(x_seq)){
      mse_matrix[i,f] <- (results_df$gh_x[i]-results_df$g_x[i])^2
    }
    
  }
  return(list(coverage_matrix=coverage_matrix,
              mse_matrix=mse_matrix, results_df=results_df, 
              CI_upper=CI_upper, CI_lower=CI_lower))
}

local_linear_empirical<-function(x, Y, kernel_function, x_seq){
  red_h_seq=seq(from=5, to=40, length.out=10)
  h_cv<-cv_local_linear(kernel_function=kernel_function, x=x, Y=Y, h=red_h_seq)
  smoothing_matrix_grid<-matrix(NA, nrow=length(x_seq), ncol=length(x))
  gh_x<- c()
  for (k in 1:length(x_seq)){
    xeval=x_seq[k]
    b<-calc_S1_S2_func(kernel_function=kernel_function, x=x, xeval=xeval, h=h_cv)$b
    l<-c()
    for (i in 1:length(x)){
    l[i]<-b[i]/sum(b)
    }
    smoothing_matrix_grid[k,] <- l
    gh_x[k]<-l%*%Y
    }
  
  smoothing_matrix<-matrix(NA, nrow=length(x_seq), ncol=length(x))
  gh_x_i <- c()
  for (k in 1:length(x)){
    xeval=x[k]
    b<-calc_S1_S2_func(kernel_function=kernel_function, x=x, xeval=xeval, h=h_cv)$b
    l<-c()
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix[k,] <- l
    gh_x_i[k]<-l%*%Y
  }
  se<-error_variance_estimation(x=x, Y=Y, gh_x_i=gh_x_i,
                                smoothing_matrix=smoothing_matrix,
                                smoothing_matrix_grid=smoothing_matrix_grid,
                                kernel_function=kernel_function)

  quantile<-qnorm(1-alpha/2)

  CI_lower<-gh_x-quantile*se
  CI_upper<-gh_x+quantile*se
  
  results_df<-data.frame('x_seq'=x_seq, 'gh_x'=gh_x, 'CI_lower'=CI_lower, "CI_upper"=CI_upper)
  
  results_plot<-data.frame('x'=x, 'Y'=Y)
  
  plot<-ggplot(results_df, aes(x=x_seq)) +
    geom_line(aes(y=gh_x), color="black") +
    geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper, fill="Confidence band"),alpha=0.3)+
    #geom_point(data=results_plot, aes(x=x,y=Y), color="black", size=0.05) + 
    theme_bw() + xlab("Potential work experience (years)") + 
    scale_fill_manual("",values="grey12")+
    theme( text = element_text(size = 26), legend.position = "none")+
    ylab("Log hourly wage")
  
  return(list(h_cv=h_cv))
}

local_linear_empirical<-function(x, Y, x_seq){
  
  #performs local linear estimation for the empirical application
  
  # Arguments: x: experience levels
  #            Y: log wage
  #            x_seq: grid points
  # Returns:   plot: empirical plot
  
  h_baseline<-sd(x)*length(x)^(-1/5)
  red_h_seq=seq(from=h_baseline*0.5, to=h_baseline*5, length.out=20)
  h_cv<-cv_local_linear(kernel_function=epa_kernel, x=x, Y=Y, h=red_h_seq)
  h_cv<-length(x)^(-1/50)*h_cv
  smoothing_matrix_grid<-matrix(NA, nrow=length(x_seq), ncol=length(x))
  gh_x<- c()
  for (k in 1:length(x_seq)){
    l<-c()
    xeval=x_seq[k]
    b<-calc_S1_S2_func(kernel_function=epa_kernel, x=x, xeval=xeval, h=h_cv)$b
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix_grid[k,] <- l
    gh_x[k]<-l%*%Y
  }
  
  smoothing_matrix<-matrix(NA, nrow=length(x), ncol=length(x))
  gh_x_i<-c()
  for (k in 1:length(x)){
    l<-c()
    xeval=x[k]
    b<-calc_S1_S2_func(kernel_function=epa_kernel, x=x, xeval=xeval, h=h_cv)$b
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix[k,] <- l
    gh_x_i[k]<-l%*%Y
  }
  #using the leave-one-out errors
  loo_squared_errors<- (((Y-gh_x_i)/(1-diag(smoothing_matrix)))^2)
  h_errors<-sd(x)*length(x)^(-1/5)
  local_constant_results<-local_constant(kernel=epa_kernel, bandwidth=h_errors, x_vector=x, evaluation_vector = x)
  denominator <-local_constant_results$denominator
  numerator <- local_constant_results$numerator
  
  sigma_x_i<-rowSums(numerator*t(replicate(length(x),loo_squared_errors)))/denominator
  v_x<-as.vector(sigma_x_i%*%t(smoothing_matrix_grid^2))
  s_x <- sqrt(v_x)
  
  
  CI_upper<-gh_x-qnorm(0.95)*s_x
  CI_lower<-gh_x+qnorm(0.95)*s_x
  
  results_df<-data.frame('x_seq'=x_seq, 'gh_x'=gh_x, 'CI_lower'=CI_lower, "CI_upper"=CI_upper)
  
  plot<-ggplot(results_df, aes(x=x_seq)) +
    geom_line(aes(y=gh_x), color="black") +
    geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper, fill="Confidence band"),alpha=0.3)+
    #geom_point(data=results_plot, aes(x=x,y=Y), color="black", size=0.05) + 
    theme_bw() + xlab("Potential work experience (years)") + 
    scale_fill_manual("",values="grey12")+
    theme( text = element_text(size = 26), legend.position = "none")+
    ylab("Log hourly wage")
  
  return(plot)
}



plot_bias_local_function<-function(df_plot){
  # Plots the bias for local constant
  # Arguments: df_plot: dataframe with the plot data
  # Returns: Plot_bias: bias plot
plot_bias<-ggplot(df_plot, aes(x=x_seq)) +
  geom_line(aes(y=g_x, colour="True regression function"))+
  
  geom_line(aes(y=E_gh, colour="Expected value of local constant estimator")) +
  scale_colour_manual("",
                      breaks = c("True regression function", "Expected value of local constant estimator"),
                      values = c("red", "blue"))+
  xlab("x")+
  theme_bw()+
  theme(legend.direction="vertical", axis.title.y = element_blank(), text = element_text(size = 20))
  return(plot_bias)
}

plot_mse_function<- function(x_seq, data1, data2, limits, s){
  # Plots MSE for sim 4 and 5
  
  #Arguments: - x_seq: grid points
  #           - data1: local constant data
  #           - data2: local linear data
  #           - limits: plot limits
  #           - s: string for simulation
  # Returns:  saves plot
  for (i in 1:3){
    plot_df_mse<-data.frame('x_seq'=x_seq, 'mse_constant'=100*data1$mse_matrix_replication[,i], 'mse_linear'=100*data2$mse_matrix_replication[,i])
    plot_bias<-ggplot(plot_df_mse, aes(x=x_seq)) +
      
      geom_line(aes(y=mse_constant, colour="Local Constant"))+
      geom_line(aes(y=mse_linear, colour="Local Linear")) +
      scale_colour_manual("",
                          breaks = c("Local Constant", "Local Linear"),
                          values = c("red", "blue"))+
      xlab("x")+
      ylab("100 x MSE")+
      ylim(limits)+
      theme_bw()+
      theme(legend.direction="horizontal", legend.position="bottom", text = element_text(size = 25))
    ggsave(filename=sprintf("plot%s_%i.png", s, i), plot=plot_bias, device="png",width=8, height=6, dpi=500)
  }
}

plot_coverage_function<- function(x_seq, data1, data2, limits, s){
  
  # Plots coverage for sim 4 and 5
  
  #Arguments: - x_seq: grid points
  #           - data1: local constant data
  #           - data2: local linear data
  #           - limits: plot limits
  #           - s: string for simulation
  # Returns:  saves plot
  
  for (i in 1:3){
    plot_df_coverage<-data.frame('x_seq'=x_seq, 'coverage_constant'=data1$coverage_matrix_replication[,i], 'coverage_linear'=data2$coverage_matrix_replication[,i])
    plot_df_coverage<-plot_df_coverage[plot_df_coverage$x_seq>0.14&plot_df_coverage$x_seq<0.86, ]
    plot_bias<-ggplot(plot_df_coverage, aes(x=x_seq)) +
      
      geom_line(aes(y=coverage_constant, colour="Local Constant"))+
      geom_line(aes(y=coverage_linear, colour="Local Linear")) +
      geom_hline(yintercept=0.95, col="green")+
      scale_colour_manual("",
                          breaks = c("Local Constant", "Local Linear"),
                          values = c("red", "blue"))+
      xlab("x")+
      ylab("Coverage")+
      ylim(limits)+
      theme_bw()+
      
      theme(legend.direction="horizontal", legend.position="bottom", text = element_text(size = 25))
    ggsave(filename=sprintf("plotthings%s_%i.png", s, i), plot=plot_bias, device="png",width=8, height=6, dpi=500)
  }
}


plot_function<-function(df_plot){
  #Plots the regression function
  # Arguments: df_plot: dataframe with true conditional mean fct
  #Returns: plot: plot with the conditional mean function
  plot<-ggplot(df_plot, aes(x=x_seq)) +
    geom_line(aes(y=g_x), colour="black")+
    xlab("x")+
    ylab("Regression function")+
    theme_bw()+
    theme(legend.direction="vertical", text = element_text(size = 20))
  return(plot)
}



plot_beta_density<-function(p1, p2){
  #plots the beta density
  #Arguments: p1: parameter 1 of the beta distribution
  #           p2: parameter 2of the beta distribution
  # returns: plot: beta distribution plot
  x_seq=seq(from=0.01, to=0.99, by=0.01)
  dens_beta<-data.frame('x_seq'=x_seq, 'density'=dbeta(x_seq, p1, p2))
  plot<-ggplot(dens_beta, aes(x=x_seq))+
    geom_line(aes(y=density), color="black")+
    xlab("x")+
    ylab("Density")+
    theme_bw()+
    theme(text = element_text(size = 15))
  return(plot)
}


plot_bias_linear_function<-function(df_plot){
  # plots bias of local linear estimator
  # arguments: df_plot: dataframe with expected value and true function
  #returns: bias plot
  plot_bias<-ggplot(df_plot, aes(x=x_seq)) +
    geom_line(aes(y=g_x, colour="True regression function"))+
    geom_line(aes(y=E_gh, colour="Expected value of local linear estimator")) +
    scale_colour_manual("",
                        breaks = c("True regression function", "Expected value of local linear estimator"),
                        values = c("red", "blue"))+
    xlab("x")+
    theme_bw()+
    theme(legend.direction="vertical", axis.title.y = element_blank(), text = element_text(size = 20))
  return(plot_bias)
}


replication_function<-function(reps, estimator, kernel_function, label, n, p1, p2, mu, var, x_seq,
                               e, r, alpha, undersmooth_seq, variance_estimator){
  
  # performs simulation runs
  
  # Arguments: - reps: number of Monte carlo simulations
  #            - kernel_function: kernel function
  #            - label: string for the simulation number
  #            - n: sample size
  #            - p1: parameter 1 of the beta distribution
  #            - p2: parameter 2 of the beta distribution
  #            - x_seq: points of evaluation
  #            - e: lower end of uniform distribution
  #            - r: upper end of uniform distribution
  #            - alpha: critical value
  #            - undersmooth_seq: undersmoothing sequence
  #            - variance estimator: string to denote the type of variance estimator
  
  # Returns:   - coverage_matrix_replication: mean of booleans for
  #                     of the individual runs
  #            - mse_matrix_replication: mean of pointwise squared loss
  
  coverage_matrix_replication<-matrix(NA, nrow=length(x_seq), ncol=length(undersmooth_seq)+1)
  mse_matrix_replication <- matrix(NA, nrow=length(x_seq), ncol=length(undersmooth_seq))
  
  if (estimator=="local_constant"){
  
  replication<-replicate(reps, pointwise_intervals_local_constant(label=label, kernel_function = kernel_function,
                                                     n=n, p1=p1, p2=p2, mu=mu, var=var, x_seq=x_seq, e=e, r=r,
                                                     alpha=alpha, undersmooth_seq=undersmooth_seq, 
                                                     variance_estimator = variance_estimator))
  }
  
  if (estimator=="local_linear"){
    
    replication<-replicate(reps, pointwise_intervals_local_linear(label=label, kernel_function = kernel_function,
                                                           n=n, p1=p1, p2=p2, mu=mu, var=var, x_seq=x_seq, e=e, r=r,
                                                           alpha=alpha, undersmooth_seq=undersmooth_seq, 
                                                           variance_estimator = variance_estimator))
  }
  
  #iterating over the coverage and mse matrices of the individual runs
  # and calculating means
  coverage<-c()
    for(i in 1:length(x_seq)){
      for(j in 1:(length(undersmooth_seq)+1)){
        for (rep in 1:reps){
        coverage[rep]<-replication['coverage_matrix',][[rep]][i,j]
        }
        coverage_matrix_replication[i, j]<-mean(coverage, na.rm=TRUE)
      }
    }
  
  mse<-c()
  for(i in 1:length(x_seq)){
    for(j in 1:(length(undersmooth_seq))){
      for (rep in 1:reps){
        mse[rep]<-replication['mse_matrix',][[rep]][i,j]
      }
      mse_matrix_replication[i, j]<-mean(mse, na.rm=TRUE)
    }
  }
  
  
  return(list(coverage_matrix_replication=coverage_matrix_replication, mse_matrix_replication=mse_matrix_replication))
}



multivariate_function<-function(n, factor_undersmooth, grid, var, cov, df){
  #multivariate local constant estimation
  
  #Arguments: n: sample size
  #           factor_undersmooth: undersmooth_factor
  #           grid: grid of evaluation points
  #          var: variance of normal distribution
  #           cov: covariance of normal distribution
  #           df: degrees of freedom of t distribution
  #Return    coverage: coverage vector
  #           mse: mse vector
data<-rmvnorm(n, mean=c(0,0), sigma=matrix(c(var, cov, cov, var), 2, 2))

g_x_function<-function(x1, x2){
  g<-exp(-(x1^2+x2^2)^2)
  return(g)
}
g_x<-g_x_function(grid[,1], grid[,2])

Y=exp(-(data[,1]^2+data[,2]^2)^2)+rt(n, df=df)#rgamma(n, shape=shape, scale=scale)-shape*scale

#getting the bandwidths

bw1=sd(data[,1])*n^(-1/6)
bw2=sd(data[,2])*n^(-1/6)
bw<-c(bw1, bw2)*n^factor_undersmooth

k<-matrix(NA, nrow=n, ncol=2)
numerator_grid<-matrix(NA, nrow=nrow(grid), ncol=nrow(data))

#iterating over the grid

for (j in 1:nrow(grid)){
 for (i in 1:n){
  for (d in 1:2){
    k[i,d]=epa_kernel(v=(data[i,d]-grid[j,d])/bw[d])
  }
 }
   numerator_grid[j,]<-rowProds(k)
}
denominator_grid<-rowSums(numerator_grid)

smoothing_matrix_grid<-matrix(NA, nrow=nrow(grid), ncol=n)

  for (i in 1:nrow(grid)){
    for (j in 1:nrow(data)){
      smoothing_matrix_grid[i,j]<-numerator_grid[i,j]/denominator_grid[i]
    }
  }

gh_x=rowSums(numerator_grid*t(replicate(nrow(grid),Y)))/denominator_grid

numerator<-matrix(NA, nrow=n, ncol=n)

k<-matrix(NA, nrow=n, ncol=2)

for (j in 1:n){
  for (i in 1:n){
    for (d in 1:2){
      k[i,d]=epa_kernel(v=(data[i,d]-data[j,d])/bw[d])
    }
  }
  numerator[j,]<-rowProds(k)
}
denominator<-rowSums(numerator)

smoothing_matrix<-matrix(NA, nrow=n, ncol=n)

for (i in 1:n){
  for (j in 1:n){
    smoothing_matrix[i,j]<-numerator[i,j]/denominator[i]
  }
}

#from here on out it is the same as done above
gh_x_i=rowSums(numerator*t(replicate(n,Y)))/denominator

trace_L<-sum(diag(smoothing_matrix))
vtilde<-sum(diag(crossprod(smoothing_matrix, smoothing_matrix)))

l_norm<-sqrt(rowSums(smoothing_matrix_grid^2, na.rm = TRUE))

  standard_error_not_robust<-sqrt(sum((Y-gh_x_i)^2)/(n-2*trace_L+vtilde))
  se<-standard_error_not_robust*l_norm
  
  quantile<-qnorm(0.95)
  
  CI_lower<-gh_x-quantile*se
  CI_upper<-gh_x+quantile*se
  coverage<-c()
  for (i in 1:nrow(grid)){
    coverage[i] <- (CI_lower[i]<g_x[i]&g_x[i]<CI_upper[i])
  }
  mse<-c()
  for (i in 1:nrow(grid)){
    mse[i] <- (gh_x[i]-g_x[i])^2
  }
  
  return(list(coverage=coverage, mse=mse))
}


replication_multivariate<-function(reps, n, factor_undersmooth, grid, var, cov, df){
  #replication function for multivariate case
  #arguments: reps: number of repetitions
  #           n: sample size
  #           factor_undersmooth: undersmooth_factor
  #           grid: grid of evaluation points
  #          var: variance of normal distribution
  #           cov: covariance of normal distribution
  #           df: degrees of freedom of t distribution
replication<-replicate(reps, multivariate_function(n=n, grid=grid, factor_undersmooth = factor_undersmooth,
                                                   var=var, cov=cov, df=df))
coverage<-colMeans(do.call(rbind, replication['coverage', ]), na.rm=TRUE)
mse<-colMeans(do.call(rbind, replication['mse', ]), na.rm=TRUE)

return(list(coverage=coverage, mse=mse))
}

