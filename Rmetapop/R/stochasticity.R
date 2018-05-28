#' Generate a spatially and temporally correlated time series 
#' @param n_iter number of observations in the multivariate time series.
#' @param n_loc number of locatons. 
#' @param site_sd univariate standard devation of each time series.
#' @param spat_sd scale of spatial correlation. 
#' @param phi AR(1) model coefficient.
#' @param cor_mat user supplied correlation matrix (if unspecified will be generated from the spat_sd parameter).
#' @description Simulates a multivarate normal, AR(1) time series that can be used for random errors.  Each time point includes vector of correlated deviations, which exhibit either Gaussian spatial correlations OR a specified correlation matrix.   If log = TRUE, the sd is adjusted such that the final expected sd is equal to that specified.
#' 
#' Step 1: Generate or supply a spatial correlation matrix (cor_mat)
#' 
#' Step 2: Generate a variance-covariance matrix using site_sd and cor_mat
#' 
#' Step 3: Create a time series of multivariate normal vectors (each vector is of length n_loc) with n_iter vector observations in the time series. 
#' 
#' Step 4: Run each of the n_loc time series through an AR(1) filter with first order correlation phi
#' 
#' Step 5: Generate observed correlation matrix and vector of first order correlations the time series
#' @example /inst/examples/spat_temp_ts_examples.R
#' @seealso \code{\link{spat_cor_mat}}
spat_temp_ts <- function(n_iter, n_loc, site_sd, spat_sd, phi, cor_mat = NULL, log = TRUE,offset=0.1,circular= TRUE) {
  
  ### correct site_sd for AR(1) structure
  site_sd2 <- sqrt(site_sd^2 * (1 - phi^2))
  
  ### create spatial covariance matrix if it is not provided
  if (is.null(cor_mat)) {
    ### generate the spatial correlation matrix among sites start by generating the
    ### spatial correlation by distance
    cor_mat <- spat_cor_mat(n_loc, spat_sd = spat_sd,offset= offset,circular= TRUE)
  } else {
    cor_mat <- cor_mat
  }
  
  ### create the spatial covariance matrix from the correlation matrix ###
  # vcov_mat <- diag(rep(site_sd2, n_loc)) %*% cor_mat %*% diag(rep(site_sd2, n_loc))
  vcov_mat <- cor_mat
  ### create the empty time series
  d <- ts(matrix(0, ncol = n_loc, nrow = n_iter))
  
  ### create the spatially correlated ts of errors
  e <- ts(mvrnorm(n_iter,d[1,],Sigma= vcov_mat))
  ### use the error vector to populate the ts
  for (i in 2:n_iter) {
    d[i, ] <- (phi * d[i - 1, ] + e[i, ]*site_sd2)
  }
  
  ### return spatial correlations and pacf are approximately correct
  fin_pacf_vec <- apply(d, 2, function(x) pacf(x, plot = FALSE, na.action= na.pass)$acf)[1, ]
  fin_cor_mat <- cor(d)
  
  ### provide temporal SD for each location ###
  fin_sd_vec <- apply(d, 2, sd)
  if (log == TRUE) {
    if(n_loc==1){
      fin_ts <- ts(apply(d, 1, "-", colMeans(d) + 0.5 * fin_sd_vec^2))
    } else {
      fin_ts <-ts(t(apply(d, 1, "-", colMeans(d) + 0.5 * fin_sd_vec^2)))
    }
  } else {
    fin_ts <- ts(t(apply(d, 1, "-", colMeans(d))))
  }
  return(list(pacf = fin_pacf_vec, cor_mat = fin_cor_mat, sd = fin_sd_vec, ts = fin_ts))
}

#' Generate a spatially correlated matrix
#' @param n_loc number of locatons. 
#' @param spat_sd scale of spatial correlation for the Gaussian case. 
#' @param spat_sd scale of spatial correlation for the Cauchy case. 
#' @param sumto1 if TRUE makes the columns sum to 1.
#' @description Generate a spatially correlated matrix where correlations decrease in space via a Gaussian process (with spat_sd) OR a Cauchy process (with spat_scale)
#' @example /inst/examples/stray_mat_example.R
spat_cor_mat <- function(n_loc, spat_sd = 1, spat_scale = NULL, 
                         sumto1 = FALSE,distance= NULL,offset= 0.1,circular= TRUE) {
  
  cor_mat <- matrix(0,ncol=n_loc,nrow= n_loc)
  
  if(n_loc>2){
    if(is.null(distance)){
      distance = dist(1:n_loc,upper= TRUE,diag= TRUE)
      distance = as.matrix(distance/max(distance+1))
    }
  } else {
    distance= matrix(c(1,0,0,1),ncol=2)
  }
  if (is.null(spat_scale)) {
    if(circular ==TRUE){
      for(i in 1:(n_loc)){
        for(j in 1:(n_loc)){
          cor_mat[i,j] <- exp(-2*(sin(pi*abs(distance[i,j]))^2)/spat_sd^2)
          cor_mat[i,i] <- 1+offset
        }
      }
    } else {
      for(i in 1:(n_loc)){
        for(j in 1:(n_loc)){
          cor_mat[i,j] <- exp(-(distance[i,j])^2/(2*spat_sd^2))
          cor_mat[i,i] <- 1+offset
        }
      }
    }
  } else {
    if(circular ==TRUE){
      for(i in 1:(n_loc)){
        for(j in 1:(n_loc)){
          cor_mat[i,j] <- exp(-2*(sin(pi*abs(distance[i,j]))^2)/spat_scale^2)
          cor_mat[i,i] <- 1+offset
        }
      }
    } else {
      for(i in 1:(n_loc)){
        for(j in 1:(n_loc)){
          cor_mat[i,j] <- (1/(1+(distance[i,j])^2)^(spat_scale))
          cor_mat[i,i] <- 1+offset
        }
      }
    }
  }
  if (sumto1 == TRUE) {
    cor_mat <- apply(cor_mat,2,function(x) x/sum(x))
  } else {
    cor_mat <- cor_mat/max(cor_mat)
  }
  return(cor_mat)
}

#' Generate a series of random stray matrices using a mean stray probability matrix
#' @param stray_mat a mean stray matrix.
#' @param n_iter number of observations.
#' @param scale higher values provide lower variability in random samples.
#' @description Generate a series of random stray matrices by sampling from the Dirichlet distribution with a mean stray probability matrix and scale parameter
#' @example  /inst/examples/stochastic_stray_matrix_example.R
ran_stray_prob <- function(stray_mat,n_iter,scale){
  Crand <- array(NA, dim = c(dim(stray_mat), n_iter))
  for (i in 1:n_iter) {
    Crand[, , i] <- apply(matrix(stray_mat * scale, ncol = ncol(stray_mat)), 2, 
                          rdirichlet, n = 1)
  }
  Crand
}

#' Generate a series of correlated, random survival matrices
#' @param stray_mat a mean stray matrix.
#' @param mean mean survival rate
#' @param corr the overdispersion parameter
#' @description Generate a series of random survival elements by sampling from the beta-binomial
#' @example  /inst/examples/stochastic_stray_matrix_example.R
#'                                                      
ran_surv_prob <- function(mean,corr,N=10000,reps= 1){
  if(!is.null(dim(mean))){
    M <- apply(mean,1:length(dim(mean)),function(x) rbetabinom(reps,size= N,prob= mean, rho= corr)/N)
  }
  else{
    M <- rbetabinom(reps,size= N,prob= mean, rho= corr)/N
  }
  return(M)
}


