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
#' Step 4: Run each of the n_log time series through an AR(1) filter with first order correlation phi
#' 
#' Step 5: Generate observed correlation matrix and vector of first order correlations the time series
#' @example /inst/examples/spat_temp_ts_examples.R
#' @seealso \code{\link{spat_cor_mat}}
spat_temp_ts <- function(n_iter, n_loc, site_sd, spat_sd, phi, cor_mat = NULL, log = TRUE) {
    
    ### correct site_sd for AR(1) structure
    site_sd2 <- sqrt(site_sd^2 * (1 - phi^2))
    
    ### create spatial covariance matrix if it is not provided
    if (is.null(cor_mat)) {
        ### generate the spatial correlation matrix among sites start by generating the
        ### spatial correlation by distance
        cor_mat <- spat_cor_mat(n_loc, spat_sd = spat_sd)
    } else {
        cor_mat <- cor_mat
    }
    
    ### create the spatial covariance matrix from the correlation matrix ###
    vcov_mat <- diag(rep(site_sd2, n_loc)) %*% cor_mat %*% diag(rep(site_sd2, n_loc))
    
    ### create the empty time series
    d <- ts(matrix(0, ncol = n_loc, nrow = n_iter))
    
    ### create the spatially correlated ts of errors
    e <- ts(rmvnorm(n_iter, sigma = vcov_mat))
    ### use the error vector to populate the ts
    for (i in 2:n_iter) {
        d[i, ] <- phi * d[i - 1, ] + e[i, ]
    }
    
    ### return spatial correlations and pacf are approximately correct
    fin_pacf_vec <- apply(d, 2, function(x) pacf(x, plot = F, lag_max = 1)$acf)[1, 
        ]
    fin_cor_mat <- cor(d)
    
    ### provide temporal SD for each location ###
    fin_sd_vec <- apply(d, 2, sd)
    if (log == TRUE) {
        fin_ts <- ts(t(apply(d, 1, "-", colMeans(d) + 0.5 * fin_sd_vec^2)))
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
#' @examples 
#' par(mfrow = c(3,2),mai= c(0.75,0.75,0.75,0.75),mgp = c(1.5,0.5,0.1))
#' 
#'matplot(spat_cor_mat(n_loc=10, spat_scale = 1),type= "l",ylab= "correlation",xlab= "site",main= "Cauchy, scale= 1")
#'matplot(spat_cor_mat(n_loc=10, spat_sd = 1),type= "l",ylab= "correlation",xlab= "site",main= "Gaussian, sd= 1")
#' matplot(spat_cor_mat(n_loc=10, spat_scale = 2),type= "l",ylab= "correlation",xlab= "site",main= "Cauchy, scale= 2")
#' matplot(spat_cor_mat(n_loc=10, spat_sd = 2),type= "l",ylab= "correlation",xlab= "site",main= "Gaussian, sd= 2")
#' matplot(spat_cor_mat(n_loc=10, spat_scale = 10),type= "l",ylab= "correlation",xlab= "site",main= "Cauchy, scale= 10")
#' matplot(spat_cor_mat(n_loc=10, spat_sd = 10),type= "l",ylab= "correlation",xlab= "site",main= "Gaussian, sd= 10")
spat_cor_mat <- function(n_loc, spat_sd = 1, spat_scale = NULL, sumto1 = FALSE) {
    if (is.null(spat_scale)) {
        dn1 <- dnorm(c(1:(2 * n_loc - 1)), mean = n_loc, sd = spat_sd)
        dn1 <- dn1/dn1[n_loc]
    } else {
        dn1 <- dcauchy(c(1:(2 * n_loc - 1)), location = n_loc, scale = spat_scale)
        dn1 <- dn1/dn1[n_loc]
    }
    
    ### now generate the correlation matrix ###
    cor_mat <- matrix(dn1[n_loc:(2 * n_loc - 1)])
    for (i in 1:(n_loc - 1)) {
        cor_mat <- cbind(cor_mat, dn1[(n_loc:(2 * n_loc - 1)) - (i)])
    }
    
    if (sumto1 == TRUE) {
        cor_mat <- apply(cor_mat, 2, "/", colSums(cor_mat))
    }
    return(cor_mat)
}

#' Generate a series of random stray matrices using a mean stray probability matrix
#' @param stray_mat a mean stray matrix.
#' @param n_iter number of observations.
#' @param scale higher values provide lower variability in random samples.
#' @description Generate a series of random stray matrices by sampling from the Dirichlet distribution with a mean stray probability matrix and scale parameter
ran_stray_prob <- function(stray_mat,n_iter,scale){
  Crand <- array(NA, dim = c(dim(stray_mat), n_iter))
  for (i in 1:n_iter) {
    Crand[, , i] <- apply(matrix(stray_mat * scale, ncol = ncol(stray_mat)), 2, 
                          rdirichlet, n = 1)
  }
  Crand
}


