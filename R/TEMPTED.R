#' @title Decomposition of temporal tensor
#' @description This is the main function of tempted.
#' @param datlist A length n list of matrices.
#' Each matrix represents a subject,
#' with columns representing samples from this subject,
#' the first row representing the sampling time points,
#' and the following rows representing the feature values.
#' @param r Number of components to decompose into, i.e. rank of the CP type decomposition.
#' Default is set to 3.
#' @param smooth Smoothing parameter for RKHS norm.
#' Larger means smoother temporal loading functions. Default is set to be 1e-8.
#' Value can be adjusted depending on the dataset by checking the smoothness of the estimated temporal loading function in plot.
#' @param interval The range of time points to ran the decomposition for.
#' Default is set to be the range of all observed time points.
#' User can set it to be a shorter interval than the observed range.
#' @param resolution Number of time points to evaluate the value of the temporal loading function.
#' Default is set to 101. It does not affect the subject or feature loadings.
#' @param maxiter Maximum number of iteration. Default is 20.
#' @param epsilon Convergence criteria for difference between iterations. Default is 1e-4.
#' @return The estimations of the loadings.
#' \describe{
#'   \item{A_hat}{Subject loading, a subject by r matrix.}
#'   \item{B_hat}{Feature loading, a feature by r matrix.}
#'   \item{Phi_hat}{Temporal loading function, a resolution by r matrix.}
#'   \item{time_Phi}{The time points where the temporal loading function is evaluated.}
#'   \item{Lambda}{Eigen value, a length r vector.}
#'   \item{r_square}{Variance explained by each component. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using individual components.}
#'   \item{accum_r_square}{Variance explained by the first few components accumulated. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using the first few components.}
#' }
#' @examples
#' # Take a subset of the samples so the example runs faster
#'
#' # Here we are taking samples from the odd months
#' sub_sample <- rownames(meta_table)[(meta_table$day_of_life%/%12)%%2==1]
#' count_table_sub <- count_table[sub_sample,]
#' processed_table_sub <- processed_table[sub_sample,]
#' meta_table_sub <- meta_table[sub_sample,]
#'
#' # for count data from longitudinal microbiome studies
#'
#' datlist <- format_tempted(count_table_sub,
#'                           meta_table_sub$day_of_life,
#'                           meta_table_sub$studyid,
#'                           pseudo=0.5,
#'                           transform="clr")
#'
#' mean_svd <- svd_centralize(datlist, r=1)
#'
#' res_tempted <- tempted(mean_svd$datlist, r=2, smooth=1e-5)
#'
#' # for preprocessed data that do not need to be transformed
#' \donttest{
#' datlist <- format_tempted(processed_table_sub,
#'                           meta_table_sub$day_of_life,
#'                           meta_table_sub$studyid,
#'                           pseudo=NULL,
#'                           transform="none")
#'
#' mean_svd <- svd_centralize(datlist, r=1)
#'
#' res_tempted <- tempted(mean_svd$datlist, r=2, smooth=1e-5)
#' }
#' # plot the temporal loading
#'
#' plot_time_loading(res_tempted, r=2)
#' @export
#' @md
tempted <- function(datlist, r = 3, smooth=1e-6,
                    interval = NULL, resolution = 101,
                    maxiter=20, epsilon=1e-4){
  # a function to iteratively update phi
  freg_rkhs <- function(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth=1e-8){
    A <- Kmat
    for (i in 1:length(Ly)){
      A[ind_vec==i,] <- A[ind_vec==i,]*a_hat[i]^2
    }
    cvec <- unlist(Ly)

    A_temp <- A + smooth*diag(ncol(A))
    beta <- solve(A_temp)%*%cvec

    phi_est <- Kmat_output %*% beta
    return(phi_est)
  }
  n <- length(datlist)
  p <- nrow(datlist[[1]])-1

  Lambda <- rep(0, r)
  A <- matrix(0, n, r)
  B <- matrix(0, p, r)
  Phi <- matrix(0, resolution, r)
  PCname <- paste0('PC', 1:r)
  colnames(A) <- PCname
  colnames(B) <- PCname
  colnames(Phi) <- PCname
  rownames(A) <- names(datlist)
  rownames(B) <- rownames(datlist[[1]])[-1]


  # Calculate range.
  timestamps_all <- do.call(c,lapply(datlist, FUN=function(u){u[1,]}))

  timestamps_all <- sort(unique(timestamps_all))
  if (is.null(interval)){
    interval <- c(timestamps_all[1], timestamps_all[length(timestamps_all)])
  }

  # rescale the time to 0-1.
  input_time_range <- c(timestamps_all[1], timestamps_all[length(timestamps_all)])
  for (i in 1:n){
    datlist[[i]][1,] <- (datlist[[i]][1,] - input_time_range[1]) / (input_time_range[2] - input_time_range[1])
  }
  interval <- (interval - input_time_range[1]) / (input_time_range[2] - input_time_range[1])

  res <- NULL
  Lambda <- rep(0, r)
  X <- NULL
  y0 <- NULL
  Rsq <- accumRsq <- rep(0, r)

  ti <- vector(mode='list', length=n)
  for (i in 1:n){
    temp <- 1 + round((resolution-1) * (datlist[[i]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] <- 0
    ti[[i]] <- temp
  }

  tipos <- vector(mode='list', length=n)
  for (i in 1:n){
    keep <- ti[[i]]>0
    tipos[[i]] <- keep
    y0 <- c(y0, as.vector(t(datlist[[i]][2:(p+1),keep])))
  }

  Lt <- list()
  ind_vec <- NULL
  for (i in 1:n){
    Lt <- c(Lt, list(datlist[[i]][1,]))
    ind_vec <- c(ind_vec, rep(i,length(Lt[[i]])))
  }

  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)

  # calculate rank-1 component sequentially.
  for (s in 1:r){
    # Step 1: initialization.
    message(sprintf("Calculate the %dth Component", s))

    # intialization of b
    data_unfold = NULL
    y <- NULL
    for (i in 1:n){
      data_unfold = cbind(data_unfold, datlist[[i]][2:(p+1),])
      y <- c(y, as.vector(t(datlist[[i]][2:(p+1),tipos[[i]]])))
    }
    b.initials <- svd(data_unfold, nu=r, nv=r)$u
    b_hat <- b.initials[,1]
    # initialization of a
    a_hat <- rep(1,n)/sqrt(n)

    # iteratively update a,b,phi
    t <- 0
    dif <- 1
    while(t<=maxiter & dif>epsilon){
      # update phi:
      Ly <- list()
      for (i in 1:n){
        Ly <- c(Ly, list(a_hat[i]*as.numeric(b_hat%*%datlist[[i]][2:(p+1),])))
      }
      phi_hat <- freg_rkhs(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
      phi_hat <- phi_hat / sqrt(sum(phi_hat^2))

      # update a:
      a_tilde <- rep(0,n)
      for (i in 1:n){
        t_temp <- tipos[[i]]
        a_tilde[i] <- b_hat %*% datlist[[i]][2:(p+1),t_temp] %*% phi_hat[ti[[i]][t_temp]]
        a_tilde[i] <- a_tilde[i] / sum((phi_hat[ti[[i]][t_temp]])^2)
      }
      a.new <- a_tilde / sqrt(sum(a_tilde^2))
      dif <- sum((a_hat - a.new)^2)
      a_hat <- a.new

      # update b:
      temp_num <- matrix(0,p,n)
      temp_denom <- rep(0,n)
      for (i in 1:n){
        t_temp <- tipos[[i]]
        temp_num[,i] <- datlist[[i]][2:(p+1),t_temp] %*% phi_hat[ti[[i]][t_temp]]
        temp_denom[i] <-sum((phi_hat[ti[[i]][t_temp]])^2)
      }
      b_tilde <- as.numeric(temp_num%*%a_hat) / as.numeric(temp_denom%*%(a_hat^2))
      b.new <- b_tilde / sqrt(sum(b_tilde^2))
      dif <- max(dif, sum((b_hat - b.new)^2))
      b_hat <- b.new

      t <- t+1
    }

    # calculate lambda
    x <- NULL
    for (i in 1:n){
      t_temp <- ti[[i]]
      t_temp <- t_temp[t_temp>0]
      x <- c(x,as.vector(t(a_hat[i]*b_hat%o%phi_hat[t_temp])))
    }
    X <- cbind(X, x)
    lm_fit <- lm(y~x-1)
    lambda <- as.numeric(lm_fit$coefficients)
    A[,s] <- a_hat
    B[,s] <- b_hat
    Phi[,s] <- t(phi_hat)
    Lambda[s] <- lambda
    Rsq[s] <- summary(lm_fit)$r.squared
    accumRsq[s] <- summary(lm(y0~X-1))$r.squared

    # update datlist
    for (i in 1:n){
      temp <- tipos[[i]]
      datlist[[i]][2:(p+1),which(temp)] <- datlist[[i]][2:(p+1),which(temp)] -
        Lambda[s] * A[i,s] * (B[,s] %*% t(Phi[ti[[i]][temp],s]))
    }
    message(paste0("Convergence reached at dif=", dif, ', iter=', t))
  }
  lm_fit <- lm(y0~X-1)
  Lambda <- as.numeric(lm_fit$coefficients)

  # revise the sign of Lambda
  for (r in 1:length(Lambda)){
    if (Lambda[r]<0){
      Lambda[r] <- -Lambda[r]
      A[,r] <- -A[,r]
    }
  }
  # revise the signs to make sure summation of phi is nonnegative
  sgn.phi <- sign(colSums(Phi))
  sgn.phi[sgn.phi==0] <- 1
  for (r in 1:ncol(Phi)){
    Phi[,r] <- sgn.phi[r]*Phi[,r]
    A[,r] <- sgn.phi[r]*A[,r]
  }
  # revise the signs to make sure summation of B is nonnegative
  sgn.B <- sign(colSums(B))
  sgn.B[sgn.B==0] <- 1
  for (r in 1:ncol(Phi)){
    B[,r] <- sgn.B[r]*B[,r]
    A[,r] <- sgn.B[r]*A[,r]
  }
  time_return <- seq(interval[1],interval[2],length.out = resolution)
  time_return <- time_return * (input_time_range[2] - input_time_range[1]) + input_time_range[1]
  results <- list("A_hat" = A, "B_hat" = B,
                 "Phi_hat" = Phi, "time_Phi" = time_return,
                 "Lambda" = Lambda, "r_square" = Rsq, "accum_r_square" = accumRsq)
  return(results)
}


#' @title Remove the mean structure of the temporal tensor
#' @description This function first average the feature value of all time points for each subject to form a subject by feature matrix.
#' Next, it performs a singular value decomposition of this matrix and construct the matrix's rank-r approximation.
#' Then, it subtracts this rank-r subject by feature matrix from the temporal tensor.
#' @param datlist A length n list of matrices.
#' Each matrix represents a subject,
#' with columns representing samples from this subject,
#' the first row representing the sampling time points,
#' and the following rows representing the feature values.
#' @param r The number of ranks in the mean structure. Default is 1.
#' @return A list of results.
#' \describe{
#'   \item{datlist}{The new temporal tensor after mean structure is removed.}
#'   \item{A_tilde}{The subject singular vector of the mean structure, a subject by r matrix.}
#'   \item{B_tilde}{The feature singular vector of the mean structure, a feature by r matrix.}
#'   \item{lambda_tilde}{The singular value of the mean structure, a length r vector.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R, Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749. \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @seealso Examples can be found in \code{\link{tempted}}.
#' @export
#' @md
svd_centralize <- function(datlist, r = 1){
  n <- length(datlist)
  p <- nrow(datlist[[1]])-1
  mean_hat <- matrix(0,n,p)
  for (i in 1:length(datlist)){
    mean_hat[i,] <- apply(datlist[[i]][-1,], 1, mean)
  }
  mean_hat_svd <- svd(mean_hat, nu=r, nv=r)
  mean_hat_svd1 <- mean_hat_svd$u %*% t(mean_hat_svd$v * mean_hat_svd$d[1:r])
  mf.new <- datlist
  for (i in 1:length(datlist)){
    mf.new[[i]][-1,] <- datlist[[i]][-1,] - mean_hat_svd1[i,]
  }
  results <- list("datlist" = mf.new, "A_tilde" = mean_hat_svd$u,
                  "B_tilde" = mean_hat_svd$v, "lambda_tilde" = mean_hat_svd$d[1:r])
  return(results)
}


#' @title Format data table into the input of tempted
#' @description This function applies a variety of transformations to the read counts and
#' format the sample by feature table and meta data into a data list
#' that can be used as the input of \code{\link{tempted}} and \code{\link{svd_centralize}}.
#' For data that are not read counts, or data that are not microbiome data,
#' the user can apply their desired transformation to the data before formatting into list.
#' @param featuretable A sample by feature matrix.
#' @param timepoint The time stamp of each sample, matched with the rows of \code{featuretable}.
#' @param subjectID The subject ID of each sample, matched with the rows of \code{featuretable}.
#' @param threshold A threshold for feature filtering for microbiome data.
#' Features with zero value percentage > threshold will be excluded. Default is 0.95.
#' @param pseudo A small number to add to all the counts before
#' normalizing into proportions and log transformation.
#' Default is 1/2 of the smallest non-zero value that is specific for each sample.
#' This pseudo count is added for \code{transform=c("logcomp", "clr", "logit")}.
#' @param transform The transformation applied to the data.
#' \code{"logcomp"} for log of compositions.
#' \code{"comp"} for compositions.
#' \code{"ast"} for arcsine squared transformation.
#' \code{"clr"} for central log ratio transformation.
#' \code{"lfb"} for log 2 fold change over baseline (first time point) transformation.
#' \code{"logit"} for logit transformation.
#' \code{"none"} for no transformation.
#' Default \code{transform="clr"} is recommended for microbiome data.
#' For data that are already transformed, use \code{transform="none"}.
#' @return A length n list of matrices as the input of \code{\link{tempted}} and \code{\link{svd_centralize}}.  Each matrix represents a subject, with columns representing samples from this subject, the first row representing the sampling time points, and the following rows representing the feature values.
#' @seealso Examples can be found in \code{\link{tempted}}.
#' @export
#' @md
format_tempted <- function(featuretable, timepoint, subjectID,
                           threshold=0.95, pseudo=NULL, transform="clr"){
  ntm <- which(table(subjectID)==1)
  if(length(ntm)>0)
    stop(paste('Please remove these subjects with only one time point:',
               paste(names(ntm), collapse=', ')))
  if (length(subjectID)!=nrow(featuretable))
    stop('length of subjectID does not match featuretable!')
  if (length(timepoint)!=nrow(featuretable))
    stop('length of timepoint does not match featuretable!')
  # get pseudo count
  if (is.null(pseudo) & (transform %in% c("clr", "logcomp", "logit"))){
    pseudo <- apply(featuretable, 1, function(x){
      min(x[x!=0])/2
    })
  }
  # keep taxon that has non-zeros in >1-threshold samples
  featuretable <- featuretable[,colMeans(featuretable==0)<=threshold]
  if(transform=='logcomp' | transform=="lfb"){
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable)))
  }else if(transform=='comp'){
    featuretable <- featuretable
    featuretable <- t(featuretable/rowSums(featuretable))
  }else if(transform=='ast'){
    featuretable <- featuretable
    featuretable <- t(asin(sqrt(featuretable/rowSums(featuretable))))
  }else if(transform=='clr'){
    featuretable <- featuretable+pseudo
    featuretable <- log(featuretable/rowSums(featuretable))
    featuretable <- t(featuretable-rowMeans(featuretable))
  }else if(transform=="lfb"){
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable), 2))
  }else if(transform=='logit'){
    featuretable <- featuretable+pseudo
    featuretable <- t(featuretable/rowSums(featuretable))
    featuretable <- log(featuretable/(1-featuretable))
  }else if(transform=='none'){
    featuretable <- t(featuretable)
  }else{
    message('Input transformation method is wrong! logcomp is applied instead')
    featuretable <- featuretable+pseudo
    featuretable <- t(log(featuretable/rowSums(featuretable)))
  }
  featuretable <- rbind(timepoint, featuretable)
  rownames(featuretable)[1] <- 'timepoint'
  subID <- unique(subjectID)
  nsub <- length(subID)

  # construct list of data matrices, each element representing one subject
  datlist <- vector("list", length = nsub)
  names(datlist) <- subID

  # Each slice represents an individual (unequal sized matrix).
  for (i in 1:nsub){
    datlist[[i]] <- featuretable[, subjectID==subID[i]]
    datlist[[i]] <- datlist[[i]][,order(datlist[[i]][1,])]
    datlist[[i]] <- datlist[[i]][,!duplicated(datlist[[i]][1,])]
    if(transform=="lfb"){
      datlist[[i]] <- datlist[[i]][,-1, drop=FALSE] - datlist[[i]][,1]
    }
  }
  return(datlist)
}




#' @title Caculate the Bernoulli kernel
#' @description This function is used to calculate the kernel matrix for the
#' RKHS regression that iteratively updates the temporal loading function.
#' @param x,y Two values between which the Bernoulli kernel is calculated.
#' @return The calculated kernel between \code{x} and \code{y}.
#' @references
#' Han, R., Shi, P. and Zhang, A.R. (2023) \emph{Guaranteed functional tensor singular value decomposition}. Journal of the American Statistical Association, pp.1-13. doi: 10.1080/01621459.2022.2153689.
#' @md
bernoulli_kernel <- function(x, y){
  k1_x <- x-0.5
  k1_y <- y-0.5
  k2_x <- 0.5*(k1_x^2-1/12)
  k2_y <- 0.5*(k1_y^2-1/12)
  xy <- abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4_xy <- 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern_xy <- k1_x %*% t(k1_y) + k2_x %*% t(k2_y) - k4_xy + 1
  return(kern_xy)
}




#' @title Calculate the de-noised temporal tensor
#' @description This function constructs a de-noised version of the temporal tensor
#' using the low-rank components obtained by \code{\link{svd_centralize}} \code{\link{tempted}} and uses the loadings to
#' @param res_tempted Output of tempted
#' @param mean_svd Output of svd_centralize
#' @return The de-noised functional tensor
#' @export
#' @md
tdenoise <- function(res_tempted, mean_svd=NULL){
  n <- nrow(res_tempted$A_hat)
  p <- nrow(res_tempted$B_hat)
  resol <- nrow(res_tempted$Phi_hat)
  tensor_est <- array(0,dim=c(n,p,resol))
  if (!is.null(mean_svd))
    tensor_est <- (mean_svd$A_tilde %*% t(mean_svd$B_tilde * mean_svd$lambda_tilde)) %o%
    rep(1, resol)
  for (i in 1:ncol(res_tempted$A_hat)){
    tensor_est <- tensor_est+res_tempted$A_hat[,i]%o%res_tempted$B_hat[,i]%o%res_tempted$Phi_hat[,i]*res_tempted$Lambda[i]
  }
  dimnames(tensor_est)[[3]] <- res_tempted$time_Phi
  return(tensor_est)
}


#' @title Estimate subject loading of testing data
#' @description This function estimates the subject loading of the testing data
#' based on feature and temporal loading from training data,
#' so that both the testing data and training data have the same dimensionality reduction.
#' @param datlist Testing data formatted into datlist in the same fashion as the training data.
#' The same transformation needs to be used for both training and testing data.
#' @param res_tempted Result from \code{\link{tempted}} ran on the training data.
#' @param mean_svd Result from \code{\link{svd_centralize}} ran on the training data.
#' @return estimated subject loading of testing data
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R, Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749. \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @examples
#' # Take a subset of the samples so the example runs faster
#'
#' # Here we are taking samples from the odd months
#' sub_sample <- rownames(meta_table)[(meta_table$day_of_life%/%12)%%2==1]
#' count_table_sub <- count_table[sub_sample,]
#' processed_table_sub <- processed_table[sub_sample,]
#' meta_table_sub <- meta_table[sub_sample,]
#'
#' # split the example data into training and testing
#'
#' id_test <- meta_table_sub$studyid=="2"
#'
#' count_train <- count_table_sub[!id_test,]
#' meta_train <- meta_table_sub[!id_test,]
#'
#' count_test <- count_table_sub[id_test,]
#' meta_test <- meta_table_sub[id_test,]
#'
#' # run tempted on training data
#'
#' datlist_train <- format_tempted(count_train,
#'                                 meta_train$day_of_life,
#'                                 meta_train$studyid,
#'                                 threshold=0.95,
#'                                 pseudo=0.5,
#'                                 transform="clr")
#'
#' mean_svd_train <- svd_centralize(datlist_train, r=1)
#'
#' res_tempted_train <- tempted(mean_svd_train$datlist,
#' r=2, smooth=1e-5)
#'
#' # get the overlapping features
#'
#' count_test <- count_test[,rownames(datlist_train[[1]])[-1]]
#'
#' datlist_test <- format_tempted(count_test,
#'                                meta_test$day_of_life,
#'                                meta_test$studyid,
#'                                threshold=1,
#'                                pseudo=0.5,
#'                                transform="clr")
#'
#' # estimate the subject loading of the testing subject
#'
#' sub_test <- est_test_subject(datlist_test, res_tempted_train, mean_svd_train)
#'
#' # train logistic regression classifier on training subjects
#'
#' metauni <- unique(meta_table_sub[,c("studyid", "delivery")])
#' rownames(metauni) <- metauni$studyid
#' Atrain <- as.data.frame(res_tempted_train$A_hat)
#' Atrain$delivery <- metauni[rownames(Atrain),"delivery"]=="Cesarean"
#' glm_train <- glm(delivery ~ PC1+PC2,
#'                  data=Atrain, family=binomial(link="logit"))
#' summary(glm_train)
#'
#' # predict the label of testing subject, whose true label is "Cesarean"
#'
#' predict(glm_train, newdata=as.data.frame(sub_test), type="response")
#'
#' @export
#' @md
est_test_subject <- function(datlist, res_tempted, mean_svd=NULL){
  B <- res_tempted$B_hat
  Phi <- res_tempted$Phi_hat
  Lambda <- res_tempted$Lambda
  time_return <- res_tempted$time_Phi
  n <- length(datlist)
  p <- nrow(B)
  r <- ncol(B)
  resolution <- length(time_return)
  A_test <- matrix(0,n,r)
  y <- NULL
  ti <- vector(mode = "list", length = n)
  # get the coordinate of observed time points in the returned time grid
  for (i in 1:n){
    ti[[i]] <- sapply(datlist[[i]][1,], function(x){which.min(abs(x-time_return))})
    y <- c(y, as.numeric(t(datlist[[i]][-1,ti[[i]]>0])))
  }
  mf.new <- datlist
  if(!is.null(mean_svd)){
    mean_hat <- matrix(0,n,p)
    for (i in 1:length(datlist)){
      mean_hat[i,] <- apply(datlist[[i]][-1,], 1, mean)
    }
    mean_hat_svd1 <- mean_hat%*%tcrossprod(mean_svd$B_tilde)
    mf.new <- datlist
    for (i in 1:length(datlist)){
      mf.new[[i]][-1,] <- datlist[[i]][-1,] - mean_hat_svd1[i,]
    }
  }

  for (s in 1:r){
    for (i in 1:n){
      t_temp <- ti[[i]]>0
      A_test[i,s] <- B[,s] %*% mf.new[[i]][2:(p+1),t_temp] %*% Phi[ti[[i]][t_temp],s]
      A_test[i,s] <- A_test[i,s] / sum((Phi[ti[[i]][t_temp],s])^2) / Lambda[s]
      mf.new[[i]][2:(p+1),t_temp] <- mf.new[[i]][2:(p+1),t_temp] -
        Lambda[s] * A_test[i,s] * (B[,s] %*% t(Phi[ti[[i]][t_temp],s]))
    }
  }
  rownames(A_test) <- names(mf.new)
  colnames(A_test) <- paste0('PC', 1:r)
  return(A_test)
}


#' @title Aggregate features using feature loadings
#' @description This function aggregate the features into "meta features" by
#' calculating a weighted summation of the features using feature loading of each component as weights.
#' It can also aggregate features by using the combination of multiple components by ranking the features
#' by a linear combination of feature loadings from multiple components.
#' @param res_tempted Output of \code{\link{tempted}}.
#' @param mean_svd Output of \code{\link{svd_centralize}}.
#' @param datlist Output of \code{\link{format_tempted}}, the original temporal tensor that will be aggregated.
#' @param pct The percent of features to aggregate,
#' features ranked by absolute value of the feature loading of each component.
#' Default is 1, which means 100% of features are aggregated.
#' Setting \code{pct=0.01} means top 1% of features is aggregated,
#' where features are ranked in absolute value of feature loading of each component.
#' @param contrast A matrix choosing how components are combined,
#' each column is a contrast of length r and used to calculate the linear combination of
#' the feature loadings of r components.
#' @return A list of results.
#' \describe{
#'   \item{metafeature_aggregate}{The meta feature obtained by aggregating the observed temporal tensor. It is a data.frame with four columns: "value" for the meta feature values, "subID" for the subject ID, "timepoint" for the time points, and "PC" indicating which component was used to construct the meta feature.}
#'   \item{metafeature_aggregate_est}{The meta feature obtained by aggregating the denoised temporal tensor. It has the same structure as \code{metafeature_aggregate}.}
#'   \item{contrast}{The contrast used to linearly combine the components from input.}
#'   \item{toppct}{A matrix of TRUE/FALSE indicating which features are aggregated in each component and contrast.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R, Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749. \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @examples
#' # Take a subset of the samples so the example runs faster
#'
#' # Here we are taking samples from the odd months
#' sub_sample <- rownames(meta_table)[(meta_table$day_of_life%/%12)%%2==1]
#' count_table_sub <- count_table[sub_sample,]
#' processed_table_sub <- processed_table[sub_sample,]
#' meta_table_sub <- meta_table[sub_sample,]
#'
#' datlist <- format_tempted(count_table_sub,
#'                           meta_table_sub$day_of_life,
#'                           meta_table_sub$studyid,
#'                           pseudo=0.5,
#'                           transform="clr")
#'
#' mean_svd <- svd_centralize(datlist, r=1)
#'
#' res_tempted <- tempted(mean_svd$datlist, r=2, smooth=1e-5)
#'
#' contrast <- matrix(c(1/2,1), 2, 1)
#'
#' res_aggregate <- aggregate_feature(res_tempted,
#'                                    mean_svd,
#'                                    datlist,
#'                                    pct=1,
#'                                    contrast=contrast)
#'
#' # plot the aggregated features
#'
#' \donttest{
#' group <- unique(meta_table[, c("studyid", "delivery")])
#'
#' plot_metafeature(res_aggregate$metafeature_aggregate, group, bws=30)
#' }
#' @export
#' @md
aggregate_feature <- function(res_tempted, mean_svd=NULL, datlist,
                               pct=1,
                              contrast=NULL){
  B_data <- as.data.frame(res_tempted$B_hat)
  r <- ncol(B_data)
  if(!is.null(contrast)){
    contrast_data <- res_tempted$B_hat%*%contrast
    colnames(contrast_data) <- paste0('Contrast', 1:ncol(contrast))
    B_data <- cbind(B_data, contrast_data)
  }
  toppct <- apply(abs(B_data), 2, function(x){x>=quantile(x, 1-pct)})
  datalist_agg <- sapply(datlist, function(x){t(B_data*toppct)%*%x[-1,]}, simplify=F)
  metafeature_aggregate <- NULL
  for (i in 1:length(datalist_agg)){
    tmp <- data.frame(value=as.vector(datalist_agg[[i]]),
                      subID=names(datalist_agg)[i],
                      timepoint=rep(datlist[[i]][1,], each=ncol(B_data)),
                      PC=rep(rownames(datalist_agg[[i]]), ncol(datalist_agg[[i]])))

    metafeature_aggregate <- rbind(metafeature_aggregate, tmp)
  }
  metafeature_aggregate <- metafeature_aggregate[,c("value", "subID", "timepoint", "PC")]

  # estimated
  tensor_est <- tdenoise(res_tempted, mean_svd)
  tensor_est_agg <- apply(tensor_est, c(1,3), function(x){(t(B_data*toppct)%*%x)})
  metafeature_aggregate_est <- NULL
  for (i in 1:r){
    tmp <- data.frame(value=as.vector(tensor_est_agg[i,,]),
                      subID=rep(dimnames(tensor_est_agg)[[2]], dim(tensor_est_agg)[3]),
                      timepoint=rep(res_tempted$time_Phi, each=dim(tensor_est_agg)[2]),
                      PC=colnames(B_data)[i])
    metafeature_aggregate_est <- rbind(metafeature_aggregate_est,tmp)
  }
  metafeature_aggregate_est <- metafeature_aggregate_est[,c("value", "subID", "timepoint", "PC")]
  metafeature_aggregate_est$type <- 'estimated'

  return(list(metafeature_aggregate=metafeature_aggregate,
              metafeature_aggregate_est=metafeature_aggregate_est,
              contrast=contrast,
              toppct=toppct))
}


#' @title Take log ratio of the abundance of top features over bottom features
#' @description Top and bottom ranking features are picked based on feature loadings (and their contrasts).
#' The log ratio abundance of the top ranking features over the bottom ranking features is produced as the main result.
#' This function and its result is designed for longitudinal microbiome data,
#' and may not be meaningful for other type of temporal data.
#' @param res_tempted Output of \code{\link{tempted}}.
#' @param datlist Output of \code{format_tempted(, transform="none")}, the temporal tensor that include the raw read counts.
#' @param absolute \code{absolute = TRUE} means features are ranked by the absolute value of feature loadings,
#' and the top \code{pct} percent of features are picked.
#' \code{absolute = FALSE} means features are ranked by the original value of feature loadings,
#' and the top and bottom \code{pct} percent of features are picked.
#' Then ratio is taken as the abundance of the features with positive loading
#' over the abundance of the features with negative loading.
#' @param pct The percent of features to sum up. Default is 0.05, i.e. 5%.
#' @param contrast A matrix choosing how components are combined,
#' each column is a contrast of length r and used to calculate the linear combination of
#' the feature loadings of r components.
#' @return A list of results:
#' \describe{
#'   \item{metafeature_ratio}{The log ratio abundance of the top over bottom ranking features. It is a data.frame with five columns: "value" for the log ratio values, "subID" for the subject ID, and "timepoint" for the time points, and "PC" indicating which component was used to construct the meta feature.}
#'   \item{contrast}{The contrast used to linearly combine the components from input.}
#'   \item{toppct}{A matrix of TRUE/FALSE indicating which features are ranked top in each component (and contrast) and used as the numerator of the log ratio.}
#'   \item{bottompct}{A matrix of TRUE/FALSE indicating which features are ranked bottom in each component (and contrast) and used as the denominator of the log ratio.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R, Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749. \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @examples
#' # Take a subset of the samples so the example runs faster
#'
#' # Here we are taking samples from the odd months
#' sub_sample <- rownames(meta_table)[(meta_table$day_of_life%/%12)%%2==1]
#' count_table_sub <- count_table[sub_sample,]
#' processed_table_sub <- processed_table[sub_sample,]
#' meta_table_sub <- meta_table[sub_sample,]
#'
#' datlist <- format_tempted(count_table_sub,
#'                           meta_table_sub$day_of_life,
#'                           meta_table_sub$studyid,
#'                           pseudo=0.5,
#'                           transform="clr")
#'
#' mean_svd <- svd_centralize(datlist, r=1)
#'
#' res_tempted <- tempted(mean_svd$datlist, r=2, smooth=1e-5)
#'
#' datalist_raw <- format_tempted(count_table_sub, meta_table_sub$day_of_life, meta_table_sub$studyid,
#' transform="none")
#'
#' contrast <- cbind(c(1,1), c(1,-1))
#'
#' res_ratio <- ratio_feature(res_tempted, datalist_raw, pct=0.1,
#' absolute=FALSE, contrast=contrast)
#'
#' group <- unique(meta_table[, c("studyid", "delivery")])
#'
#' # plot the log ratios
#' \donttest{
#' plot_metafeature(res_ratio$metafeature_ratio, group, bws=30)
#' }
#' @export
#' @md
ratio_feature <- function(res_tempted, datlist,
                              pct=0.05, absolute=FALSE, contrast=NULL){
  B_data <- as.data.frame(res_tempted$B_hat)
  if (!is.null(contrast)){
    contrast_data <- res_tempted$B_hat%*%contrast
    colnames(contrast_data) <- paste0('Contrast', 1:ncol(contrast))
    B_data <- cbind(B_data, contrast_data)
  }
  if(!absolute){
    toppct <- apply(B_data, 2, function(x){x>quantile(x, 1-pct) & x>0})
    bottompct <- apply(-B_data, 2, function(x){x>quantile(x, 1-pct) & x>0})
  }else{
    toppct <- apply(B_data, 2, function(x){abs(x)>quantile(abs(x), 1-pct) & x>0})
    bottompct <- apply(B_data, 2, function(x){abs(x)>quantile(abs(x), 1-pct) & x<0})
  }
  pseudo <- min(sapply(datlist, function(x){
    y<-x[-1,]
    return(min(y[y!=0]))
    }))/2
  datlist_ratio <- sapply(datlist,
                          function(x){
                            tt <- (t(toppct)%*%x[-1,])
                            bb <- (t(bottompct)%*%x[-1,])
                            return(log((tt+pseudo)/(bb+pseudo)))},
                          simplify=F)
  metafeature_ratio <- NULL
  for (i in 1:length(datlist_ratio)){
    tmp <- data.frame(value=as.vector(datlist_ratio[[i]]),
                      subID=names(datlist_ratio)[i],
                      timepoint=rep(datlist[[i]][1,], each=ncol(B_data)),
                      PC=rep(rownames(datlist_ratio[[i]]), ncol(datlist_ratio[[i]])))

    metafeature_ratio <- rbind(metafeature_ratio, tmp)
  }
  metafeature_ratio <- metafeature_ratio[,c("value", "subID", "timepoint", "PC")]
  return(list(metafeature_ratio=metafeature_ratio,
              contrast=contrast,
              toppct=toppct, bottompct=bottompct))
}

#' @title Run all major functions of tempted
#' @description This function wraps functions \code{\link{format_tempted}}, \code{\link{svd_centralize}}, \code{\link{tempted}},
#' \code{\link{ratio_feature}}, \
#' and \code{\link{aggregate_feature}}.
#' @param r Number of components to decompose into, i.e. rank of the CP type decomposition.
#' Default is set to 3.
#' @param smooth Smoothing parameter for RKHS norm.
#' Larger means smoother temporal loading functions. Default is set to be 1e-8.
#' Value can be adjusted depending on the dataset by checking the smoothness of the estimated temporal loading function in plot.
#' @param featuretable A sample by feature matrix. It is an input for \code{\link{format_tempted}}.
#' @param timepoint The time stamp of each sample, matched with the rows of \code{featuretable}. It is an input for \code{\link{format_tempted}}.
#' @param subjectID The subject ID of each sample, matched with the rows of \code{featuretable}. It is an input for \code{\link{format_tempted}}.
#' @param threshold A threshold for feature filtering for microbiome data.
#' Features with zero value percentage >= threshold will be excluded. Default is 0.95.
#' It is an input for \code{\link{format_tempted}}.
#' @param pseudo A small number to add to all the counts before
#' normalizing into proportions and log transformation.
#' Default is 1/2 of the smallest non-zero value that is specific for each sample.
#' This pseudo count is added for \code{transform=c("logcomp", "clr", "logit")}.
#' It is an input for \code{\link{format_tempted}}.
#' @param transform The transformation applied to the data.
#' \code{"logcomp"} for log of compositions.
#' \code{"comp"} for compositions.
#' \code{"ast"} for arcsine squared transformation.
#' \code{"clr"} for central log ratio transformation.
#' \code{"logit"} for logit transformation.
#' \code{"none"} for no transformation.
#' Default \code{transform="clr"} is recommended for microbiome data.
#' For data that are already transformed, use \code{transform="none"}.
#' It is an input for \code{\link{format_tempted}}.
#' @param r Number of components to decompose into, i.e. rank of the CP type decomposition.
#' Default is set to 3.
#' It is an input for \code{\link{tempted}}.
#' @param smooth Smoothing parameter for RKHS norm.
#' Larger means smoother temporal loading functions. Default is set to be 1e-8.
#' Value can be adjusted depending on the dataset by checking the smoothness of the estimated temporal loading function in plot.
#' It is an input for \code{\link{tempted}}.
#' @param interval The range of time points to ran the decomposition for.
#' Default is set to be the range of all observed time points.
#' User can set it to be a shorter interval than the observed range.
#' It is an input for \code{\link{tempted}}.
#' @param resolution Number of time points to evaluate the value of the temporal loading function.
#' Default is set to 101. It does not affect the subject or feature loadings. It is an input for \code{\link{tempted}}.
#' @param maxiter Maximum number of iteration. Default is 20. It is an input for \code{\link{tempted}}.
#' @param epsilon Convergence criteria for difference between iterations. Default is 1e-4. It is an input for \code{\link{tempted}}.
#' @param r_svd The number of ranks in the mean structure. Default is 1. It is an input for \code{\link{svd_centralize}}.
#' @param pct_ratio The percent of features to sum up. Default is 0.05, i.e. 5%.
#' It is an input for \code{\link{ratio_feature}}.
#' @param absolute \code{absolute = TRUE} means features are ranked by the absolute value of feature loadings,
#' and the top \code{pct_ratio} percent of features are picked.
#' \code{absolute = FALSE} means features are ranked by the original value of feature loadings,
#' and the top and bottom \code{pct_ratio} percent of features are picked.
#' Then ratio is taken as the abundance of the features with positive loading
#' over the abundance of the features with negative loading.
#' It is an input for \code{\link{ratio_feature}}.
#' @param pct_aggregate The percent of features to aggregate,
#' features ranked by absolute value of the feature loading of each component.
#' Default is 1, which means 100% of features are aggregated.
#' Setting \code{pct_aggregate=0.01} means top 1% of features is aggregated,
#' where features are ranked in absolute value of feature loading of each component.
#' It is an input for \code{\link{aggregate_feature}}.
#' @param contrast A matrix choosing how components are combined,
#' each column is a contrast of length r and used to calculate the linear combination of
#' the feature loadings of r components.
#' It is an input for \code{\link{ratio_feature}} and It is an input for \code{\link{aggregate_feature}}.
#' @param do_ratio Whether to calculate the log ratio of features.
#' @return A list including all the input and output of functions \code{\link{format_tempted}}, \code{\link{svd_centralize}}, \code{\link{tempted}},
#' \code{\link{ratio_feature}}, and \code{\link{aggregate_feature}}.
#' \describe{
#'   \item{input}{All the input options of function \code{\link{tempted_all}}.}
#'   \item{datalist_raw}{Output of \code{\link{format_tempted}} with option \code{transform="none"}.}
#'   \item{datlist}{Output of \code{\link{format_tempted}}.}
#'   \item{mean_svd}{Output of \code{\link{svd_centralize}}.}
#'   \item{A_hat}{Subject loading, a subject by r matrix.}
#'   \item{B_hat}{Feature loading, a feature by r matrix.}
#'   \item{Phi_hat}{Temporal loading function, a resolution by r matrix.}
#'   \item{time_Phi}{The time points where the temporal loading function is evaluated.}
#'   \item{Lambda}{Eigen value, a length r vector.}
#'   \item{r_square}{Variance explained by each component. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using individual components.}
#'   \item{accum_r_square}{Variance explained by the first few components accumulated. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using the first few components.}
#'   \item{metafeature_ratio}{The log ratio abundance of the top over bottom ranking features. It is a data.frame with five columns: "value" for the log ratio values, "subID" for the subject ID, and "timepoint" for the time points, and "PC" indicating which component was used to construct the meta feature.}
#'   \item{toppct_ratio}{A matrix of TRUE/FALSE indicating which features are ranked top in each component (and contrast) and used as the numerator of the log ratio.}
#'   \item{bottompct_ratio}{A matrix of TRUE/FALSE indicating which features are ranked bottom in each component (and contrast) and used as the denominator of the log ratio.}
#'   \item{metafeature_aggregate}{The meta feature obtained by aggregating the observed temporal tensor. It is a data.frame with four columns: "value" for the meta feature values, "subID" for the subject ID, "timepoint" for the time points, and "PC" indicating which component was used to construct the meta feature.}
#'   \item{toppct_aggregate}{A matrix of TRUE/FALSE indicating which features are aggregated in each component and contrast.}
#'   \item{contrast}{The contrast used to linearly combine the components from input.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R, Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749. \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @examples
#' # Take a subset of the samples so the example runs faster
#'
#' # Here we are taking samples from the odd months
#' sub_sample <- rownames(meta_table)[(meta_table$day_of_life%/%12)%%2==1]
#' count_table_sub <- count_table[sub_sample,]
#' processed_table_sub <- processed_table[sub_sample,]
#' meta_table_sub <- meta_table[sub_sample,]
#'
#' # for preprocessed data that do not need to be transformed
#'
#' \donttest{
#' res.processed <- tempted_all(processed_table_sub,
#'                              meta_table_sub$day_of_life,
#'                             meta_table_sub$studyid,
#'                              threshold=1,
#'                              transform="none",
#'                              r=2,
#'                              smooth=1e-5,
#'                              do_ratio=FALSE)
#'
#' # for count data that will have pseudo added and clr transformed
#'
#' res.count <- tempted_all(count_table_sub,
#'                          meta_table_sub$day_of_life,
#'                          meta_table_sub$studyid,
#'                          threshold=0.95,
#'                          transform="clr",
#'                          pseudo=0.5,
#'                          r=2,
#'                          smooth=1e-5,
#'                          pct_ratio=0.1,
#'                          pct_aggregate=1)
#' }
#' # for proportional data that will have pseudo added and clr transformed
#'
#' res.proportion <- tempted_all(count_table_sub/rowSums(count_table_sub),
#'                               meta_table_sub$day_of_life,
#'                               meta_table_sub$studyid,
#'                               threshold=0.95,
#'                               transform="clr",
#'                               pseudo=NULL,
#'                               r=2,
#'                               smooth=1e-5,
#'                               pct_ratio=0.1,
#'                               pct_aggregate=1)
#'
#' # plot the temporal loading and subject trajectories grouped by delivery mode
#'
#' plot_time_loading(res.proportion, r=2)
#'
#' group <- unique(meta_table[,c("studyid", "delivery")])
#'
#' # plot the aggregated features
#' \donttest{
#' plot_metafeature(res.proportion$metafeature_aggregate, group, bws=30)
#' }
#' @export
#' @md
tempted_all <- function(featuretable, timepoint, subjectID,
                        threshold=0.95, pseudo=NULL, transform="clr",
                        r = 3, smooth=1e-6,
                        interval = NULL, resolution = 51,
                        maxiter=20, epsilon=1e-4, r_svd=1,
                        do_ratio=TRUE, pct_ratio=0.05, absolute=FALSE,
                        pct_aggregate=1, contrast=NULL){
  datlist <- format_tempted(featuretable=featuretable, timepoint=timepoint, subjectID=subjectID,
                            threshold=threshold, pseudo=pseudo, transform=transform)
  datalist_raw <- format_tempted(featuretable=featuretable, timepoint=timepoint, subjectID=subjectID,
                                threshold=threshold, pseudo=pseudo, transform="none")
  mean_svd <- svd_centralize(datlist, r_svd)
  res_tempted <- tempted(datlist=mean_svd$datlist, r = r, smooth=smooth,
                         interval = interval, resolution = resolution,
                         maxiter=maxiter, epsilon=epsilon)
  if (do_ratio)
    res_ratio <- ratio_feature(res_tempted=res_tempted, datlist=datalist_raw,
                             pct=pct_ratio, absolute=absolute, contrast=contrast)
  res_aggfeat <- aggregate_feature(res_tempted=res_tempted, mean_svd=mean_svd, datlist=datlist,
                                   pct=pct_aggregate, contrast=contrast)
  res_all <- list(input=as.list(match.call()))
  res_all$datalist_raw <- datalist_raw
  res_all$datlist <- datlist
  res_all <- append(res_all, res_tempted)
  if (do_ratio){
    res_all$metafeature_ratio <- res_ratio$metafeature_ratio
    res_all$toppct_ratio <- res_ratio$toppct
    res_all$bottompct_ratio <- res_ratio$toppct
  }
  res_all$metafeature_aggregate <- res_aggfeat$metafeature_aggregate
  res_all$toppct_aggregate <- res_aggfeat$toppct
  res_all$contrast <- contrast
  return(res_all)
}



#' @title Plot nonparametric smoothed mean and error bands of features versus time
#' @description This is a handy function to plot the smoothed mean and error bands for multiple features.
#' @param feature_mat A sample by feature matrix. Each feature will be plotted separately as a facet.
#' The features can be original features, meta features, log ratios, or any variables of interest.
#' @param time_vec A vector of time points matched to the rows of \code{feature_mat}.
#' @param group_vec A vector of factor variable indicating the group membership
#' of samples matched to the rows of \code{feature_mat}.
#' @param coverage The coverage rate for the error band. Default is 0.95.
#' @param bws The smoothness parameter for the smoothing lines and error bands.
#' A larger value means a smoother line.
#' Default is NULL and calculated by function \code{np::npreg()}.
#' @param nrow The number of rows to plot the features used in function \code{ggplot2::facet_wrap()}.
#' @return A ggplot2 object.
#' @examples
#' # plot the summary of selected features
#'
#' feat.names <- c("OTU4447072", "OTU4467447")
#'
#' proportion_table <- count_table/rowSums(count_table)
#'
#' plot_feature_summary(proportion_table[,feat.names],
#'                      meta_table$day_of_life,
#'                      meta_table$delivery,
#'                      bws=30)
#' @export
#' @md
plot_feature_summary <- function(feature_mat, time_vec, group_vec,
                                 coverage=0.95, bws=NULL, nrow=1){
  nfeature <- ncol(feature_mat)
  if(!is(group_vec, 'factor')) group_vec <- as.factor(group_vec)
  group_level <- levels(group_vec)
  time_all <- NULL
  mean_all <- NULL
  merr_all <- NULL
  feature_all <- NULL
  group_all <- NULL
  if(is.null(colnames(feature_mat))){stop('feature_mat needs to have column names!')}
  CI_length <- -qnorm((1-coverage)/2)
  for (jj in 1:nfeature){
    for (ii in 1:length(group_level)){
      ind <- group_vec==group_level[ii]
      if(is.null(bws)){
        model_np <- npreg(feature_mat[ind,jj]~time_vec[ind],
                          regtyle="ll", bwmethod="cv.aic")
      }else{
        model_np <- npreg(feature_mat[ind,jj]~time_vec[ind], bws=bws,
                          regtyle="ll", bwmethod="cv.aic")
      }
      time_eval <- as.vector(t(model_np$eval))
      mean_eval <- model_np$mean[order(time_eval)]
      merr_eval <- model_np$merr[order(time_eval)]
      time_eval <- sort(time_eval)

      time_all <- c(time_all, time_eval)
      mean_all <- c(mean_all, mean_eval)
      merr_all <- c(merr_all, merr_eval)
      feature_all <- c(feature_all,
                       rep(colnames(feature_mat)[jj], length(time_eval)))
      group_all <- c(group_all,
                     rep(group_level[ii], length(time_eval)))
    }
  }
  group_all <- factor(group_all, levels=group_level)
  tab_summary <- data.frame(time_all=time_all, mean_all=mean_all, merr_all=merr_all,
                            group_all=group_all, feature_all=feature_all)
  .data <- NULL
  p_summary <- ggplot(data=tab_summary,
                      aes(x=.data$time_all, y=.data$mean_all, group=.data$group_all, color=.data$group_all)) +
    geom_line() +
    geom_ribbon(aes(ymin=.data$mean_all-CI_length*.data$merr_all, ymax=.data$mean_all+CI_length*.data$merr_all,
                    color=.data$group_all, fill=.data$group_all), linetype=2, alpha=0.3) +
    ylab(paste0('mean +/- ', round(CI_length,2), '*se')) + facet_wrap(.~.data$feature_all, scales="free", nrow=nrow)
  return(p_summary)
}



#' @title Plot nonparametric smoothed mesan and error bands of meta features versus time
#' @description This function plot the smoothed mean and error band of meta features
#' grouped by a factor variable provided by the user.
#' @param metafeature \code{metafeature_ratio} from the output of \code{\link{ratio_feature}} and
#' \code{\link{tempted_all}},
#' \code{metafeature_aggregate} from the output of
#' \code{\link{ratio_feature}} and \code{\link{tempted_all}},
#' or \code{metafeature_aggregate_est} from the output of \code{\link{ratio_feature}}.
#' @param group A subject by 2 data.frame with the first column for subject ID and second column for group membership.
#' @param coverage The coverage rate for the error band. Default is 0.95.
#' @param bws The smoothness parameter for the smoothing lines and error bands.
#' A larger value means a smoother line.
#' Default is NULL and calculated by function \code{np::npreg()}.
#' @param nrow The number of rows to plot the features used in function \code{ggplot2::facet_wrap()}.
#' @return A ggplot2 object.
#' @seealso Examples can be found in \code{\link{tempted_all}}, \code{\link{ratio_feature}} and \code{\link{aggregate_feature}}.
#' @export
#' @md
plot_metafeature <- function(metafeature, group,
                             coverage=0.95, bws=NULL, nrow=1){
  colnames(group) <- c("subID", "group")
  tab_feat_ratio <- merge(metafeature, group, by="subID")
  ## summed up, by mean and sd
  reshape_feat_ratio <- reshape(tab_feat_ratio,
                                idvar=c("subID","timepoint") ,
                                v.names=c("value"), timevar="PC",
                                direction="wide")
  CC <- grep('value',colnames(reshape_feat_ratio))
  colnames(reshape_feat_ratio)[CC] <- substr(colnames(reshape_feat_ratio)[CC], start=7, stop=100)
  feature_mat_ratio <- reshape_feat_ratio[,CC]
  colnames(feature_mat_ratio)
  time_vec_ratio <- reshape_feat_ratio$timepoint
  group_vec_ratio <- factor(reshape_feat_ratio$group)
  p_feat_summary <- plot_feature_summary(feature_mat_ratio,
                                               time_vec_ratio,
                                               group_vec_ratio, bws=bws, nrow=nrow)
  return(p_feat_summary)
}


#' @title Plot the temporal loading functions
#' @description This function uses \code{ggplot2::geom_line()} in ggplot2 to plot the temporal loading functions from \code{\link{tempted}}.
#' @param res Output of function \code{\link{tempted}}.
#' @param r The number of components to plot. By default all the components estimated by \code{\link{tempted}} will be plotted.
#' @param ... Arguments to put in \code{ggplot2::geom_line(aes(...))}.
#' @return An ggplot2 object.
#' @seealso Examples can be found in \code{\link{tempted_all}} and \code{\link{tempted}}.
#' @export
#' @md
plot_time_loading <- function(res, r=NULL, ...){
  Phi_data <- res$Phi_hat
  if(is.null(r)) r <- ncol(Phi_data)
  Phi_data <- Phi_data[,1:r]
  ntime <- nrow(Phi_data)
  Phi_data <- data.frame(timepoint=res$time_Phi, value=as.vector(Phi_data),
                         component=as.factor(rep(1:r,each=ntime)))
  .data <- NULL
  ptime <- ggplot(data=Phi_data, aes(x=.data$timepoint, y=.data$value, color=.data$component)) + geom_line(aes(...))
  return(ptime)
}



#' @title Reconstruct tensor from low dimensional components
#' @description This function reconstructs the temporal tensor from the low dimensional components extracted by \code{\link{tempted}} and \code{\link{svd_centralize}}. 
#' @param res_tempted Output of function \code{\link{tempted}}.
#' @param res_svd Output of function \code{\link{svd_centralize}} if mean subtraction \code{\link{svd_centralize}} was performed before \code{\link{tempted}}. If mean subtraction was not performed, specify \code{res_svd=NULL} and provide \code{datlist}.
#' @param datlist Output of function \code{\link{format_tempted}} if mean subtraction was not performed and \code{res_svd=NULL}, or left as \code{NULL} if mean subtraction was performed and \code{res_svd} is provided.
#' @param r_reconstruct The number of components from TEMPTED to be used for the tensor reconstruction.
#' @return datlist_est The reconstructed tensor stored in a datlist format same as the output of \code{\link{format_tempted}}.
#' @examples
#' # Take a subset of the samples so the example runs faster
#'
#' # Here we are taking samples from the odd months
#' sub_sample <- rownames(meta_table)[(meta_table$day_of_life%/%12)%%2==1]
#' count_table_sub <- count_table[sub_sample,]
#' processed_table_sub <- processed_table[sub_sample,]
#' meta_table_sub <- meta_table[sub_sample,]
#' # reconstruct with mean subtraction
#' datlist <- format_tempted(processed_table_sub,
#'                           meta_table_sub$day_of_life,
#'                           meta_table_sub$studyid,
#'                           pseudo=NULL,
#'                           transform="none")
#'
#' mean_svd <- svd_centralize(datlist, r=1)
#'
#' res_tempted <- tempted(mean_svd$datlist, r=2, smooth=1e-5)
#' 
#' datlist_est <- reconstruct(res_tempted, mean_svd, datlist=NULL, r_reconstruct=2)
#' vec_est <- unlist(sapply(datlist_est, function(x){x[-1,]}))
#' vec_obs <- unlist(sapply(datlist, function(x){x[-1,]}))
#' R2 <- 1-sum((vec_est-vec_obs)^2)/sum(vec_obs^2)
#' R2
#'
#' # reconstruct without mean subtraction
#' datlist <- format_tempted(processed_table_sub,
#'                           meta_table_sub$day_of_life,
#'                           meta_table_sub$studyid,
#'                           pseudo=NULL,
#'                           transform="none")
#'
#' res_tempted <- tempted(datlist, r=2, smooth=1e-5)
#' 
#' datlist_est <- reconstruct(res_tempted, res_svd=NULL, datlist=datlist, r_reconstruct=2)
#' vec_est <- unlist(sapply(datlist_est, function(x){x[-1,]}))
#' vec_obs <- unlist(sapply(datlist, function(x){x[-1,]}))
#' R2 <- 1-sum((vec_est-vec_obs)^2)/sum(vec_obs^2)
#' R2
#' 
#' @export
#' @md
reconstruct <- function(res_tempted, res_svd=NULL, datlist=NULL, r_reconstruct = NULL){
  n <- max(length(res_svd$datlist), length(datlist))
  if (n==0) {
    stop("One of res_svd or datlist needs to be feeded.")
  }
  datlist_est <- vector(mode='list', length=n)
  
  r <- length(res_tempted$Lambda)
  if (r < r_reconstruct){
    r_reconstruct <- min(r, r_reconstruct)
    print("r_reconstruct is larger than the fitted rank. Used the available ranks instead.")
  }
  
  rec1 <- lapply(seq(1, r_reconstruct), function(s) {
    res_tempted$Lambda[s] * res_tempted$A_hat[, s] %o% res_tempted$B_hat[, s] %o% res_tempted$Phi[, s]
  })
  Yhat <- Reduce("+", rec1)
  
  if (!is.null(res_svd)){
    rec2 <- lapply(seq(1, length(res_svd$lambda_tilde)), function(s) {
      res_svd$lambda_tilde[s] * res_svd$A_tilde[, s] %o% res_svd$B_tilde[, s]
    })
    Ymean <- Reduce("+", rec2)
    for (ii in 1:n){
      datlist_est[[ii]] <- res_svd$datlist[[ii]]
      ti <- 1 + round((length(res_tempted$time_Phi)-1) * (res_svd$datlist[[ii]][1,] - res_tempted$time_Phi[1]) / (res_tempted$time_Phi[length(res_tempted$time_Phi)] - res_tempted$time_Phi[1]))
      datlist_est[[ii]][-1,] <- Yhat[ii,,ti] + Ymean[ii,]
    }
  }else{
    for (ii in 1:n){
      datlist_est[[ii]] <- datlist[[ii]]
      ti <- 1 + round((length(res_tempted$time_Phi)-1) * (datlist[[ii]][1,] - res_tempted$time_Phi[1]) / (res_tempted$time_Phi[length(res_tempted$time_Phi)] - res_tempted$time_Phi[1]))
      datlist_est[[ii]][-1,] <- Yhat[ii,,ti]
    }
  }
  return(datlist_est)
}




#' Meta data table from the ECAM data
#' @format A data.frame with rows representing samples and matching with data.frame \code{count_table} and \code{processed_table}
#' and three columns:
#' \describe{
#'     \item{studyid}{character denoting the subject ID of the infants.}
#'     \item{delivery}{character denoting the delivery mode of the infants.}
#'     \item{day_of_life}{character denoting the age of infants measured in days when microbiome sample was taken.}
#' }
#' @references Bokulich, Nicholas A., et al. "Antibiotics, birth mode, and diet shape microbiome maturation during early life." Science translational medicine 8.343 (2016): 343ra82-343ra82.
#' @md
"meta_table"


#' OTU read count table from the ECAM data
#' @format A data.frame with rows representing samples and matching with data.frame \code{meta_table}
#' and columns representing microbial features (i.e. OTUs). Each entry is a read count.
#' @references Bokulich, Nicholas A., et al. "Antibiotics, birth mode, and diet shape microbiome maturation during early life." Science translational medicine 8.343 (2016): 343ra82-343ra82.
#' @md
"count_table"


#' Central-log-ratio (clr) transformed OTU table from the ECAM data
#' @format A data.frame with rows representing samples and matching with data.frame \code{meta_table}
#' and columns representing microbial features (i.e. OTUs).
#' Entries do not need to be transformed, and will be directly used by \code{\link{tempted}}.
#' This data.frame is used to illustrate how \code{\link{tempted}} can be used for
#' general form of multivariate longitudinal data already preprocessed by user.
#' @references Bokulich, Nicholas A., et al. "Antibiotics, birth mode, and diet shape microbiome maturation during early life." Science translational medicine 8.343 (2016): 343ra82-343ra82.
#' @md
"processed_table"

