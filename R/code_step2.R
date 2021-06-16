#' STEP 2 IN PNC ANALYSIS - CONVERT TO NETWORK
#' Obtain data in format time x 264 for each subject and convert to network.
#' @param images list of images containing region level summaries of nifti data, one object per subject in the list.
#' @param lambda vector of length n (number of subjects) of lambda used in graphical lasso per subject .
#' @return List of adjacency matrices for each subject.
#' @export
network.estimate <- function(images, lambda = FALSE){
  # library(huge)
  file.names <- names(images)
  list_adj_mat <- list()
  list_lambda = NULL
  for(i in 1:length(images)){
    subj_data <- t(images[[file.names[i]]])
    trans_data = huge::huge.npn(subj_data)
    mb_fit = huge::huge(trans_data)
    mb_graph = huge::huge.select(mb_fit)
    adj_mat = as.matrix(mb_graph$refit)
    list_adj_mat[[i]] <- adj_mat # create list of adjacency matrices
    list_lambda = c(list_lambda, mb_graph$opt.lambda)
  }
  names(list_adj_mat) <- file.names
  if(!lambda){
    return(list_adj_mat)
  }else{
    return(list(adj_mat = list_adj_mat, lambda = list_lambda))
  }
}
