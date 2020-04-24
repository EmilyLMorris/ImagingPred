#' STEP 2 IN PNC ANALYSIS - CONVERT TO NETWORK
#' Obtain data in format time x 264 for each subject and convert to network.
#' @param files list of files containing region level summaries of nifti data, one object per subject in the list.
#' @param path location of files
#' @return List of adjacency matrices for each subject.
#' @export
network.estimate <- function(files, path = NULL){
  # library(huge)
  file.names <- names(files)
  list_adj_mat <- list()
  for(i in 1:length(files)){
    subj_data <- t(files[[file.names[i]]])
    trans_data = huge::huge.npn(subj_data)
    mb_fit = huge::huge(trans_data)
    mb_graph = huge::huge.select(mb_fit)
    adj_mat = as.matrix(mb_graph$refit)
    list_adj_mat[[i]] <- adj_mat # create list of adjacency matrices
  }
  names(list_adj_mat) <- file.names
  return(list_adj_mat)
}
