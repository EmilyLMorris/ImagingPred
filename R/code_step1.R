#' STEP 1 IN PNC ANALYSIS - READ IN NIFTI FILE
#' Function summarizes imaging from nifti files into 264 regions with corresponding time series
#' @param file.names provide names of files
#' @param path optional path of files, needs to be included if path to files is not included in file.names
#' @param method method for summarizing within region. Default is averaging all voxels in a given region, "pca" also available if using the first prinicipal component.
#' @param resolution resolution of images. Template provided for 2mm and 3mm resolution, default is 2.
#' @param template option to provide an alternative template if not using standard 2mm or 3mm ones provided.
#' @param PNC option to specify is using PNC raw data, which needs to be converted to standard coordinates.
#' @return List of files containing region level voxel summaries for each subject.
#' @export
read.nifti.files <- function(file.names, path = NULL, method = "average", resolution = 2, template = NULL, PNC = FALSE){
  # library(oro.nifti)
  # library(neurobase) # has readnii function
  subj.list <- list()
  # if(is.null(file.names)){
  #   file.names <- list.files(path)
  # }
  if(PNC){
    node_264 = neurobase::readnii("inst/extdata/total_sphere_3mm.nii") #only need to run this once
    transform = rbind(node_264@srow_x,
                      node_264@srow_y,
                      node_264@srow_z)
    coord_mean <- matrix(NA, nrow=264, ncol=3)
    for(r in 1:264){
      std_coords = t(t(which(node_264==r, arr.ind = TRUE)%*%transform[,1:3])+transform[,4])
      coord_mean[r,] <- apply(std_coords,2,mean)
    }
    for(i in 1:length(file.names)){
      if(!is.null(path)){
        img=neurobase::readnii(file.path(path = path, file.names[i]))
      }else{
        img=neurobase::readnii(file.names[i])
      }
      cat("done \n")
      img_trans =  rbind(img@srow_x,
                         img@srow_y,
                         img@srow_z)
      img_std_coords = t(t(which(img[,,,1]!=100,arr.ind = TRUE)%*%img_trans[,1:3])+img_trans[,4])
      loc <- rep(NA, nrow(img_std_coords))
      for(reg in 1:264){
        k = length(which(node_264 == reg))
        dist <- apply(img_std_coords, 1, function(x) dist(rbind(x,coord_mean[reg,])))
        loc[which(dist %in% sort(dist)[1:k])] <- reg # location of 20 nearest neighbors
      } # fairly slow
      cat("done with nearest neighbor \n")
      idx.img <- which(img[,,,1]!=100, arr.ind = TRUE) # reduce to subset that are actually part of brain
      nifti.mat.red <- matrix(NA, nrow = 264, ncol = n.time)
      for(region in 1:264){
        idx <- which(loc == region)
        idx.new <- idx.img[idx,]
        temp <- t(sapply(1:nrow(idx.new), function(x) img[idx.new[x,1],idx.new[x,2],idx.new[x,3],]))
        nifti.mat.red[region,] <- colMeans(temp) # only implemented the average for PNC data
      }
      subj.list[[i]] <- nifti.mat.red
    }
    names(subj.list) <- file.names
    return(subj.list) # use this loop for PNC data only - need to convert to std coordinates
  }
  # read in 2mm or 3mm map
  if(resolution == 2){
    mask <- neurobase::readnii("inst/extdata/total_sphere.nii")
  }else if(resolution == 3){ # option for 3mm resolution files - use different mask here
    mask <- neurobase::readnii("inst/extdata/total_sphere_3mm.nii")
  }
  if(!is.null(template)){
    mask <- template
  }
  # loop through file names - assuming 1 file per subject
  for(i in 1:length(file.names)){
    if(!is.null(path)){
      nifti.file <- neurobase::readnii(file.path(path = path, file.names[i]))
    }else{
      nifti.file <- neurobase::readnii(file.names[i])
    }
    if(!is.null(dim(nifti.file)[4])){n.time <- dim(nifti.file)[4]}else{n.time <- 1}
    nifti.mat.red <- matrix(NA, nrow = 264, ncol = n.time)
    # 2 options: average and PCA (1st PC) to get the region time series for each of 264 regions - combining across voxels
    if(method == "average"){
      for(j in 1:264){
        idx <- which(mask == j, arr.ind = TRUE)
        temp.mat <- sapply(1:nrow(idx), function(r) nifti.file[idx[r,1],idx[r,2],idx[r,3],])
        temp.mat.avg <- colMeans(t(temp.mat))
        nifti.mat.red[j,] <- temp.mat.avg
      }
    }else if(method == "PCA"){
      for(j in 1:264){
        idx <- which(mask == j, arr.ind = TRUE)
        temp.mat <- sapply(1:nrow(idx), function(r) nifti.file[idx[r,1],idx[r,2],idx[r,3],])
        pca.results <- stats::prcomp(t(temp.mat))
        nifti.mat.red[j,] <- pca.results$rotation[,"PC1"]
      }
    }
    subj.list[[i]] <- nifti.mat.red
  }
  names(subj.list) <- file.names
  return(subj.list)
}

