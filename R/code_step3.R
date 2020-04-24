#' STEP 3 IN PNC ANALYSIS - ANALYSIS ON NETWORK
#' Perform prediction of connectivity using provided covariates.
#' @param networks list of files containing estimated adjacency matrix for each subject.
#' @param covariate data.frame containing all covariates.
#' @param method used to specify which approach to use for prediction; default option is "SVM" and"randomForest" is also available.
#' @param missing.prop double value used to specify the threshold for filtering out covariates if they are missing more than this specified percentage.
#' @param ID vector of subject IDs that corresponds to the covariate matrix, option to specify if not included in the covariate data.frame.
#' @param newdata option to provide new data to perform the prediction.
#' @return List of two items: predicted variables_selected; predicted is a vector of predicted values and variables_selected
#' is a matrix where a 1 indicates variables (columns) included in the model for a give edge (rows).
#' @export
prediction.analysis <- function(networks, covariates, method = "SVM", missing.prop = 0.05, ID = NULL, newdata = NULL){
  # library(glmnet)
  # library(e1071)
  # library(randomForest)
  # convert to y matrix - # subjects X # edges
  resolution <- dim(networks[[1]])[1]
  y.save <- matrix(0, nrow=length(networks), ncol= resolution^2)
  for(i in 1:resolution){
    for(j in 1:resolution){
      y <- NULL
      for(k in 1:length(networks)){ # loop through every subject
        mat.temp <- networks[[k]]
        y <- c(y, mat.temp[i,j])
      }
      y.save[,((i-1)*resolution+j)] <- y
    }
  }
  # filter out those with too little variability
  idx_remove <- which(colMeans(y.save) < 0.05)
  y.reduced <- y.save[,-idx_remove]

  # match covariates to imaging IDs
  fmri.id <- names(networks) # list of networks is named by ID
  if(!is.null(ID)){
    covariates$subj.id <- ID
  }
  cov.id <- covariates$subj.id # ID is named subj.id - change this later
  if(all(!fmri.id %in% cov.id)){
    print("Error: none of the IDs match")
  }
  common_id = intersect(fmri.id,cov.id)
  idx_cov = match(common_id,cov.id)

  # remove covariates with too many missing values
  missing <- apply(covariates, 2, function(col)sum(is.na(col))/length(col))
  idx_missing <- which(missing > missing.prop)  # remove if missing in more than 5% of subjects
  if(!all(missing < missing.prop)){
    cov_red <- covariates[idx_cov,-idx_missing]
  }else{
    cov_red <- covariates[idx_cov,]
  }
  cov.red <- cov_red[stats::complete.cases(cov_red),] # reduce to complete cases of remaining covariates
  if(is.null(newdata)){newdata <- cov.red}

  # create formula for elastic net
  formula1 <- paste(names(cov.red), collapse = "+")
  fmla <- stats::as.formula(paste( "~", formula1))

  common_id.glm = intersect(fmri.id,cov.red$subj.id)
  idx_fmri = match(common_id.glm,fmri.id)

  # loop through each edge and do elastic net + SVM/random forest
  pred.mat <- matrix(NA, nrow=dim(y.reduced)[2], ncol=nrow(newdata))
  variables_selected <- matrix(NA, nrow=dim(y.reduced)[2], ncol=ncol(cov.red)) # nrow(cov.red)
  colnames(variables_selected) <- colnames(cov.red)
  for(i in 1:dim(y.reduced)[2]){
    y.glm <- as.factor(y.reduced[idx_fmri,i]) # reduce to one edge
    data.fit <- data.frame(y.glm=y.glm, cov.red)

    # pred.save <- rep(NA, length=nrow(cov.red))
    X.glmnet <- stats::model.matrix(fmla, data=data.fit)
    # colnames(variables_selected) <- colnames(X.glmnet) # may be different from colnames(cov.red) because of factors
    X.new <- stats::model.matrix(fmla, data=newdata) # not sure if this will work
    fit_glmnet <- glmnet::cv.glmnet(X.glmnet, data.fit$y.glm, alpha=0.5, family="binomial")

    tmp_coeffs <- stats::coef(fit_glmnet, s = "lambda.1se") # use this lambda to control FDR
    selected <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    for(s in 1:length(selected[,1])){
      if(selected$name[s] %in% colnames(variables_selected)){
        col.sel <- which(selected$name[s] == colnames(variables_selected))
        variables_selected[i,col.sel] <- 1
      }
    }

    # fit SVM or random forest using results of glmnet
    if(dim(selected)[1] > 1){
      var.sel <- selected$name[2:dim(selected)[1]]
      var.temp <- paste(selected$name[2:dim(selected)[1]], collapse="+")
      fmla.glmnet <- stats::as.formula( paste("y.glm ~", var.temp))
      data.fit$y.glm <- as.factor(data.fit$y.glm)

      if(method == "randomForest"){ # random forest for prediction
        p = length(var.sel)
        rf.model <- randomForest::randomForest(fmla.glmnet, data = data.fit, importance = TRUE, ntree = 500, mtry=p/2, proximity=TRUE)
        pred_y = predict(rf.model, newdata, type="vote")[,2]
      }else{ # SVM for prediction
        svm_model <- e1071::svm(fmla.glmnet, data=data.fit, probability = TRUE)
        pred_y <- predict(svm_model, newdata, decision.values = TRUE, probability = TRUE)
      }
    }
    pred.mat[i,] <- pred_y
  }
  predicted.values = t(pred.mat) # subjects X edges matrix of predicted yes/no connectivity - or probability of connected?
  return(list(predicted = predicted.values, variables_selected = variables_selected))
}
