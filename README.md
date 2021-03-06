ImagingPred
===========

This package performs three main functions.

1.  Reads multiple NIfT1 files and converts them into a dataframe with a
    summary of each region of interst.

2.  Performs network estimation using the R package huge - currently
    implemented with the Meinshausen-Buhlmann option for undirected
    graph estimation.

3.  Performs variable selection on covariates provided to predict
    connectivity at each edge of the network.

### Installation

In order to install the package, several R packages are needed. Some
functions rely on oro.nifti, neurobase, huge, glmnet, e1071, and
randomForest.

The following line of code installs the package from GitHub:

    library(devtools)
    devtools::install_github("EmilyLMorris/ImagingPred")

### Example 1

Here is a simple example for each of the functions. First read in a
NIfTI file, this example may take a few minutes to run.

    library(oro.nifti)
    library(neurobase)
    set.seed(2020)
    # generate an example nifti file
    img = array(rnorm(100*100*100*120), dim = c(100,100,100,120))
    nifti.img = oro.nifti::nifti(img) # generate an example nifti file
    nifti.img
    writenii(nifti.img, "example1") 
    img1 = read.nifti.files("example1")
    network1 = network.estimate(img1)
    table(network1)
    # repeat this for multiple "subjects"
    for(n in 1:10){
      set.seed(2020 + n)
      img = array(rnorm(100*100*100*120), dim = c(100,100,100,120))
      nifti.img = oro.nifti::nifti(img) # generate an example nifti file
      writenii(nifti.img, paste("example", n+1, sep = "")) 
    }

    # function read.nifti.files reads the .nii files and summarizes into 264 regions of interest
    file.names = list.files(getwd())[grep(".nii", list.files(getwd()))]
    img.summary = read.nifti.files(file.names)
    length(img.summary)
    dim(img.summary[[1]])

Function used to generate network estimate for each of the subjects:

    networks.list = network.estimate(img.summary)

Finally we demonstrate how to use the function, prediction.analysis, in
two ways to perform variable selection and to predict connectivity.

    # generate covsariate matrix 
    n = length(networks.list)
    X <- data.frame(matrix(rnorm(n*10), nrow = n), subj.id = file.names)

    results.ex1 = prediction.analysis(networks = networks.list, # network estimates for each subject
                                   covariates = X, # covariate matrix
                                   method = "SVM", # can choose between SVM and random forest 
                                   missing.prop = 0, # can specify a threshold to remove covariates if missing for too many subjects
                                   ID = NULL, # specify subject ID to match networks to covariates 
                                   newdata = NULL) 
    results.ex1$predicted.values
    results.ex1$variables_selected

### Example 2

Example using simulated data with some true signal:

    load("example_pkg.RData")

    results.ex2 = prediction.analysis(networks = networks.list, 
                                   covariates = data.frame(X), 
                                   method = "SVM", 
                                   missing.prop = 0, 
                                   ID = names(networks.list), 
                                   newdata = NULL) 
    dim(results.ex2$predicted)
    length(which(!is.na(results.ex2$variables_selected)))
    results.ex2$index_edges # use to compare to true connections in simulated data
