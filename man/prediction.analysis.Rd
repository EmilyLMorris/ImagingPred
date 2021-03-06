% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/code_step3.R
\name{prediction.analysis}
\alias{prediction.analysis}
\title{STEP 3 IN PNC ANALYSIS - ANALYSIS ON NETWORK
Perform prediction of connectivity using provided covariates.}
\usage{
prediction.analysis(
  networks,
  covariates,
  method = "SVM",
  missing.prop = 0.05,
  ID = NULL,
  newdata = NULL
)
}
\arguments{
\item{networks}{list of files containing estimated adjacency matrix for each subject.}

\item{covariates}{data.frame containing all covariates.}

\item{method}{used to specify which approach to use for prediction; default option is "SVM" and"randomForest" is also available.}

\item{missing.prop}{double value used to specify the threshold for filtering out covariates if they are missing more than this specified percentage.}

\item{ID}{vector of subject IDs that corresponds to the covariate matrix, option to specify if not included in the covariate data.frame.}

\item{newdata}{option to provide new data to perform the prediction.}
}
\value{
List of three items: predicted variables_selected; predicted is a vector of predicted values; variables_selected
is a matrix where a 1 indicates variables (columns) included in the model for a give edge (rows); index_edges is the index of which edges were modeled,
some edges may be removed due to lack of variability across subjects.
}
\description{
STEP 3 IN PNC ANALYSIS - ANALYSIS ON NETWORK
Perform prediction of connectivity using provided covariates.
}
