% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean_map.R
\name{mean_map}
\alias{mean_map}
\title{Mean Map algorithm}
\usage{
mean_map(data_bags, feature_mat, bag_proportions, weight = NULL, ...)
}
\arguments{
\item{data_bags}{A factor vector of the bag factor for the training dataset.
Filled in by llp}

\item{feature_mat}{A matrix of features for each training instance
Filled in by llp}

\item{bag_proportions}{A numerical vector, with each entry representing the proportion
of positive instances for each bag. Must be in the order of the
levels of the data_bags factor.
Filled in by llp}

\item{weight}{an optional vector of weights to be used in the fitting
process. Should be ‘NULL’ or a numeric vector.}

\item{...}{ignored}
}
\description{
The algorithm expects the data organized in data$bag and data$label (the proportions).
The rest is the features vector
}
