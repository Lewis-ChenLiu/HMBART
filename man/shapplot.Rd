% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shapplot.R
\name{shapplot}
\alias{shapplot}
\title{Generate SHAP Summary Plots for NDE and NIE Effects}
\usage{
shapplot(hmbart_obj)
}
\arguments{
\item{hmbart_obj}{An object returned by the \code{hmbart} function, containing data, covariates, and effect estimates.}
}
\value{
A combined plot displaying SHAP summary plots for NDE and NIE effects, showing the relative importance
of features and their impact on predictions.
}
\description{
This function creates SHAP (SHapley Additive exPlanations) summary plots to visualize feature importance
for Natural Direct Effects (NDE) and Natural Indirect Effects (NIE) based on XGBoost models.
}
\details{
The function fits separate XGBoost models for NDE and NIE estimates using covariates as input features.
It computes SHAP values for each model and visualizes them using summary plots.
SHAP values indicate how much each feature contributes to the prediction, providing interpretability
of the model outputs.
}
