\name{moderationpredict}
\alias{moderationpredict}
\title{Predict Moderated Causal Mediation Effects using BART}
\description{
  This function predicts the Natural Direct Effect (NDE) and Natural Indirect Effect (NIE) under different moderator values using Bayesian Additive Regression Trees (BART).
}
\usage{
moderationpredict(hmbart_obj, NDE_moderator_names = NULL, NDE_moderator_values = NULL,
                  NIE_moderator_names = NULL, NIE_moderator_values = NULL, 
                  num_trees = 100, n_burn_in = 2000, n_after_burn_in = 500, seed = 42)
}
\arguments{
  \item{hmbart_obj}{An object of class \code{hmbart}, which contains the results from a Heterogeneous Causal Mediation Analysis using BART.}
  \item{NDE_moderator_names}{A character vector specifying the names of moderators for the Natural Direct Effect (NDE).}
  \item{NDE_moderator_values}{A data frame specifying the values of the NDE moderators at which predictions should be made.}
  \item{NIE_moderator_names}{A character vector specifying the names of moderators for the Natural Indirect Effect (NIE).}
  \item{NIE_moderator_values}{A data frame specifying the values of the NIE moderators at which predictions should be made.}
  \item{num_trees}{An integer specifying the number of trees in the BART model. Default is 100.}
  \item{n_burn_in}{An integer specifying the number of burn-in iterations for the BART model. Default is 2000.}
  \item{n_after_burn_in}{An integer specifying the number of post burn-in iterations for the BART model. Default is 500.}
  \item{seed}{An integer specifying the random seed for reproducibility. Default is 42.}
}
\details{
  This function fits separate BART models to estimate the Natural Direct Effect (NDE) and Natural Indirect Effect (NIE) given moderator values. 
  The posterior means and credible intervals (95\%) are returned.
}
\value{
  A list containing two data frames:
  \item{NDE_pred}{A data frame containing the posterior mean and 95\% credible interval for NDE.}
  \item{NIE_pred}{A data frame containing the posterior mean and 95\% credible interval for NIE.}
}
\seealso{
  \code{\link[bartMachine]{bartMachine}}, \code{\link[bartMachine]{bart_machine_get_posterior}}
}
\examples{
# Example usage:
\dontrun{
hmbart_obj <- some_hmbart_model()
moderationpredict(hmbart_obj, 
                  NDE_moderator_names = c("x2", "x5"),
                  NDE_moderator_values = some_dataframe,
                  NIE_moderator_names = c("x3", "x6"),
                  NIE_moderator_values = another_dataframe)
}
}
\author{
  Your Name or Maintainer Name
}
\keyword{bayesian}
\keyword{causal inference}
\keyword{mediation}
