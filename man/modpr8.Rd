\name{modpr8}
\alias{modpr8}
\docType{data}
\title{
Modified PR8 Flu Strain Vaccine Dataset
}
\description{
Titer results from a 1945 flu vaccine efficacy study were modified for use as a demonstration
of this tool.  Continuous
antibody response values were generated from the originally discrete data as 
\code{log2(R) + rnorm(1,sd=0.5)}, i.e., by
adding random noise to the base-2 log of the categorical response level.  This file contains the
results for the PR8 strain of the flu.
}
\usage{data(modpr8)}
\format{
  A data frame with 1777 observations on the following 3 variables.
  \describe{
    \item{\code{Z}}{Vaccine Status (1 = vaccinated)}
    \item{\code{S}}{Modified PR8 Antibody Response}
    \item{\code{Y}}{Infection Status (1 = infected)}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Salk J, Wilbur J, Menke J (1945).  "A Clinical, Epidemiological and Immunological
  Evaluation of Vaccine against Epidemic Influenza."  \emph{Am. J. Hyg.}, 42, 57 - 93.
}
\examples{
data(modpr8)
}
\keyword{datasets}
