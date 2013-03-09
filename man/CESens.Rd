\name{CESens}
\alias{CESens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sensitivity Analysis For Causal Effect Estimation}
\description{
  An interactive tool for estimating principal stratification-based causal effects 
  complete with sensitivity analysis.  User can define distributions for factual and counterfactual
  trial data and compute causal effect and 95\% confidence intervals.  Data set can be loaded at
  function call or loaded later by interacting with the tool.  Tool provides feedback regarding
  distributional assumptions that would give rise to invalid causal effect estimates.  By varying
  the sensitivity parameter (S(1)|Y(0)=1 distribution), the user can examine the robustness of
  the causal effect estimate.
}
\usage{
CESens(data = NULL, meanS1=0, sdS1=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{  Optional dataset.  File format requirements are:
  \itemize{
    \item Each observation on a separate line.  No header.
    \item Treatment indicator (1/0) first, then biomarker value, then outcome indicator(1/0).
    \item Comma, tab, or space separated.
  }
  }
  \item{meanS1}{  Optional estimated mean for biomarker value for treatment group.  Ignored if a dataset is assigned to \command{data}.}
  \item{sdS1}{  Optional estimated standard deviation for biomarker value for treatment group.  Ignored if a dataset is assigned to \command{data}.}
}
\details{
Currently supported distributions for the factual S(1) biomarker data in the treatment group:
\itemize{
\item Normal:  Parameters are estimated from the data's mean and standard deviation.
\item Gamma:  Parameters are estimated from the data's mean and standard deviation.  If the data
would yield an invalid estimate for either the gamma distributions scale or shape parameters, an
error message is given.
\item Non-parametric:  Kernel density estimating using \command{density}.  
}
For the sensitivity parameter S(1)|Y(0)=1, slider bars are provided to allow the user to 
vary the mean and the standard deviation of the sensitivity parameter's distribution.  
The slider bars are in reference to the mean and standard deviation for S(1).
Supported distributions include:
\itemize{
\item Normal:  Parameters are estimated from the parameter's mean and standard deviation, which
are scaled from those of S(1) using the slider bars.
\item Gamma:  Parameters are estimated from the parameter's mean and standard deviation, which
are scaled from those of S(1) using the slider bars.  If the data
would yield an invalid estimate for either the gamma distributions scale or shape parameters, an
error message is given.
\item Non-parametric:  This option uses the non-parametric density estimation for S(1) offset by the
estimated mean for S(1)|Y(0) from the slider bars.  S(1) must be set to \dQuote{non-parametric} for this
option to be valid.
}
The tool provides a histogram of the S(1) data and a density sketch for S(1) using current assumptions.

The tools also displays the density assumption for the sensitivity parameter S(1)|Y(0) and its
reference distribution, S(1), in the plot marked \dQuote{Distribution of S(1)|Y(0)=1}.  

Causal effect is estimated and displayed in the plot marked \dQuote{Causal Effect (Given S(1))}.  
Confidence intervals (95\%) optionally be calculated using \command{boot} and displayed on the 
causal effect plot.  

The tool also provides the capability to zoom in on a section of the causal effect plot.  When it
zooms, only the causal effect is displayed.  Any confidence intervals must be recalculated.

The widgets that comprise the tool are from \pkg{rpanel}.

Known bug:
\enumerate{
\item The \dQuote{Load Data File} function does not always exit cleanly when attempting to load a file
with an unexpected format.
}
}
\value{
No value is returned by the function.
}
\references{
Bowman, Adrian et al., \dQuote{rpanel:  Simple Interactive Controls 
  for R Functions Using the tcltk Package}, \emph{Journal of Statistical 
  Software}, Jan 2007, Vol 17, Issue 9.
 See also \url{http://www.stats.gla.ac.uk/~adrian/rpanel/}.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
