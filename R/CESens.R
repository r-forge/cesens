#Sensitivity analysis tool for causal effect estimation.
#Reference:  Bowman, Adrian et al., "rpanel:  Simple Interactive Controls 
#  for R Functions Using the tcltk Package", _Journal of Statistical 
#  Software_, Jan 2007, Vol 17, Issue 9.
# see also http://www.stats.gla.ac.uk/~adrian/rpanel/


#Load libraries
library(rpanel)
library(tkrplot)
library(boot)

#Define sensitivity analysis tool as a function.
CESens <- function(data=NULL, meanS1=0, sdS1=1) {

#Make the plot comparing the S(1) distribution to that of the sensitivity
#parameter, S(1)|Y(0).
mkplt <- function(x, y, ylbl, refx, refy, tlstr, py0, y1s1d, s1gam, s1y0gam) {
  ylim <- c(min(y, refy), (max(y, refy)*1.6))
#limitu and limitl are the constraints on the density assumptions to keep the
#causal effect estimate within legitimate range.
  limitl <- refy/py0*(y1s1d-1)
  limitu <- refy/py0*(y1s1d+1)
  plot(refx, refy, type="l", ylim=ylim, ylab=ylbl, 
       xlab="S(1)", lty=1)
  lines(x, y, type = "l", lty=2)
  lines(refx, limitl, lty=1, col="red")
  lines(refx, limitu, lty=1, col="red")
#s1gam and s1y0gam refer to flags set in S1vals and S1Y0vals respectively to
#indicate that the gamma distribution cannot be estimated under the current
#assumptions.
  if (s1gam | s1y0gam) {
    charscale <- 1.25
    msglen <- strwidth("Gamma Disabled.", cex=charscale)
    xposn <- 0.95*min(refx)+0.5*diff(range(refx))-0.65*msglen
    yposn <- 0.5*diff(ylim)+min(ylim)
    text(xposn, yposn, pos=4, "Gamma Disabled.", col="blue", cex=charscale)
  }
  title(tlstr)
  title("Keep S(1)|Y(0)=1 inside constraint", cex.main=0.8, line=0.8)
  legend("topleft", #inset=c(0, -0.4), ,
         c("S(1)", "S(1)|Y(0)=1", "Constraint"), lty=c(1,2,1), xpd=TRUE,
         col=c("black", "black","red"))
}

#Get values for the distribution of S(1)|Y(0) based on current assumptions.
S1Y0vals <- function(setdist, setmean=0, setsd=1, dfree=20, ngrid, data=NULL) {
  msg <- NULL
  nogamma <- FALSE
#Normal Distribution
  if (setdist == "n") {
    if (is.null(data)) {
      x <- seq((setmean-3*setsd), (setmean+3*setsd), length=ngrid)
    } else {
      x <- data
    }
    dens <- dnorm(x, mean=setmean, sd=setsd)
  } else {
#t Distribution (currently not enabled)
    if (setdist == "t") {
      if (is.null(data)) {
        if (abs(setmean) < .Machine$double.eps ^ 0.5) {
          x <- seq(qt(0.25, df=dfree), qt(0.975, df=dfree), length=ngrid)
        } else {
          x <- seq(qt(0.025, df=dfree, ncp=setmean), 
                   qt(0.975, df=dfree, ncp=setmean), length=ngrid)
        }
      } else {
        x <- data
      }
      if (abs(setmean) < .Machine$double.eps ^ 0.5) {
        dens <- dt(x, df=dfree)
      } else {
        dens <- dt(x, df=dfree, ncp=setmean)
      }
    } else {
#Gamma Distribution
      if (setdist == "g") {
        if ((setsd < 1) | (setmean <= 0)) {
          nogamma <- TRUE
          count <- 0
          if (setsd < 1) {
            msg <-paste(msg, "S(1)|Y(0) variance is ", sprintf("%3.2f", setsd),
                              " < 1.  ", sep="")
            count <- count + 1
          }
          if (setmean <= 0) {
            msg <- paste(msg, "S(1)|Y(0) mean is ", sprintf("%3.2f", setmean),
                              " <= 0.  ", sep="")
            count <- count + 1
          }
          countmsg <- c("this requirement.", "these requirements.")
          msg <- paste(msg, "\nCannot fit a gamma distribution to ", 
                       countmsg[count], "\n", sep="")
          x <- seq((setmean-3*setsd), (setmean+3*setsd), length=ngrid) 
        } else {
            var <- setsd*setsd
            alpha <- setmean*setmean/var
            beta <- setmean/var
          }
        if (is.null(data)) {
           x <- seq(qgamma(0.025, alpha, beta), 
                   qgamma(0.975, alpha, beta), length=ngrid)
        } else {
            x <- data
          }
        if (nogamma) {
          dens <- rep(-1e6, ngrid)
        } else {
          dens <- dgamma(x, alpha, beta)
        }
      } else {
#Empty else clause in case we want to add other distributions in the future.
      }
    }
  }
  list(x=x, dens=dens, msg=msg, nogamma=nogamma)
}

#Infer values for S(1) based on available data and/or current assumptions.
S1vals <- function(setdist, setmean=0, setsd=1, ngrid, data) {
  msg <- NULL
  nogamma <- FALSE
#Normal Distribution
  if (setdist == "n") {
    if (is.null(data)) {
      x <- seq((setmean-3*setsd), (setmean+3*setsd), length=ngrid)
    } else {
      x <- data
    }
    dens <- dnorm(x, mean=setmean, sd=setsd)
  } else {
#t Distribution (currently not enabled)
    if (setdist == "t") {
      var <- setsd*setsd
      dfree <- 2*var/(var-1)
      if (is.null(data)) {
        if (abs(setmean) < .Machine$double.eps ^ 0.5) {
          x <- seq(qt(0.25, df=dfree), qt(0.975, df=dfree), length=ngrid)
        } else {
          x <- seq(qt(0.025, df=dfree, ncp=setmean), 
                   qt(0.975, df=dfree, ncp=setmean), length=ngrid)
        }
      } else {
        x <- data
      }
      if (abs(setmean) < .Machine$double.eps ^ 0.5) {
        dens <- dt(x, df=dfree)
      } else {
        dens <- dt(x, df=dfree, ncp=setmean)
      }
    } else {
#Gamma Distribution
      if (setdist == "g") {
        if ((setsd < 1) | (setmean <= 0)) {
          count <- 0
          if (setsd < 1) {
            msg <-paste(msg, "Variance is ", sprintf("%3.2f", setsd),
                              " < 1.  ", sep="")
            count <- count + 1
          }
          if (setmean <= 0) {
            msg <- paste(msg, "Mean is ", sprintf("%3.2f", setmean),
                              " <= 0.  ", sep="")
            count <- count + 1
          }
          countmsg <- c("this requirement.", "these requirements.")
          msg <- paste(msg, "Cannot fit a gamma distribution to ", 
                       countmsg[count], "\n", sep="")
          x <- seq((setmean-3*setsd), (setmean+3*setsd), length=ngrid)
          nogamma <- TRUE
        } else {
            var <- setsd*setsd
            alpha <- setmean*setmean/var
            beta <- setmean/var
          }
        if (is.null(data)) {
            x <- seq(qgamma(0.025, alpha, beta), 
                     qgamma(0.975, alpha, beta), length=ngrid)
          } else {
              x <- data
          }
        if (nogamma) {
          dens <- rep(1e6, ngrid)
        } else {
          dens <- dgamma(x, alpha, beta)
        }
      } else {
#Non-parametric Density estimation
        if (setdist == "np") {
          if (is.null(data)) {
            x <- seq(0, 1, length=ngrid)
            dens <- rep(0.5, ngrid)
            msg <- "No data on which to calculate S(1)|Y(0) nonparametric kernel."
          } else {
            np.fit <- density(x=data, bw="sj", kernel="gaussian", n=ngrid,
                              na.rm=TRUE)
            dens <- np.fit$y
            x <- np.fit$x
          }
        }
      }
    }
  }
  list(x=x, dens=dens, msg=msg, nogamma=nogamma)
}

#Estimate values for P(Y(1)|S(1)) based on logistic regression.
#betas are set with logistic regression in read.data whenever a dataset 
#becomes available.
Y1S1vals <- function(beta, setmean, setsd, ngrid, data) {
  if (!is.null(data)) {
    x <- data
  } else {
    x <- seq((setmean-3*setsd), (setmean+3*setsd), length=ngrid)
  }
  dens <- exp(beta[1] + beta[2]*x)/(1+exp(beta[1] + beta[2]*x))
  list(x=x, dens=dens)
}

#Read data
read.data <- function(data) {
#Mean, SD of S(Z=1)=S(1)
  meanS1 <- mean(data[which(data[,1]==1), 2])
  sdS1 <- sd(data[which(data[,1]==1), 2])
#P(Y(Z=0)=1) = P(Y(0)=1) = E(Y(0))
  PY0 <- mean(data[which(data[,1]==0), 3])
#Logistic Regression fit to estimate betas for Y1S1 distribution.
  fit <- glm(data[,3] ~ data[,2], family=binomial(link='logit'), 
             data=data[which(data[,1]==1),])
  beta <- coef(fit)
  S1x.read <- sort(data[which(data[,1]==1),2])
  list(meanS1=meanS1, sdS1=sdS1, PY0=PY0, beta=beta, S1x.read=S1x.read)
}


#List file requirements
filerqmt <- function(object) {
  txt1 <- "Each observation on a separate line.  No header."
  txt2 <- "\nTreatment indicator (1/0) first, then biomarker value, then outcome indicator(1/0)."
  txt3 <- "\nComma, tab, or space separated."
  rp.messagebox(title="File Requirements", paste(txt1, txt2, txt3, sep=""))
  object
}


#Read a data file
readfile <- function(object) {
  read <- FALSE
#Change "\" to "/".  \ is a special character in regular expressions, hence
# the need for multiples of them.
  work <- gsub("\\\\", "/", object$filenm)
#Can't have a trailing "/" in the file name.
  if (grepl("/", substr(work, nchar(work), nchar(work)))) {
    work <- substr(work, 1, (nchar(work)-1))
  }
#Confirm file exists and is not just a directory.
  if (!file_test("-f", object$filenm)) {
    pieces <- strsplit(work, "/")
    npc <- length(pieces[[1]])
    if (npc > 2) {
      tmpstr <- pieces[[1]][1]
      for (i in 2:(npc-1)) {
        tmpstr <- paste(tmpstr, "/", pieces[[1]][i], sep="")
      }
    } else {
      if (npc == 2) {
        tmpstr <- pieces[[1]][1]
      } else {
        tmpstr <- getwd()
      }
    }
    msgtxt <- paste("Cannot find file.  Looked in directory ", tmpstr, 
                    ".", sep="")
    rp.messagebox(title="File Not Found", msgtxt)
  } else {
#Try comma-separated
    new.data <- read.table(work, header=FALSE, stringsAsFactors=FALSE,
                     sep=",")
    if (ncol(new.data) != 3) {
#Try white-space separated.
      new.data <- read.table(work, header=FALSE, stringsAsFactors=FALSE,
                     sep="")
      if (ncol(new.data) != 3) {
        rp.messagebox(title="Unknown Format", 
           "Unknown file format.  Please check file requirements and try again.")
      } else {
        read <- TRUE
      }
    } else {
      read <- TRUE
    }
  }
  if (read) {
#Read the file, put all the data where it belongs, and update the density
#plot.
    update <- read.data(new.data)
    object$meanS1 <- update$meanS1
    object$sdS1 <- update$sdS1
    object$PY0 <- update$PY0
    object$beta <- update$beta
    object$S1$x <- seq(min(na.omit(update$S1x.read)), 
                       max(na.omit(update$S1x.read)), 
                       length=ngrid)
    object$data <- new.data
    object$S1x.read <- update$S1x.read
    object <- replotS1Y0(object)
  }
  object
}

file.nav <- function(panel) {
  require(tcltk)
  fileName <- tclvalue(tkgetOpenFile()) # Very simple, isn't it?
  if (!nchar(fileName)) {
      tkmessageBox(message = "No file was selected!")
  } else {
#      tkmessageBox(message = paste("The file selected was", fileName))
    panel$filenm <- fileName
    panel <- readfile(panel)
  }
  panel
}


#Function for use with integrate.  Returns density value for a given x value.
get.area <- function(x, y) {
  temp.area <- 0
  for (i in 1:(length(x)-1)) {
    y1 <- y[i]
    y2 <- y[(i+1)]
    x1 <- x[i]
    x2 <- x[(i+1)]
    delta <- diff(c(x1, x2))
    temp.area <- temp.area + (min(y1, y2) + 0.5*abs(y2-y1))*delta
  }
  temp.area
}

#Governing function for P(S(1)|Y(0)=1) plot
#Initialize
initS1Y0 <- function(panel) {
  npmsg <- NULL
  ordlab <- "Density"
  titlestr <- "Distribution of S(1)|Y(0)=1"
  newx <- seq(min(na.omit(panel$S1$x)), 
              max(na.omit(panel$S1$x)), length=panel$ngrid)
  npmsg0 <- "No data with which to calculate a nonparametric kernel."
#Calculate teh quantities that go into the causal effect estimate.
#Special cases for non-parametric density handled here so that the right
#error messages are coordinated on the different displays.
  if (panel$distnmS1 == "np") {
    if (is.null(panel$data)) {
      npmsg <- npmsg0
      panel$S1$x <- newx
      panel$S1$dens <- rep(0, panel$ngrid)
    } else {
      temp <- S1vals(panel$distnmS1, panel$meanS1, panel$sdS1,
                         panel$ngrid, panel$S1x.read)
#Interpolate at established grid.
      panel$S1$x <- seq(min(na.omit(panel$S1x.read)), 
                        max(na.omit(panel$S1x.read)), length=panel$ngrid)
      newx <- panel$S1$x
      inter.temp <- approx(x=temp$x, y=temp$dens, xout=panel$S1$x)
      panel$S1$dens <- inter.temp$y
    }
  } else {
    if (is.null(panel$data)) {
      panel$S1 <- S1vals(panel$distnmS1, panel$meanS1, panel$sdS1,
                         panel$ngrid, newx)
    } else {
      newx <- seq(min(na.omit(panel$S1x.read)), 
                  max(na.omit(panel$S1x.read)), length=panel$ngrid)
      panel$S1 <- S1vals(panel$distnmS1, panel$meanS1, panel$sdS1,
                         panel$ngrid, newx)
    }
#For debugging purposes
#rp.messagebox(panel$S1x.read)
panel$debug$x <- panel$data[which(panel$data[,1]==1),2]
#panel$debug$x <- panel$S1$x
panel$debug$y <- panel$S1$dens
#replotDeb(panel)
  }
  panel$sdS1Y0 <- panel$sdS1Y0rto*panel$sdS1
  panel$meanS1Y0 <- panel$meanS1Y0rto*panel$sdS1 + panel$meanS1
  if (panel$distnmS1Y0 == "np") {
    if (is.null(panel$data)) {
      npmsg <- npmsg0
      panel$S1Y0$x <- newx
      panel$S1Y0$dens <- rep(1e6, panel$ngrid)
    } else {
      if (panel$distnmS1 != "np") {
        npmsg <- "Non-parametric density for S(1) must be calculated first."
        panel$S1Y0$x <- newx
        panel$S1Y0$dens <- rep(1e6, panel$ngrid)
      } else {
        dilS1Y0x <- (panel$S1$x-mean(panel$S1$x))*panel$sdS1Y0rto
        shiftS1Y0x <- dilS1Y0x + mean(panel$S1$x) + panel$meanS1Y0rto*panel$sdS1
#Re-normalize density values
        rescaleS1Y0 <- get.area((dilS1Y0x+mean(panel$S1$x)), panel$S1$dens)
        inter.shf <- approx(shiftS1Y0x, panel$S1$dens/rescaleS1Y0, xout=newx, 
                            rule=2)
        panel$S1Y0$dens <- inter.shf$y
      }
    }
  } else {
    panel$S1Y0 <- S1Y0vals(panel$distnmS1Y0, panel$meanS1Y0, panel$sdS1Y0, 
                    panel$dfS1Y0, panel$ngrid, newx)
  }
  panel$Y1S1 <- Y1S1vals(panel$beta, panel$meanY1S1, panel$sdY1S1, 
                        panel$ngrid, newx)
#Update density plot
  mkplt(panel$S1Y0$x, panel$S1Y0$dens, ordlab, panel$S1$x, panel$S1$dens, 
         titlestr, panel$PY0, panel$Y1S1$dens, panel$S1$nogamma,
         panel$S1Y0$nogamma)
  msgs <- paste(panel$S1$msg, panel$S1Y0$msg, npmsg, sep="")
  panel$dspMsg <- msgs
#Update message pane.
  panel <- replotMsg(panel)
#Update S(1) sketch
  panel <- replotData(panel)
#Update causal effect plot
  panel <- replotEffS1(panel)
  panel
}


#Update first panel (technical requirement of rp.panel)
replotS1Y0 <- function(panel) {
  rp.tkrreplot(panel, plotS1Y0)
  panel
}

evals <- function(data) {
#Calculate causal effect.
#  Original expression:
#    panel$Y1S1$dens - panel$S1Y0$dens/panel$S1$dens*panel$PY0
  data[,1] - data[,2]/data[,3]*data[,4]
}


initEffS1 <- function(panel) {
#Plot causal effect.
  ordlab <- "E[Y(1) - Y(0)|S(1)]"
  titlestr <- "Causal Effect (Given S(1))"
  evals.data <- cbind(panel$Y1S1$dens, panel$S1Y0$dens, panel$S1$dens, 
                       panel$PY0)
#Calculate causal effect.
  Evals <- evals(evals.data)
#Handle zoom on plot.
  if (panel$default) {
    ylim <- c(-1,1)
    xlim <- c(min(panel$S1$x), max(panel$S1$x))
  } else {
    ylim <- c(panel$newlim[2], panel$newlim[4])
    xlim <- c(panel$newlim[1], panel$newlim[3])
  }
  plot(panel$S1$x, Evals, ylab=ordlab, main=titlestr, xlab="S(1)", 
          type='l', xlim=xlim, ylim=ylim)
  abline(a=0, b=0, lty=5)
#Warning if causal effect exceeds its support.
  if (length(which(Evals < -1)) + length(which(Evals > 1))) {
    xposn <- min(xlim)-0.05*diff(range(xlim))
    yposn <- min(ylim)+0.5*diff(range(ylim))
    text(xposn, yposn, pos=4,
         "Warning: Chosen parameters make 
          \nCausal Effect exceed constraint.",
         col="red")
  }
#Display bootstrap results if requested.
  if (panel$booting) {
    lines(panel$S1$x, panel$bootCI[,1], lty=3)
    lines(panel$S1$x, panel$bootCI[,2], lty=3)
    lines(panel$S1$x, panel$bootCI[,3], col="green")
    legend("topleft", 
           legend=c("Causal Effect", "95% Conf. Int.", "Bootstrap Est."), 
           lty=c(1,3,1), col=c("black", "black", "green"))
    title(paste("Number of Bootstrap Samples:  ", panel$R, sep=""), line=0.8,
          cex.main=0.8)
  }
  panel
}

#Update panel (technical requirement of rp.panel)
replotEffS1 <- function(panel) {
  rp.tkrreplot(panel, plotEffS1)
  panel
}


initMsg <- function(panel) {
#Display error (and possibly other) messages from the program.
  par(oma=c(0,0,0,0), mar=c(0,4.1,0,0))
  plot(0:1,0:1, type='n', bty='o', ylab="", xlab="", xaxt='n',yaxt='n') 
  mtext("Message: ", side=2, las=1, adj=1)
  if (length(panel$dspMsg)) {
    text(0,0.5, panel$dspMsg, adj=c(0,0.5))
  }
  panel
}


#Update panel (technical requirement of rp.panel)
replotMsg <- function(panel) {
  rp.tkrreplot(panel, plotMsg)
  panel
}

initSensLbl <- function(panel) {
#Display label on plot.
  par(oma=c(0,0,0,0), mar=c(0,0.0,0,0))
  plot(0:1,0:1, type='n', bty='o', ylab="", xlab="", xaxt='n',yaxt='n') 
  text(0.5,0.5, "Sensitivity", font=2)
  panel
}


#Update panel (technical requirement of rp.panel)
replotSensLbl <- function(panel) {
  rp.tkrreplot(panel, lblSens)
  panel
}

initFilesLbl <- function(panel) {
#Display label on plot.
  par(oma=c(0,0,0,0), mar=c(0,0.0,0,0))
  plot(0:1,0:1, type='n', bty='o', ylab="", xlab="", xaxt='n',yaxt='n') 
  text(0.5,0.5, "Files", font=2)
  panel
}


#Update panel (technical requirement of rp.panel)
replotFilesLbl <- function(panel) {
  rp.tkrreplot(panel, lblFiles)
  panel
}

initCIsLbl <- function(panel) {
#Display label on plot.
  par(oma=c(0,0,0,0), mar=c(0,0.0,0,0))
  plot(0:1,0:1, type='n', bty='o', ylab="", xlab="", xaxt='n',yaxt='n') 
  text(0.5,0.5, "Bootstrap Confidence Intervals", font=2)
  panel
}


#Update panel (technical requirement of rp.panel)
replotCIsLbl <- function(panel) {
  rp.tkrreplot(panel, lblCIs)
  panel
}

initPlotLbl <- function(panel) {
#Display label on plot.
  par(oma=c(0,0,0,0), mar=c(0,0.0,0,0))
  plot(0:1,0:1, type='n', bty='o', ylab="", xlab="", xaxt='n',yaxt='n') 
  text(0.5,0.5, "Plot View", font=2)
  panel
}


#Update panel (technical requirement of rp.panel)
replotPlotLbl <- function(panel) {
  rp.tkrreplot(panel, lblPlot)
  panel
}


initData <- function(panel) {
#Sketch the S(1) data in a histogram and add the current assumed density.
  par(oma=c(0,0,0,0), mar=c(0,0,1.1,0))
  if (is.null(panel$data)) {
    plot(0:1,0:1, type='n', main="S(1) Sketch", cex.main=0.8)
    text(0.5,0.5, "No data", adj=c(0.5,0.5))
  } else {
    hist(panel$S1x.read, bty='o', ylab="", xlab="", xaxt='n',yaxt='n', 
         main="S(1) Sketch", cex.main=0.8, freq=FALSE) 
    lines(panel$S1$x, panel$S1$dens)
  }
  panel
}

#Update panel (technical requirement of rp.panel)
replotData <- function(panel) {
  rp.tkrreplot(panel, plotData)
  panel
}

zoom <- function(panel) {
#Zoom on causal effect plot.
#Error checking.  If an illegal value is encountered, reset the parameter to
#the default value.
  errlvl <- 0
  errtxt <- ""
  errspot <- NULL
  deflim <- c(min(panel$S1$x)*0.9, -1.1*(max(panel$PY0, (1-panel$PY0))),
            max(panel$S1$x)*1.1, 1.1*(max(panel$PY0, (1-panel$PY0))))
  comp <- deflim
  spots <- c("lower left x", "lower left y", "upper right x", "upper right y")
  work <- gsub(" *", "", panel$newlim)
  for (i in 1:length(work)) {
    if(is.na(as.numeric(work[i]))) {
      testval <- ""
      if (nchar(work[i]) > 0) {
        for (j in 1:nchar(work[i])) {
          testval <- paste(testval, " ", sep="")
        }
      }
      if (work[i] != testval) {
        errlvl <- errlvl + 1
        errspot <- c(errspot, spots[i])
        errtxt <- "Non-blank character(s) encountered.  If you want to use the default bound, enter blanks only.  Otherwise, enter numerical bounds for the zoom.\n"
      } 
    } else {
      comp[i] <- as.numeric(work[i])
    }
  }
  errtxta <- ""
  if (errlvl > 0) {
    errtxta <- "For "
  }
  count <- errlvl - 1
  while (count > 0) {
    errtxta <- paste(errtxta, errspot[(errlvl - count)], ", ", sep="")
    count <- count - 1
  }
  if (errlvl > 0) {
    errtxta <- paste(errtxta, errspot[errlvl], ":\n", sep="")
  }
  errtxtb <- ""
  if ((comp[1]>comp[3]) | (comp[2]>comp[4])) {
    errlvl <- errlvl + 1
    errtxtb <- "Lower left corner must be lower and to the left of upper right corner.\n"
  }
  errtxtc <- ""
  if ((comp[1]<deflim[1]) | (comp[2]<deflim[2]) | (comp[3]>deflim[3]) | 
      (comp[4]>deflim[4])) {
    errtxtc <- "Cannot zoom out beyond range of data.\n"
    errlvl <- errlvl + 1
  }
  if (errlvl > 0) {
    rp.messagebox(title="Illegal Value", 
      paste(errtxtb, errtxtc, errtxta, errtxt, sep=""))
  }
#If all is well, reset the limits for the plot.
  if (errlvl == 0) {
    panel$default <- FALSE
    panel$newlim <- comp
    panel <- replotEffS1(panel)
  }
  panel
}

reset <- function(panel) {
#Reset the causal effect plot limits to the default values and replot.
  panel$default <- TRUE
  panel <- replotEffS1(panel)
  panel
}

#This function checks the validity of the user-supplied number of bootstrap 
#samples, displays error messages if necessary, and converts it to a number.
checkR <- function(panel) {
  err <- FALSE
  if (is.na(as.numeric(panel$Rtxt))) {
    err <- TRUE
  } else {
    if (as.numeric(panel$Rtxt) <= 0) {
      err <- TRUE
    }
  }
  if (err) {
    rp.messagebox(title="Invalid Value", "Please enter a positive integer.")
  } else {
    panel$R <- as.numeric(panel$Rtxt)
    panel$dspMsg <- paste("Number of bootstrap samples reset to ", panel$R, 
                          ".", sep="")
    panel <- replotMsg(panel)
  }
  panel
}

#This function calculates the Causal Effect statistic and the values it 
#requires for use by the the boot() bootstrap function in R.  
boot.fn <- function(set, indices, ngrid, sdS1Y0rto, meanS1Y0rto, distnmS1Y0, 
                       distnmS1, meanY1S1, sdY1S1, xvals, object) {
  bootrun <- read.data(set[indices,])
  bootx <- xvals
#Calculate required values for causal effect calculation, analogous to what
#Y1S0vals() does.
  bootY1S1 <- Y1S1vals(bootrun$beta, meanY1S1, sdY1S1, ngrid,
                       bootx)
  temp <- S1vals(distnmS1, bootrun$meanS1, bootrun$sdS1, ngrid, 
                 bootrun$S1x.read)
  inter.temp <- approx(x=temp$x, y=temp$dens, xout=bootx, rule=2)
  bootS1d <- inter.temp$y
  boot.sdS1Y0 <-  sdS1Y0rto*bootrun$sdS1
  boot.meanS1Y0 <- meanS1Y0rto*bootrun$sdS1 + bootrun$meanS1
  if (distnmS1Y0 == "np") {
    dilbootS1Y0 <- (bootx-mean(bootx))*sdS1Y0rto
    shiftS1Y0x <- dilbootS1Y0 + mean(bootx) + meanS1Y0rto*bootrun$sdS1
    rescaleS1Y0boot <- get.area((dilbootS1Y0+mean(bootx)), bootS1d)
    inter.shf <- approx(shiftS1Y0x, bootS1d/rescaleS1Y0boot, xout=bootx, 
                        rule=2)
    bootS1Y0 <- list(x=inter.shf$x, dens=inter.shf$y)
  } else {
    bootS1Y0 <- S1Y0vals(distnmS1Y0, boot.meanS1Y0, boot.sdS1Y0, 20, ngrid,
                         bootx)
  }
  bootdens <- cbind(bootY1S1$dens, bootS1Y0$dens, bootS1d, bootrun$PY0)
  bootevals <- evals(bootdens)
#Capture some information for debugging purposes
  object$debug <- list(xs=bootrun$S1x.read, x=bootx, y=bootS1d)
#object <- replotDeb(object)
#Return the causal effect estimate for this bootstrap sample.
  bootevals
}

catchfn <- function(e) {
  rp.messagebox(e)
}

boot.Ints <- function(panel) {
#Governing function for bootstrapping.  Performs the bootstrap and then 
#calculates confidence intervals.
#Error checking
  if (is.null(panel$data)) {
    panel$dspMsg <- "No data set loaded; cannot calculate confidence intervals."
    panel <- replotMsg(panel)
  } else {
    if ((panel$distnmS1Y0 == "np") & (panel$distnmS1 != "np")) {
      panel$dspMsg <- "Non-parametric density for S(1) must be calculated first."
      panel <- replotMsg(panel)
    } else {
#Perform the bootstrap
      bootvals <- boot(data=panel$data, statistic=boot.fn, R=panel$R, 
                        ngrid=panel$ngrid, sdS1Y0rto=panel$sdS1Y0rto, 
                        meanS1Y0rto=panel$meanS1Y0rto, 
                        distnmS1Y0=panel$distnmS1Y0, distnmS1=panel$distnmS1,
                        meanY1S1=panel$meanY1S1, sdY1S1=panel$sdY1S1,
                        xvals=panel$S1$x, object=panel)

#Calculate the confidence intervals.
      intvec <- NULL
      for (i in 1:ncol(bootvals$t)) {
#        citemp <- boot.ci(bootvals, index=i, type="basic")
        citemp <- boot.ci(bootvals, index=i, type="perc")
#        intvec <- rbind(intvec, c(citemp$t0, citemp$basic[,4], citemp$basic[,5]))
        intvec <- rbind(intvec, c(citemp$t0, citemp$perc[,4], citemp$perc[,5]))
      }

#Put everything in the global space so other functions can access it.
      panel$bootCI <- cbind(intvec[,2:3], bootvals$t0)
      panel$booting <- TRUE
#Update density and causal effect plots
      panel <- replotS1Y0(panel)
      panel <- replotEffS1(panel)
      panel$booting <- FALSE
    }
  }
  panel
}

initDeb <- function(panel) {
#display used for debugging purposes
  par(oma=c(0,0,0,0), mar=c(3,3,1.1,0))
  hist(panel$debug$x, freq=FALSE)
  lines(panel$debug$x, panel$debug$y)
  plot(panel$debug$x, panel$debug$y, type='l')
  panel
}

#Update panel (technical requirement of rp.panel)
replotDeb <- function(panel) {
  rp.tkrreplot(panel, plotDebug)
  panel
}




debug.box <- function(panel) {
#Function for debugging purposes.
#  ce <- panel$Y1S1$dens - panel$S1Y0$dens/panel$S1$dens*panel$PY0
#  limitl <- panel$S1$dens/panel$PY0*(panel$Y1S1$dens-1)
#  limitu <- panel$S1$dens/panel$PY0*(panel$Y1S1$dens+1)
#  rp.messagebox(paste("\nY1S1", panel$Y1S1$dens[which(panel$S1$x >10)], 
#                      "\nS1Y0", panel$S1Y0$dens[which(panel$S1$x >10)],
#                      "\nS1(dens)", panel$S1$dens[which(panel$S1$x >10)], 
#                      "\nS1", panel$S1$x[which(panel$S1$x >10)], 
#                      "\nresult", ce[which(panel$S1$x >10)], 
#                      "\nupper", limitu[which(panel$S1$x >10)], 
#                      "\nlower", limitl[which(panel$S1$x >10)]))
rp.messagebox(panel$debug$bootevals)
  panel
}


#Grid length.  All quantities used in the causal effect estimation have to
#be calculated at the same points.
ngrid <- 1024


#Set up initial values.  If the function is given a dataset, calculate the
#values based on it.  Otherwise, use some defaults.
if(!is.null(data)) {
  init <- read.data(data)
} else {
  init <- list(meanS1=meanS1, sdS1=sdS1, PY0=0.5, beta=c(0, 0.5), 
                S1x.read=NULL)
}
init <- c(meanY1S1=0, sdY1S1=1, meanS1Y0=0, sdS1Y0=1, init, distnmS1Y0="n", 
           distnmS1="n", meanS1Y0rto=0, sdS1Y0rto=1, default=TRUE, ngrid=ngrid) 
init$newlim=c(1, -0.5, 100, 0.5)
init$dspMsg <- ""
init$R <- 100


Y1S1init <- Y1S1vals(init$beta, init$meanY1S1, init$sdY1S1, init$ngrid, init$S1x.read)
S1Y0init <- S1Y0vals(init$distnmS1Y0, init$meanS1Y0, init$sdS1Y0, dfree=20, 
                    init$ngrid, init$S1x.read)
S1init <- S1vals(init$distnmS1, init$meanS1, init$sdS1, init$ngrid, init$S1x.read)

debug <- list(xs=init$S1x.read, x=S1init$x, y=S1init$dens)


#Most of the widgets created write the output of the widget to a global 
#variable.  The R command checker thinks these variables have no binding.
#So it is necessary to NULLify them prior to creating the control panel to
#assuage the checker.  That's the purpose of this next command.

plotS1Y0 <- plotEffS1 <- plotMsg <- plotData <- plotDebug <- plotEffS1 <- NULL
plotMsg <- plotData <- plotS1Y0 <- distnmS1 <- distnmS1Y0 <- NULL
meanS1Y0rto <- sdS1Y0rto <- filenm <- newlim <- Rtxt <- NULL
replotSensLbl <- replotFilesLbl <- replotCIsLbl <- replotPlotLbl <- NULL
lblSens <- lblFiles <- lblCIs <- lblPlot <- NULL

#Initialize the window.
master.panel <- rp.control("Simple Sensitivity Analysis", ngrid=init$ngrid, 
                  meanY1S1=init$meanY1S1, sdY1S1=init$sdY1S1, 
                  beta=init$beta, Y1S1=Y1S1init,
                  meanS1Y0=init$meanS1Y0, sdS1Y0=init$sdS1Y0, dfS1Y0=20, 
                  distnmS1Y0=init$distnmS1Y0, S1Y0=S1Y0init,
                  meanS1Y0rto=init$meanS1Y0rto, sdS1Y0rto=init$sdS1Y0rto,
                  meanS1=init$meanS1, sdS1=init$sdS1, 
                  distnmS1=init$distnmS1, S1=S1init,
                  PY0=init$PY0, S1x.read=init$S1x.read, data=data, 
                  default=init$default, newlim=init$newlim,
                  bootCI=NULL, booting=FALSE, dspMsg=init$dspMsg,
                  debug=debug, R=init$R)

#Labels
rp.tkrplot(master.panel, lblSens, plotfun=initSensLbl, vscale=0.05, hscale=0.9,
           row=0, column=0, columnspan=2, sticky="ew", rowspan=1)

rp.tkrplot(master.panel, lblFiles, plotfun=initFilesLbl, vscale=0.05, hscale=0.9,
           row=0, column=2, columnspan=2, sticky="ew", rowspan=1)

rp.tkrplot(master.panel, lblCIs, plotfun=initCIsLbl, vscale=0.05, hscale=0.9,
           row=2, column=2, columnspan=2, sticky="sew", rowspan=1)

rp.tkrplot(master.panel, lblPlot, plotfun=initPlotLbl, vscale=0.05, hscale=0.9,
           row=5, column=2, columnspan=2, sticky="ew", rowspan=1)

#Summary Measure Panel
rp.tkrplot(master.panel, plotEffS1, plotfun=initEffS1, 
           row=8, column=2, columnspan=2, sticky="ns", rowspan=1)

#Message Panel
rp.tkrplot(master.panel, plotMsg, plotfun=initMsg, vscale=0.2, hscale=1.8,
           row=9, column=0, columnspan=4, sticky="ew", rowspan=1)

#Sketch the data
#vscale was 0.33
rp.tkrplot(master.panel, plotData, plotfun=initData, vscale=0.4, hscale=0.36,
           row=5, column=0, sticky="n", rowspan=3)

#Plot area for debugging
#rp.tkrplot(master.panel, plotDebug, plotfun=initDeb, row=8, column=4)

#Draw the first panel.
rp.tkrplot(master.panel, plotS1Y0, plotfun=initS1Y0,
           row=8, column=0, columnspan=2, sticky="ns", rowspan=1)



#Rest of the widgets
#Select S1 distribution
rp.radiogroup(master.panel, distnmS1,
                   labels=c("Normal", "Gamma", "Non-Parametric"),
                   title="S(1) Distribution Options",
                   values=c("n", "g", "np"), action=replotS1Y0, 
                   row=1, column=0, sticky="ns", rowspan=4)

#Select S1Y0 distribution
rp.radiogroup(master.panel, distnmS1Y0,
                   labels=c("Normal", "Gamma", "Shifted S(1) Non-parm."),
                   title="S(1)|Y(0)=1 (Sensitivity) Options",
                   values=c("n", "g", "np"), action=replotS1Y0, 
                   row=1, column=1, sticky="ns", rowspan=4)

#Offest S1Y0 mean by a value set by the slider that is scaled by the value
#set for the S1Y0 standard deviation.
rp.slider(master.panel, meanS1Y0rto, from=-10, to=10, action=replotS1Y0, 
                      title="S(1)|Y(0) Center Shift (in S(1) sd's)", 
                      showvalue=TRUE,
                      row=6, column=1, sticky="ew", resolution=0)
#Set the S1Y0 standard deviation to a value scaled by the S(1) standard
#deviation.
rp.slider(master.panel, sdS1Y0rto, from=0.1, to=10, action=replotS1Y0, 
                      title="S(1)|Y(0) Range Stretch", 
                      resolution=0,
                      showvalue=TRUE, row=7, column=1, sticky="ew")

#Degrees of freedom no longer needed since t distribution no longer
# supported.
#rp.doublebutton(master.panel, dfS1Y0, step=1, 
#                   title="DoF (for t)", action=replotS1Y0,
#                   showvalue=TRUE, row=3, column=1, sticky="n", range=c(1,NA),
#                   showvaluewidth=2)

#Text entry box for file name to read in.
#rp.textentry(master.panel, filenm, action=readfile, 
#                   labels=c("File Name: "), initval="", row=0, column=2, 
#                   columnspan = 1, sticky = "n", width=15, height=5)

#Button to request file format requirements.
rp.button(master.panel, action=filerqmt, title="File Requirements",
                   row=1, column=3, sticky="nw")

#Button to initiate navigation to a file
rp.button(master.panel, action=file.nav, title="Load Data File", row=1, 
                   column=2, sticky="n")

#Button to request bootstrap confidence intervals.
rp.button(master.panel, action=boot.Ints, 
                   title="Calculate", row=3, column=2, 
                   sticky="s", columnspan=2)

#Text entry for number of bootstrap samples
rp.textentry(master.panel, Rtxt, action=checkR, 
                    labels=c(paste("Boot Samples:  (Default ", init$R, ")",
                                   sep="")), 
                    width=6, height=105, initval=init$R, row=4, column=2, 
                    columnspan=2, sticky="s")

#Button for debugging
#rp.button(master.panel, action=debug.box, title="Debug Info", row=6, 
#                   column=4, sticky="sw")

#Text entry boxes for zoom box boundary definitions (for causal effect plot)
rp.textentry(master.panel, newlim, action=zoom, labels=c("lower left x", 
                     "lower left y", "upper right x", "upper right y"), 
              title="Zoom Box Corners\n'Enter' to Apply", 
              initval=c("", "", "", ""), row=6, column=2, columnspan=1,
              sticky="e", width=5, height=5, rowspan=2)

#Button to reset causal effect plot limits.
rp.button(master.panel, action=reset, title="Default Plot Size", row=6, column=3,
              sticky="w")

}






