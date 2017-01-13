#set working directory (no trailing slash)
setwd("C:/path/to/working/directory")

############ settings - change with care ##############
in_file <- "your_data.csv" #name of file to import

meancent <- 1 #set to 0 to disable mean centering of predictors
del_outs <- 1 #set to 0 to disable removal of outliers
heteros <- 1  #set to 0 to disable heterscedasticity-consistent SEs
jntech <- 1 #set to 0 to disable johnson-neyman table
img_type <- "png" #type of moderation image to create (svg, emf or png)
img_alias <- "cairo" #create antialiased png with Cairo, or none
img_dim <- 10 #height and width of png
img_fontsize <- 24 #size of font for moderation plot (comes out much smaller)
use_covars <- 1 #set to 0 to disable use of covariates in model

#set to names of columns in your dataset
iv_var <- "name_of_iv_column"   #x
mod_var <- "name_of_mod_column" #m
dv_var <- "name_of_dv_column"   #y

#set to names of covariates in your dataset
#leave empty if you have no covariates
cov1 <- "name_of_cov1_var"
cov2 <- "name_of_cov2_var"
cov3 <- "name_of_cov3_var"

#set descriptive names of iv, dv and moderator variables 
#(will be used in graphs instead of column variable names)
iv_var_desc <- "Description IV"
mod_var_desc <- "Description Moderator"
dv_var_desc <- "Descrition DV"
mod_line_desc <- c("high MOD","mean","low MOD")

#set to names of control variables
cov1_desc <- "Description Covariate 1"
cov2_desc <- "Description Covariate 2"
cov3_desc <- "Description Covariate 3"

################### dont change below ###################

#install missing libraries
list.of.packages <- c("foreign", "stringr", "lmtest", "zoo", "sandwich", 
                      "Cairo", "devEMF", "psych", "rockchalk", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load required libraries
library(foreign)
library(stringr)
library(lmtest) #needed for bp
library(zoo)
library(sandwich)
library(Cairo)
library(psych)
library(MASS)

############# Functions ##############

options(scipen=999) #turn off scientific notation
#Apa rule: If a value has the potential to exceed 1.0, use the leading zero, 
#if a value can never exceed 1.0, do not use the leading 
num_dec <- function(val) { formatC(round(as.numeric(val), 4), format='f', digits=4) }
num_apa <- function(val) { sub("^(-?)0.|^(-?) 0.", "\\1.", num_dec(val)) }
no_zero <- function(val) { sub("^(-?)0.|^(-?) 0.", "\\1.", val) }
print_gap <- 5 #default gap between columns in output tables

add_stars <- function(p, cutoffs = c(0.1, 0.05, 0.01, 0.001)) {
  stopifnot(length(cutoffs) == 4)
  if (inherits(p, c('matrix', 'data.frame')) && length(dim(p)) == 2) {
    apply(p, c(1,2), add_stars, cutoffs = cutoffs)
  } else {
    if (length(p) > 1) { 
      sapply(p, add_stars, cutoffs = cutoffs)
    } else {
      ifelse(p > cutoffs[1], paste(no_zero(p),'   '), 
             ifelse(p > cutoffs[2], paste(no_zero(p),'.  '), 
                    ifelse(p > cutoffs[3], paste(no_zero(p),'*  '), 
                           ifelse(p > cutoffs[4], paste(no_zero(p),'** '), paste(no_zero(p),'***')))))
    }
  }
}

img_save <- function(num, img_height) {
  imgname <- paste(filename,"_",num,sep="")
  if (img_type == "svg") {
    svg(paste(paste(imgname,sep="_"), "svg", sep="."))
  } else if (img_type == "emf") {
    library(devEMF)
    emf(paste(paste(imgname,sep="_"), "emf", sep="."), 
        width=img_dim, height=img_height, 
        family="Times", pointsize=img_fontsize)
  } else {
    #reducing image heigt by percentage gets rid of top margin after squaring
    png(paste(paste(imgname,sep="_"), "png", sep="."), 
        units='in', res=96, type=img_alias, width=img_dim, 
        height=img_height, pointsize=img_fontsize
    )
  }
}

#plots the moderation
modplot<-function(plotdata){
  #set parameter for outer margins of plot
  #bottom,left,top,right
  op <- par(mar=c(3,2.9,0,0), family="serif")
  
  #adjustments for legend outside of plot
  plot.new() #must be called to draw empty test legend
  l <- legend(0, 0, bty='n', mod_line_desc, plot=F, lwd=0, cex=0.9, 
              title=paste(mod_var_desc,"(M)"))
  #l$rect, dimensions of rectange, w=width, h=height
  #calculate right margin width in ndc
  w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
  w = w + 0.05
  #omd values 1 and 3 define in percentages of the device region
  #the starting points of the x and y axes, while values 
  #2 and 4 define the end points (left and right), 1 = 100%
  op <- par(omd=c(0, 1-w, 0, 1), pty="s")
  
  #assign the plotting data from the moderation outcome
  numplot_data1 <- as.numeric(plotdata[,1])
  numplot_data3 <- as.numeric(plotdata[,3])
  plottest <- cbind(numplot_data3)
  
  if (iv_dicot < 1) {
    lowline <- plottest[1:3,]
    medline <- plottest[4:6,]
    hiline <- plottest[7:9,]
    lty_legend <- c(1,2,3)
  } else {
    lowline <- c(plottest[1],plottest[3],plottest[5])
    hiline <- c(plottest[2],plottest[4],plottest[6])
    lty_legend <- c(1,2)
    mod_line_desc <- mod_line_desc[-2] #remove mean
  }
  
  #plot an empty plot
  plot(NULL, ylim=c(min(numplot_data3)-0.1,max(numplot_data3)+0.1), 
       xlim=c(0.9,3.1), type="n", xaxt="n", yaxt="n", ylab="", xlab="")
  
  #draw legend outside of box
  op <- par(xpd=TRUE) #allows legend to be outside of plot
  legend(par('usr')[2], par('usr')[4], mod_line_desc, lty=lty_legend, lwd=2, cex=0.9,
         title="", inset=c(-0.3,0), bty='n', xpd=NA, seg.len=2.2)
  #make only the legend title bold
  #op <- par(font=2)
  legend(par('usr')[2], par('usr')[4], "", lwd=0, cex=0.9, 
         title=paste("  ",mod_var_desc,"(M)"), bty='n', xpd=NA)
  
  #draw horizontal grid lines, axis and moderation lines
  op <- par(xpd=FALSE, mgp=c(1,0.6,0))
  grid(NA, NULL, col = "#999999", lty = "dotted")
  axis(1, at=1:3, labels=c("-1 SD", "Mean", "+1 SD"), tck=-0.025, cex.axis=0.9)
  axis(2, cex.axis=0.9, las=1, tck=-0.025)
  mtext(1, text=iv_var_desc, line=2, font=1)
  mtext(2, text=dv_var_desc, line=2, font=1)
  lines(lowline, lty="dotted", lwd=2)
  if (iv_dicot < 1) {
    lines(medline, lty="longdash", lwd=2)
  }
  lines(hiline, lty="solid", lwd=2)
  
  #superimpose scatterplot points
  if (iv_dicot < 1) {
    op <- par(new=TRUE)
    plot(noout$IV,noout$DV, 
         xlim=c(min(numplot_data1),max(numplot_data1)),
         #ylim=c(min(plotdata[,4])-0.1,max(plotdata[,4])+0.1),
         xaxt="n", yaxt="n", ylab="", xlab="", col="#555555")
  }
  par(op)
}

#rockchalk (not called, only for testing purposes)
jnplot <- function() {
  library(rockchalk)
  op <- par(family="serif", pty="s", mar=c(2,4.2,0.1,0.1))
  jn <- plotSlopes(model, plotx="IV", modx="MOD")
  jn_psats <- testSlopes(jn)
  plot(jn_psats)
  par(op)
}

#plots the jn regions of significance
custom_jnplot <- function(jn_root1, jn_root2, jnvals) {
  #bottom,left,top,right
  op <- par(xpd=FALSE, mar=c(2.9,3.4,0.1,0.1), mgp=c(1.5,0.6,0), family="serif")
  fillcolour <- rgb(0, 0, 0, max=255, alpha=40)
  jn_roots <- cbind(jn_root1, jn_root2)
  
  x_range = range(jnvals[,1])
  y_range = range(c(jnvals[,6],jnvals[,7]))
  
  jn_llci_data <- as.data.frame(cbind(x = jnvals[,1], y = jnvals[,6]))
  jn_ulci_data <- as.data.frame(cbind(x = jnvals[,1], y = jnvals[,7]))
  jn_meff_data <- as.data.frame(cbind(x = jnvals[,1], y = jnvals[,2]))
  
  #figure out if what significant regions we have hi/low
  #lower only
  if (!is.na(jn_root1) && is.na(jn_root2)) {
    jn_llci_sub <- subset(jn_llci_data, x <= jn_root1)
    jn_ulci_sub <- subset(jn_ulci_data, x <= jn_root1)
    plot_regions <- 1
    #upper only
  } else if (is.na(jn_root1) && !is.na(jn_root2)) {
    jn_llci_sub <- subset(jn_llci_data, x >= jn_root2)
    jn_ulci_sub <- subset(jn_ulci_data, x >= jn_root2)
    plot_regions <- 2
    #both lower and upper
  } else if (!is.null(jn_root1) && !is.null(jn_root2)) {
    jn_llci_sub_lo <- subset(jn_llci_data, x <= jn_root1)
    jn_ulci_sub_lo <- subset(jn_ulci_data, x <= jn_root1)
    jn_llci_sub_hi <- subset(jn_llci_data, x >= jn_root2)
    jn_ulci_sub_hi <- subset(jn_ulci_data, x >= jn_root2)
    plot_regions <- 3
  }
  
  #smooth the lines
  lw1 <- loess(y ~ x, data=jn_llci_data)
  lw2 <- loess(y ~ x, data=jn_ulci_data)
  lw3 <- loess(y ~ x, data=jn_meff_data)
  if (plot_regions < 3) {
    lw4 <- loess(y ~ x, data=jn_llci_sub)
    lw5 <- loess(y ~ x, data=jn_ulci_sub)
  } else {
    lw4 <- loess(y ~ x, data=jn_llci_sub_lo)
    lw5 <- loess(y ~ x, data=jn_ulci_sub_lo)
    lw6 <- loess(y ~ x, data=jn_llci_sub_hi)
    lw7 <- loess(y ~ x, data=jn_ulci_sub_hi)
  }
  
  jn_ylab_desc <- paste("Conditional effect of", iv_var_desc, "on", dv_var_desc, "\n")
  #empty plot with axis set up properly
  plot(y ~ x, data=jn_llci_data, ylim=y_range, xlim=x_range, 
       xlab=paste(mod_var_desc,"(M)"), ylab=jn_ylab_desc, 
       cex.axis=0.9, las=1, tck=-0.025, type="n")
  
  #sort values to join dots
  j1 <- order(jn_llci_data$x)
  j2 <- order(jn_ulci_data$x)
  j3 <- order(jn_meff_data$x)
  
  if (plot_regions < 3) {
    j4 <- order(jn_llci_sub$x)
    j5 <- order(jn_ulci_sub$x)
    
    polygon(x = c(jn_llci_sub$x, rev(jn_llci_sub$x)), 
            y = c(lw5$fitted[j5], rev(lw4$fitted[j4])), 
            col=fillcolour, border=0)
  } else {
    j4 <- order(jn_llci_sub_lo$x)
    j5 <- order(jn_ulci_sub_lo$x)
    j6 <- order(jn_llci_sub_hi$x)
    j7 <- order(jn_ulci_sub_hi$x)
    
    polygon(x = c(jn_llci_sub_lo$x, rev(jn_llci_sub_lo$x)), 
            y = c(lw5$fitted[j5], rev(lw4$fitted[j4])), 
            col=fillcolour, border=0)
    
    polygon(x = c(jn_llci_sub_hi$x, rev(jn_llci_sub_hi$x)), 
            y = c(lw7$fitted[j7], rev(lw6$fitted[j6])), 
            col=fillcolour, border=0)
  }
  
  lines(jn_llci_data$x[j1], lw1$fitted[j1], lty="dashed")
  lines(jn_ulci_data$x[j2], lw2$fitted[j2], lty="dashed")
  lines(jn_meff_data$x[j3], lw3$fitted[j3], lty="solid", lwd=2)
  
  abline(h =0, v=jn_roots, lty="dotted")
  text(jn_roots, min(y_range), labels=round(jn_roots,2), adj=-0.4, cex=0.9)
  
  legend("topleft", legend=c("Point estimate", "95% Conf. intervals", "Region(s) of significance"),
         lty=c(1,2,NA), lwd=c(2,1,NA), col=c(1,1,fillcolour), 
         bg="white", cex=0.9, pch=c(NA,NA,NA), fill=c(NA,NA,fillcolour), 
         border = c(NA,NA,1))
  par(op)
}

################# test assumptions first ####################
#read data from csv file
if (file.exists(in_file)) {
  if (use_covars > 0) {
    master <- read.csv(in_file, sep=",")[ ,c(iv_var, mod_var, dv_var, cov1, cov2, cov3)]
  } else {
    master <- read.csv(in_file, sep=",")[ ,c(iv_var, mod_var, dv_var)]
  }
  
  #define file name from column names for later saving of files
  var_filename <- paste(colnames(master[1]),colnames(master[2]),colnames(master[3]))
  filename <- str_replace_all(var_filename, " ", "--")
  
  if (use_covars > 0) {
    names(master) <- c("IV", "MOD", "DV", "Z1", "Z2", "Z3")
  } else {
    names(master) <- c("IV", "MOD", "DV")
  }
  
  #test for dicothimous dv (y)
  if (length(unique(master$DV)) == 2) {
    dv_dicot <- 1 
  } else {
    dv_dicot <- 0
  }
  
  #test for dicothimous dv (y)
  if (length(unique(master$IV)) == 2) {
    iv_dicot <- 1 
  } else {
    iv_dicot <- 0
  }
  
  ##################### outlier analysis #####################
  # run regression for outlier check using master data (exclude covars here)
  outlier_master <- as.data.frame(cbind(master$IV, master$DV, master$MOD))
  lm_out <- lm(DV ~ IV*MOD, na.action=na.exclude, data=master)
  
  # read https://www3.nd.edu/~rwilliam/stats2/l24.pdf
  # mahalanobis distance
  mahal <- mahalanobis(outlier_master, 
                       colMeans(outlier_master, na.rm = T), 
                       cov(outlier_master, use="complete"))
  cutoff <- qchisq(1-.001, ncol(outlier_master))
  badmahal <- as.numeric(mahal > cutoff)
  
  # leverage
  k <- length(lm_out$coefficients)-1 #nubmer of variables (iv, mod and dv)
  leverage <- hatvalues(lm_out)
  cutleverage <- (2*k+2) / nrow(outlier_master)
  badleverage <- as.numeric(leverage > cutleverage)
  
  # cooks distnace
  cooks <- cooks.distance(lm_out)
  cutcooks <- 4 / (nrow(outlier_master) - k - 1)
  badcooks <- as.numeric(cooks > cutcooks)
  
  #find out which cases are outliers, add 1 to account for header
  outlier_mat <- cbind(badmahal,badleverage,badcooks)
  outlier_id <- which(grepl(3, rowSums(outlier_mat)))
  
  # overall outliers, add them up!
  totalout <- badmahal + badleverage + badcooks
  
  ##################### check assumptions #####################
  # create new dataset without outliers
  if (del_outs > 0) {
    noout <- subset(master, totalout < 3)
  } else {
    noout <- master
  }
  
  # means center predictors
  if (meancent > 0) {
    noout$IV <- as.numeric(scale(noout$IV, scale=F, center=T))
    noout$MOD <- as.numeric(scale(noout$MOD, scale=F, center=T))
  }
  
  # run regression on new outlier-free dataset (noout) -- factor for binary variables
  if (use_covars > 0) {
    model <- lm(DV ~ IV*MOD+Z1+Z2+Z3, na.action=na.exclude, data=noout)
  } else {
    model <- lm(DV ~ IV*MOD, na.action=na.exclude, data=noout)
  }
  
  pdf(paste(paste(filename,sep="_"), "pdf", sep="."), paper="a4", 
      title="Assumption checks", family="serif")
  par(mfrow = c(2,2))
  
  # create standardized residuals plot
  unstandPred <- predict(model)
  unstandResid <- resid(model)
  standPred <- (unstandPred-mean(unstandPred))/sd(unstandPred)
  standResid <- (unstandResid-mean(unstandResid))/sd(unstandResid)
  
  # create residuals histogram
  tmp <- density(standResid)
  hist(standResid, freq=FALSE, 
       main="Histogram of Standardized Residuals",
       xlab="Regression Standardized Residual", ylim=c(0, max(tmp$y)+0.1), 
       xlim=c(min(tmp$x),min(tmp$x)*-1))
  curve(dnorm, add=TRUE, col="blue", lty=1, lwd=1)
  
  # get probability distribution for residuals and plot PP
  probDist <- pnorm(standResid)
  plot(ppoints(length(standResid)), sort(probDist), 
       main = "Normal P-P Plot of Standardized Residual",
       sub = dv_var_desc,
       xlab = "Observed Probability", 
       ylab = "Expected Probability")
  # add diagonal line
  abline(0,1)
  
  # visual inspection for heteroscedasticity
  plot(standResid, standPred,
       main = "Standardized Residuals Plot", 
       xlab = "Standardized Residual values", 
       ylab = "Standardized Predictor values")
  abline(0,0)
  abline(v=0)
  dev.off()
  
  # do Breusch-Pagan test for heteroscedasticity
  # if significant you have heteroscedasticity
  bp <- bptest(model)
  
  ########### Graphical moderation data #############
  
  #estimated covariate means
  mean_iv <- mean(noout$IV)
  mean_mod <- mean(noout$MOD)
  
  if (use_covars > 0) {
    means_Z1 <- mean(noout$Z1)
    means_Z2 <- mean(noout$Z2)
    means_Z3 <- mean(noout$Z3)
  }
  
  sd_mod <- sd(noout$MOD)
  sd_iv <- sd(noout$IV)
  
  #calculate the plotting data the hard way
  intercept <- model$coefficients[1]
  beta_iv <- model$coefficients[2]
  beta_mod <- model$coefficients[3]
  if (use_covars > 0) {
    beta_Z1 <- model$coefficients[4]
    beta_Z2 <- model$coefficients[5]
    beta_Z3 <- model$coefficients[6]
    beta_int <- model$coefficients[7]
  } else {
    beta_int <- model$coefficients[4]
  }
  
  #low modsd and low ivsd (-1SD)
  lowlow <- beta_iv * (mean_iv-sd_iv) + beta_mod * (mean_mod-sd_mod) + 
    beta_int * ((mean_iv-sd_iv) * ((sd_mod-mean_mod)*-1)) + intercept
  
  #low modsd and med ivsd (-1SD)
  lowmed <- (beta_iv * mean_iv) + beta_mod * (mean_mod-sd_mod) + 
    beta_int * (mean_iv * ((sd_mod-mean_mod)*-1)) + intercept
  
  #low modsd and hi ivsd (-1SD)
  lowhi <- beta_iv * (mean_iv+sd_iv) + beta_mod * (mean_mod-sd_mod) + 
    beta_int * ((mean_iv+sd_iv) * ((sd_mod-mean_mod)*-1)) + intercept
  
  #med modsd and low ivsd (mean)
  medlow <- beta_iv * (mean_iv-sd_iv) + (beta_mod * mean_mod) + 
    beta_int * ((mean_iv-sd_iv) * mean_mod) + intercept
  
  #med modsd and med ivsd (mean)
  medmed <- (beta_iv * mean_iv) + (beta_mod * mean_mod) + 
    beta_int * (mean_iv * mean_mod) + intercept
  
  #med modsd and hi ivsd (mean)
  medhi <- beta_iv * (mean_iv+sd_iv) + (beta_mod * mean_mod) + 
    beta_int * ((mean_iv+sd_iv) * mean_mod) + intercept
  
  #hi modsd and low ivsd (+1SD)
  hilow <- beta_iv * (mean_iv-sd_iv) + beta_mod * (mean_mod+sd_mod) + 
    beta_int * ((mean_iv-sd_iv) * (mean_mod+sd_mod)) + intercept
  
  #hi modsd and med ivsd (+1SD)
  himed <- (beta_iv * mean_iv) + beta_mod * (mean_mod+sd_mod)  + 
    beta_int * (mean_iv * ((sd_mod-mean_mod)*-1)) + intercept
  
  #hi modsd and hi ivsd (+1SD)
  hihi <- beta_iv * (mean_iv+sd_iv) + beta_mod * (mean_mod+sd_mod) + 
    beta_int * ((mean_iv+sd_iv) * (mean_mod+sd_mod)) + intercept
  
  #build 2 matrices and then multiply them, to adjust the plot values for covars
  #to get estimates based on setting covariates to their sample means
  #see also http://www.calcul.com/show/calculator/matrix-multiplication
  if (iv_dicot < 1) {
    plotmat <- cbind(lowlow, lowmed, lowhi, medlow, medmed, medhi, hilow, himed, hihi)
  } else {
    plotmat <- cbind(lowlow, lowmed, medlow, medmed, hilow, himed)
  }
  
  if (use_covars > 0) {
    betamat <- matrix(rbind(intercept, beta_mod, beta_iv, beta_int, beta_Z1, beta_Z2, beta_Z3), ncol = 1)
  } else {
    betamat <- matrix(rbind(intercept, beta_mod, beta_iv, beta_int), ncol = 1)
  }
  
  #create and fill matrix for plot, +/1 1SD of moderator  (low/mean/high)
  if (meancent > 0) {
    modvals <- c(sd(noout$MOD)*-1,0,sd(noout$MOD))
    if (iv_dicot < 1) {
      ivvals <- c(sd(noout$IV)*-1,0,sd(noout$IV))
    } else {
      ivvals <- c(min(noout$IV),max(noout$IV))
    }
  } else {
    modvals <- c((mean(noout$MOD)-sd(noout$MOD)),mean(noout$MOD),(mean(noout$MOD)+sd(noout$MOD)))
    if (iv_dicot < 1) {
      ivvals <- c((mean(noout$IV)-sd(noout$IV)),mean(noout$IV),(mean(noout$IV)+sd(noout$IV)))
    } else {
      ivvals <- c(min(noout$IV),max(noout$IV))
    }
  }
  
  ivmat <- matrix(ivvals, nrow=ncol(plotmat))
  modmat <- matrix(modvals, nrow=ncol(plotmat))
  modmat <- matrix(modmat[order(modmat[,1]),], nrow=ncol(plotmat))
  
  #add re-calculated coeffs to plotdata matrix
  plotmat1 <- matrix(0, nrow=ncol(plotmat), ncol=nrow(betamat))
  for (i in 1:ncol(plotmat)){
    if (use_covars > 0) {
      plotmat1[i,] <- c(1, modmat[i], ivmat[i], (modmat[i]*ivmat[i]), means_Z1, means_Z2, means_Z3)
    } else {
      plotmat1[i,] <- c(1, modmat[i], ivmat[i], (modmat[i]*ivmat[i]))
    }
  }
  
  coeff_cov_adj <- (plotmat1 %*% betamat)
  plotdata <- cbind(num_dec(plotmat1[,3]), num_dec(plotmat1[,2]), num_dec(coeff_cov_adj))
  colnames(plotdata) <- c("iv_sd","mod_sd","dv")
  if (iv_dicot < 1) {
    rownames(plotdata) <- c("-1SD","-1SD","-1SD","Mean","Mean","Mean","+1SD","+1SD","+1SD")
  } else {
    rownames(plotdata) <- c("-1SD","+1SD","-1SD","+1SD","-1SD","+1SD")
  }
  
  ############### Calculate simples slopes ###################
  #Effect of Low SS = beta IV + beta IV:MOD * -1SD
  lowss <- (beta_iv + beta_int * ((sd_mod-mean_mod)*-1))
  medss <- (beta_iv + beta_int * (mean_mod))
  hiss <- (beta_iv + beta_int * (sd_mod+mean_mod))
  ssmat <- rbind(lowss,medss,hiss)
  
  #calculate standard error for simple slopes (this is were hell begins and never ends!!)
  if (use_covars > 0) {
    xy <- cbind(1, noout$MOD, noout$IV, (noout$MOD*noout$IV), noout$Z1, noout$Z2, noout$Z3)
  } else {
    xy <- cbind(1, noout$MOD, noout$IV, (noout$MOD*noout$IV))
  }
  invxtx <- solve(t(xy) %*% xy)
  coeff <- invxtx %*% t(xy) %*% noout$DV
  k3 <- nrow(coeff)
  resid <- data.matrix(noout$DV) - xy %*% coeff
  sse <- sum(resid ^ 2)
  mse <- sse/(nrow(noout)-ncol(xy))
  h <- xy[,1, drop=F]
  ss_output <- matrix(unique(modmat),nrow(unique(modmat)),7)
  
  #this is were magic happens, hc3 adjustment
  if (heteros > 0) {
    for (i in 1:nrow(noout)) {
      h[i,1] <- xy[i,, drop=F] %*% invxtx %*% t(xy[i,, drop=F]);
    }
    for (i in 1:k3) {
      xy[,i] <- (resid[,ncol(resid), drop=F]/(1-h))*xy[,i, drop=F]
    }
  } else {
    for (i in 1:k3) {
      xy[,i] <- sqrt(mse)*xy[,i, drop=F];
    }
  }
  
  covmat <- invxtx %*% t(xy) %*% xy %*% invxtx
  degfree <- nrow(noout)-ncol(xy)
  qtest <- qt(.975,degfree)
  
  for (m in 1:nrow(unique(modmat))) {
    #final standard error matrix
    temp1 <- 1 * covmat[3,3]
    temp2 <- (unique(modmat)[m,1]^2) * covmat[4,4]
    temp3 <- (unique(modmat)[m,1]*2) * covmat[3,4]
    #ss_output
    ss_cfse <- sqrt(temp1+temp2+temp3)
    #t-test = beta or effect of ss / standard error
    ss_ttest <- (ssmat[m,1] / ss_cfse)
    #p-value
    ss_pval <- 2*(1-pt(abs(ss_ttest),(nrow(noout)-ncol(xy))));
    #LLCI and ULCI
    ss_loci <- ssmat[m,1]-qtest*ss_cfse;
    ss_upci <- ssmat[m,1]+qtest*ss_cfse;
    ss_output[m,] <- cbind(
      num_dec(unique(modmat)[m,1]), #1sd change in predictor
      num_apa(ssmat[m,1]), #effect, not beta, but B, as it is unstandardized
      num_apa(ss_cfse), #standard error
      num_dec(ss_ttest), #t-value
      add_stars(num_dec(ss_pval)),
      num_dec(ss_loci), #lower confidence interval
      num_dec(ss_upci) #upper confidence interval
    )
  }
  
  ss_output <- as.data.frame(ss_output)
  colnames(ss_output) <- c(mod_var, "Effect", "se", "t", "p    ", "LLCI", "ULCI")
  rownames(ss_output) <- c("-1SD","Mean","+1SD")
  
  #t-test for simple slopes
  #for comparision of ss function which does not adjust for hc3
  #library(pequod)
  #test <- lmres(DV ~ IV*MOD+Z1+Z2+Z3, data=noout)
  #simpleSlope(test, pred="IV", mod1="MOD")
  
  ########### Moderation output #############
  sink(paste(filename, "txt", sep="."))
  
  star_head <- function(string, breaks=1, char="*") {
    total_stars = 84
    len_str = nchar(string)+2
    num_stars = round((total_stars - len_str)/2,0)
    l_stars = strrep(char,num_stars)
    str_breaks = strrep("\n",breaks)
    r_stars = l_stars
    str_len_total = (2 * nchar(r_stars)) + len_str
    if (str_len_total < total_stars) {
      r_stars = strrep(char,num_stars+1)
    } else if (str_len_total > total_stars) {
      r_stars = strrep(char,num_stars-1)
    } else {
      r_stars = strrep(char,num_stars)
    }
    cat(str_breaks,l_stars," ",toupper(string)," ",r_stars,str_breaks,sep="")
  }
  
  # format variables for output (add square brackets)
  iv_var_var <- paste("[",iv_var,"]", sep="")
  mod_var_var <- paste("[",mod_var,"]", sep="")
  dv_var_var <- paste("[",dv_var,"]", sep="")
  if (use_covars > 0) {
    cov1_var <- paste("[",cov1,"]", sep="")
    cov2_var <- paste("[",cov2,"]", sep="")
    cov3_var <- paste("[",cov3,"]", sep="")
  }
  
  star_head("Model 1 - Moderation",0," ")
  star_head("Variables",2)
  cat("Dependent variable   (DV)  = ", dv_var_desc, dv_var_var, "\n")
  cat("Independent variable (IV)  = ", iv_var_desc, iv_var_var,"\n")
  cat("Moderating variable  (MOD) = ", mod_var_desc, mod_var_var)
  if (use_covars > 0) {
    cat("\n", "Control variable 1   (Z1)  =  ", cov1_desc, cov1_var, "\n", sep="")
    cat("Control variable 2   (Z2)  = ", cov2_desc, cov2_var, "\n")
    cat("Control variable 3   (Z3)  = ", cov3_desc, cov3_var)
  }
  
  # overview of descriptives (insanity check)
  star_head("Descriptives",2)
  print(describe(master)) #psych lib
  
  star_head("Outliers and Data screening")
  if (sum(outlier_id) > 0) {
    cat("\n",sum(totalout == 3, na.rm = TRUE), " outlier(s) found and removed\n", sep="")
    cat("Case ID(s) removed:", outlier_id, "\n")
    cat("N =", nrow(noout),", after outlier removal\n\n")
  } else {
    cat("\n(+) No outliers found\n")
  }
  
  #check that after outlier removal skew and kurtosis are acceptable
  #A normal distribution has a skewness of 0. Skew to the left, e.g. when the 
  #hump is on the RIGHT side, has a negative negative value.
  #Kurtosis is a measure of the heaviness of the tails of a distribution. 
  #A normal distribution has kurtosis 0. Extremely nonnormal distributions 
  #may have high positive or negative kurtosis values, while nearly normal 
  #distributions will have kurtosis values close to 0. 
  desc_noout <- describe(noout)
  ku_var <- which(abs(desc_noout$kurtosis) > 2)
  sk_var <- which(abs(desc_noout$skew) > 2)
  
  if (length(ku_var) > 0) {
    cat("(-) Kurtosis in: ")
    ku_data <- num_dec(desc_noout[ku_var,]$kurtosis)
    ku_out <- paste(rownames(desc_noout[ku_var,]), ku_data, sep=" = ", collapse = ", ")
    write.matrix(ku_out)
  }
  if (length(sk_var) > 0) {
    cat("(-) Skewness in: ")
    sk_data <- num_dec(desc_noout[sk_var,]$skew)
    sk_out <- paste(rownames(desc_noout[sk_var,]), sk_data, sep=" = ", collapse = ", ")
    write.matrix(sk_out)
  }
  
  #When the variance of errors differs, at different values of the IV,
  #heteroscedasticitz is inidicated. If extreme it can lead to serious
  #distortion of findings, becasue posibility of Type I error incerases.
  #The scatterplot should be examined for heteroskedasticity, however
  #this little test will also help to determine if this is the case.
  #The studentized Breusch-Pagan test was proposed by R. Koenker in his 
  #1981 article "A Note on Studentizing a Test for Heteroscedasticity"
  #If the test is significant, then data is heteroscedastic, which is bad!
  #Read more here https://analysights.wordpress.com/tag/goldfeld-quandt-test/
  if (bp$p.value < 0.1) {
    cat("(-) Data shows hereroskedasticity, p =",add_stars(num_dec(bp$p.value)))
  } else {
    cat("(+) No hereroskedasticity, p =",num_apa(bp$p.value))
  }
  
  if (length(ku_var) | length(sk_var)>0) {
    #Trochim & Donnelly, 2006; Field, 2000 & 2009; Gravetter & Wallnau, 2014
    cat("\n\n* Skew and Kurtosis values +/-2 are considered acceptable.")
  }
  
  #print test of overall model
  star_head("Moderation analyses",2)
  #this gets the new model summary
  # use Heteroscedasticity-consistent SEs
  if (heteros > 0) {
    wald <- waldtest(model, vcov = vcovHC(model, "HC3"))
    model_out <- coeftest(model, vcov = vcovHC(model, "HC3"))
  } else {
    wald <- waldtest(model)
    model_out <- coeftest(model)
  }
  
  wald_result <- cbind(
    num_apa(sqrt(summary(model)$r.squared)), #R
    num_apa(summary(model)$r.squared), #R squared
    num_apa(mse), #mean square error
    num_dec(wald$F[2]), #f-value
    abs(wald$Df[2]), #degrees of freedom 1
    model$df.residual, #degrees of freedom 2
    add_stars(num_dec(wald$`Pr(>F)`[2]))
  )
  cat("Model Summary\n")
  colnames(wald_result) <- c("R", "R-sq", "MSE", "F", "df1", "df2", "p    ")
  rownames(wald_result) <- ""
  print.table(wald_result, right=T, quote=F, print.gap=5)
  
  cat("\nModel (coeff = beta)\n")
  model_out.coeff <- model_out[,1]
  model_out.se <- model_out[,2]
  model_out.tval <- model_out[,3]
  model_out.pval <- format.pval(model_out[,4], digits=2, eps=0.001)
  model_out.llci <- model_out.coeff-qtest*model_out.se
  model_out.ulci <- model_out.coeff+qtest*model_out.se
  model_out.df <- cbind(num_dec(model_out.coeff),
                        num_apa(model_out.se),
                        num_dec(model_out.tval),
                        add_stars(num_dec(model_out[,4])),
                        num_dec(model_out.llci),
                        num_dec(model_out.ulci)
  )
  model_out.df <- as.data.frame(model_out.df)
  if (use_covars > 0) {
    rownames(model_out.df) <- c("constant",iv_var,mod_var,cov1_var,cov2_var,cov3_var,"int1")
  } else {
    rownames(model_out.df) <- c("constant",iv_var,mod_var,"int1")
  }
  colnames(model_out.df) <- c("coeff", "se", "t", "p    ", "LLCI", "ULCI")
  print(model_out.df, quote=F, right=T, print.gap=print_gap)
  
  cat("\nProduct term key:\n")
  cat("int1   ",iv_var,"   x   ",mod_var,"\n")
  
  #Calculate r-squared change due to interactions
  #Read here http://goo.gl/HI4VnQ
  mnv <- sum(noout$DV/nrow(noout))
  ssty <- sum((noout$DV-mnv) ^ 2)
  r2full <- 1-(sse/ssty)
  model_int_tval = tail(model_out.tval, n=1)
  r2_chg <- ((model_int_tval * model_int_tval)*(1-r2full))/degfree
  int_fval <- (model_int_tval * model_int_tval)
  int_pval <- tail(model_out.pval, n=1)
  cat("\nR-square increase due to interaction(s):\n")
  chg_result <- cbind(num_dec(r2_chg),
                      num_dec(int_fval),
                      1,
                      model$df.residual,
                      add_stars(int_pval)
  )
  colnames(chg_result) <- c("R2-chng", "F", "df1", "df2", "p    ")
  rownames(chg_result) <- "int1"
  print(chg_result, quote=F, right=T, print.gap=print_gap)
  
  star_head("Simple slopes",1)
  cat("\nConditional effect of X on Y at values of the moderator(s):\n\n")
  print(ss_output, right=T, quote=F, print.gap=(print_gap-1))
  
  cat("\nEffect is not Beta, but B, as it is unstandardized here.\n")
  cat("Values for quantitative moderators are the mean and plus/minus one SD from mean.") 
  
  star_head("Plotting data",2)
  cat("Data for visualising conditional effext of X on Y\n\n")
  print.table(plotdata, right=T, quote=F, print.gap=print_gap)
  cat("\n* Estimates are based on setting covariates to their sample means.")
  
  ################ Johnson-Neyman Technique #################
  if (jntech > 0) {
    
    if (as.numeric(int_pval) <= 0.1) { #only do if int is sig.
      jn_results = matrix()
      qtest2 <- (qtest*qtest)
      ajn <- (qtest2*covmat[4,4])-(coeff[4,1]*coeff[4,1])
      bjn <- 2*((qtest2*covmat[3,4])-(coeff[3,1]*coeff[4,1]))
      cjn <- (qtest2*covmat[3,3])-(coeff[3,1]*coeff[3,1])
      radarg <- (bjn*bjn)-(4*ajn*cjn)
      den <- 2*ajn
      nrts <- 0
      
      if ((radarg >= 0) && (den != 0)) {
        x21 <- (-bjn+sqrt(radarg))/den;
        x22 <- (-bjn-sqrt(radarg))/den;
        roots <- 0;
        if ((x21 >= min(noout$MOD)) && (x21 <= max(noout$MOD)))  {
          nrts <- 1
          roots <- rbind(roots,x21)
        }
        if ((x22 >= min(noout$MOD)) && (x22 <= max(noout$MOD)))  {
          nrts <- nrts+1
          roots <- rbind(roots,x22)
        }
        roots <- cbind(as.matrix(roots),matrix(0,nrow(as.matrix(roots)),2))
        
        if (nrts > 0)  {
          roots <- roots[2:nrow(roots), 1:3, drop=F] #value
          rootsum <- (noout$MOD < roots[1,1]) #true/false
          roots[1,2] <- (sum(rootsum) / nrow(noout))*100 #% below
          rootsum <- (noout$MOD > roots[1,1])
          roots[1,3] <- (sum(rootsum) / nrow(noout))*100 #% above
          
          if (nrow(roots) == 2) {
            rootsum <- (noout$MOD < roots[2,1]);
            roots[2,2] <- (sum(rootsum) / nrow(noout))*100 #% below
            rootsum <- (noout$MOD > roots[2,1]);
            roots[2,3] <- (sum(rootsum) / nrow(noout))*100 #% above
          }
          
          #calculate conditional effect of x on y at values of moderator
          jnvals <- matrix(0,(21+nrts),7); #define empty matrix, 22 rows
          
          #fill first column of jn matrix (moderator values)
          for (i in 0:20) {
            jnvals[(i+1),1] <- min(noout$MOD) + (i*((max(noout$MOD)-min(noout$MOD))/20));
          }
          
          for (i in 1:nrts) {
            for (j in 2:nrow(jnvals)) {
              if ((roots[i,1, drop=F] > jnvals[(j-1),1, drop=F]) && 
                  (roots[i,1, drop=F] < jnvals[j,1, drop=F])) {
                jnvals[(j+1):(21+i),1] <- jnvals[j:(20+i),1, drop=F]
                jnvals[j,1] <- roots[i,1, drop=F]
              }
            }
          }
          
          #fill rest of jn matrix
          for (i in 1 : nrow(jnvals)) {
            jnvals[i,2] <- coeff[3,1] + coeff[4,1] * jnvals[i,1, drop=F] #effect
            jnvals[i,3] <- sqrt(covmat[3,3] + 2 * jnvals[i,1, drop=F] * covmat[3,4] + 
                                  (jnvals[i,1, drop=F]^2) * covmat[4,4]) #se
            jnvals[i,4] <- jnvals[i,2, drop=F] / jnvals[i,3, drop=F] #t-val
            if (dv_dicot > 0)  {
              jnvals[i,5] <- 2 * (1-pnorm(abs(jnvals[i,4, drop=F]), degfree)) #p-val
            } else {
              jnvals[i,5] <- 2 * (1-pt(abs(jnvals[i,4, drop=F]), degfree)) #p-val
            }
            jnvals[i,6] <- jnvals[i,2, drop=F] - qtest * jnvals[i,3, drop=F] #llci
            jnvals[i,7] <- jnvals[i,2, drop=F] + qtest * jnvals[i,3, drop=F] #ulci
          }
        }
        
        star_head("Johnson-Neyman Technique",2)
        cat("Moderator values(s) defining Johnson-Neyman significance region(s)\n")
        if (nrow(roots) == 2) {
          root_row1 <- cbind(num_dec(roots[1,1]),num_dec(roots[1,2]),num_dec(roots[1,3]))
          root_row2 <- cbind(num_dec(roots[2,1]),num_dec(roots[2,2]),num_dec(roots[2,3]))
          root_results <- rbind(root_row1, root_row2)
          rownames(root_results) <- c("","")
        } else {
          root_results <- cbind(num_dec(roots[1,1]),num_dec(roots[1,2]),num_dec(roots[1,3]))
          rownames(root_results) <- c("")
        }
        
        colnames(root_results) <- c("Value","% below","% above")
        print(root_results, quote=F, right=T, print.gap=print_gap);
        
        cat("\nConditional effect of X on Y at values of the moderator (M)\n")
        jn_results <- cbind(num_dec(jnvals[,1]), #mod
                            num_dec(jnvals[,2]), #effect
                            num_dec(jnvals[,3]), #se
                            num_dec(jnvals[,4]), #t
                            add_stars(num_dec(jnvals[,5])), #p
                            num_dec(jnvals[,6]), #llci
                            num_dec(jnvals[,7]) #ulci
        )
        
        jn_results <- as.data.frame(jn_results)
        rownames(jn_results) <- c()
        colnames(jn_results) <- c(mod_var,"Effect","se","t","p    ","LLCI","ULCI")
        print(jn_results, quote=F, right=T, print.gap=print_gap, row.names=F)
      }
      
      #call johnson-neyman analysis, if interaction was signficant
      if (!is.na(jn_results[1,1])==TRUE) { #no NA
        img_save(2, img_dim)
        #jnplot()
        if (!exists("roots[2,1]")) {
          jn_root1 <- roots[1,1]
          jn_root2 <- NA
        } else {
          jn_root1 <- roots[1,1]
          jn_root2 <- roots[2,1]
        }
        suppressWarnings(custom_jnplot(jn_root1, jn_root2 ,jnvals))
        dev.off()
      } 
    } else {
      star_head("Johnson-Neyman Technique",2)
      cat("There are no statistical significance transition points within the observed\nrange of the moderator.\n")
    }
  }
  
  ################ notes ###############
  star_head("Analysis Notes")
  cat("\nSignficance levels: . p < 0.1, * p < 0.5, ** p < 0.1, *** p < 0.001.\n")
  cat("Level of confidence for all confidence intervals in output: 95.00.\n")
  if (meancent > 0) {
    cat("Predictor and moderator variables were mean centered prior to analysis.\n")
  }
  if(heteros > 0) {
    cat("All standard errors for continuous outcome models are based on the HC3 estimator.\n")
  }
  
  sink() #return to output to terminal
  
  ################ plot moderation function #################
  
  
  #save moderation plot
  img_height <- (img_dim-(img_dim*0.265))
  img_save(1,img_height)
  modplot(plotdata)
  dev.off()
  
} else {
  stop("Cant find the data file!")
}
#reset parameters
options(scipen=0)