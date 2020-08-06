#
# Clear the workspace
#
rm(list = ls())
gc()

#
# Close graphic devices
#
graphics.off()

#
# Stop automatic conversion of string to factors
#
options(stringsAsFactors = FALSE)

#
# Function for Bland & Altman style plot
#
ba.plot <- function(a, b, main, xlab, ylab, xlim, ylim)
  {
  true <- b; diff.two <- a - b; mean.diff <- mean(diff.two)
  plot(true, diff.two, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, pch = ".", col = "gray", cex.lab = 0.8, cex.axis = 0.8, frame.plot = FALSE)
  x <- data.frame(cbind(true, diff.two))
  LLA <- ULA <- NULL
  i <- unique(x$true)
  for(j in i)
    {
    y <- subset(x, true == j)
    LLA <- c(LLA, quantile(y$diff.two, probs = 0.025, names = FALSE))
    ULA <- c(ULA, quantile(y$diff.two, probs = 0.975, names = FALSE))
    }
  lines(lowess(i, LLA, f = 1), lty = 3)
  lines(lowess(i, ULA, f = 1), lty = 3)
  lines(x = c(0, max(true)), y = c(mean.diff, mean.diff), lty = 1)
  }

#
# Function to estimate prevalence by PROBIT method (2 variables)
#
#                                               +----------+
# +---------------------------------------------+  Sample  +-----+
# |                                             +----------+     |
# |                                                              |
# |             .-----------.          .-----------.             |
# |         _.-'             `--.  _.-'             `--.         |
# |       ,'                     ,'                     `.       |
# |     ,'                     ,'  `.                     `.     |
# |    /                      /      \                      \    |
# |   /                      /        \                      \   |
# |  ;                      ;          :                      :  |
# |  |   (A) MUAC < 115     |    X     |    (C) OEDEMA        |  |
# |  |                      |          |                      |  |
# |  :                      :          ;                      ;  |
# |   \                      \        /                      /   |
# |    \                      \      /                      /    |
# |     \                      \    /                      /     |
# |      `.                     `.,'                     ,'      |
# |        `.                   ,'`.                   ,'        |
# |          `--.           _.-'    `--.           _.-'          |
# |              `---------'            `---------'              |
# |                                                              |
# +--------------------------------------------------------------+
#
#
#
#                           P = A + B - X
#
#
prevalencePROBIT2 <- function(data, replicates = 39)
  {
  analysis <<- "PROBIT.SAM2"
  p <- NULL
  for(r in 1:replicates)
  {
  rdata <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
  ## MUAC < 115
  x <- median(rdata$muac)
  s <- IQR(rdata$muac) / 1.34898
  pMUAC   <- pnorm(q = 125, mean = x, sd = s)
  ## Prevalence of Oedema
  pOEDEMA <- 0 ## default to non oedema
  if(!is.na(table(rdata$oedema)["1"])) ## There is oedema!
    {
    pOEDEMA <- prop.table(table(rdata$oedema))["1"]
    }
  ## Intersection between MUAC < 115 and oedema
  if(any(dim(table(rdata$muac < 115, rdata$oedema))) < 2) ## It *not* is a 2-by-2 table
  	{
  	pINTERSECTION <- rep(0, 3)
  	} else {
    temp <- prop.test(table(rdata$muac < 115, rdata$oedema)[2,1], sum(table(rdata$muac < 115, rdata$oedema)))
  	pINTERSECTION <- temp$estimate
  	}
  p <- c(p, pMUAC + pOEDEMA - pINTERSECTION)
  }
  return(list(p = median(p) * 100 , lci = quantile(p, probs = 0.025) * 100, uci = quantile(p, probs = 0.975) * 100))
  }

#
# Function to estimate prevalence by PROBIT method (3 variables)
#
#                                                +----------+
# +----------------------------------------------+  Sample  +----+
# |                                              +----------+    |
# |                                                              |
# |                                                              |
# |              .-----------.          .-----------.            |
# |          _.-'             `--.  _.-'             `--.        |
# |        ,'                     ,'                     `.      |
# |      ,'                     ,'  `.                     `.    |
# |     /                      /      \                      \   |
# |    /                      /        \                      \  |
# |   ;                      ;    X     :                      : |
# |   |   (A) MUAC < 115     |          |     (B) WHZ < -3     | |
# |   |                 .    |          |                      | |
# |   :                      .-----------.                     ; |
# |    \                 _.-' \        /  `--.                /  |
# |     \              ,'      \  Z   /       `.             /   |
# |      \           ,'         \    /          `.          /    |
# |       `.        /     W      `.,'      Y      \       ,'     |
# |         `.     /             ,'`.              \    ,'       |
# |           `--.;          _.-'    `--.           :.-'         |
# |               |---------'            `---------'|            |
# |               |                                 |            |
# |               :                                 ;            |
# |                \          (C) OEDEMA           /             |
# |                 \                             /              |
# |                  \                           /               |
# |                   `.                       ,'                |
# |                     `.                   ,'                  |
# |                       `--.           _.-'                    |
# |                           `---------'                        |
# |                                                              |
# +--------------------------------------------------------------+
#
#
#
#            P = A + B + C - X - Y - Z
#
#            where:
#
#                A = P(MUAC < 115 mm)
#                B = P(WHZ < -3)
#                C = P(OEDEMA)
#                W = P(MUAC < 115 mm & OEDEMA)
#                X = P(MUAC < 115 mm & WHZ = -3)
#                Y = P(WHZ < -3 & OEDEMA)
#                Z = P(MUAC < 115 & WHZ < -3 & OEDEMA)
#
prevalencePROBIT <- function(data, replicates = 39)
  {
  analysis <<- "PROBIT.GAM3"
  N <- nrow(data)
  p <- NULL
  for(r in 1:replicates)
  {
  rdata <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
  ## MUAC < 125
  x <- median(rdata$muac)
  s <- IQR(rdata$muac) / 1.34898
  A <- pnorm(q = 125, mean = x, sd = s)
  ## WHZ < -2
  x <- median(rdata$whz)
  s <- IQR(rdata$whz) / 1.34898
  B <- pnorm(q = -2, mean = x, sd = s)
  ### Now add in oedema
  C <- sum(data$oedema == 1) / N
  ## W, X, Y, Z
  W <- sum(data$muac < 125 & data$oedema == 1) / N
	X <- sum(data$muac >= 125 & data$whz < -2 ) / N
	Y <- sum(data$whz < -2 & data$oedema == 1)  / N
	Z <- sum(data$muac < 125 & data$whz < -2 & data$oedema == 1)  / N
  ## Prevalence
  prevalence <- A + B + C - W - X - Y - Z
  p <- c(p, prevalence)
  }
  return(list(p = median(p) * 100 , lci = quantile(p, probs = 0.025) * 100, uci = quantile(p, probs = 0.975) * 100))
  }

#
# Function to get prevalence by FREQUENCY method (use for "populations" only)
#
prevalenceFREQUENCY <- function(population)
  {
  p <- prop.table(table(population$case))["1"] * 100
  p <- ifelse(is.na(p), 0, p)
  names(p) <- NULL
  return(p)
  }

#
# Function to sample from a population and estimate prevalence
#
samplePopulation <- function(name = "Not Specified", population, size = 384)
  {
  true <- prevalenceFREQUENCY(population)
  sampled.rows <- sample(x = 1:nrow(population), size = size, replace = FALSE)
  data <- population[sampled.rows, ]
  estimate <- prevalencePROBIT(data)
  result <- list(name = name, n = size, true = true, estimate = estimate$p, lci = estimate$lci, uci = estimate$uci)
  return(result)
  }

#
# Simulation parameters
#
POP.SIZE <- 17000
SAMPLE.SIZES <- 384
REPLICATES <- 5

#
# Run sampling simulation
#
#files.2.process <- dir(path = ".", pattern = "\\.csv$")
files.2.process <- dir(path = "probitGAM", pattern = "\\.csv$")
N <- length(files.2.process) * length(SAMPLE.SIZES) * REPLICATES
surveyResults <- data.frame(name = character(N), n = numeric(N), true = numeric(N), estimate = numeric(N), lci = numeric(N), uci = numeric(N))
rowIndex <- 0
for(file.name in files.2.process)
  {
  #
  # Report name of current population and retrieve data
  #
  cat("\n", file.name, " : ", sep = "")
  #x <- read.table(file.name, header = TRUE, sep = ",")
  x <- read.table(paste("probitGAM", file.name, sep = "/"), header = TRUE, sep = ",")
  #
  # Make a population with N = c. POP.SIZE
  #
   multiples <- ceiling(POP.SIZE / nrow(x))
   y <- data.frame()
   for(i in 1:multiples)
     {
     y <- rbind(y, x)
     }
  #
  # Apply case-defintions to population (change as required)
  #
  ## y$case <- ifelse(y$muac < 115 | y$oedema == 1, 1, 2)
  ## y$case <- ifelse(y$muac < 125 | y$oedema == 1, 1, 2)
  ##y$case <- ifelse(y$muac < 115 | y$whz < -3 | y$oedema == 1, 1, 2)
  y$case <- ifelse(y$muac < 125 | y$whz < -2 | y$oedema == 1, 1, 2)
  #
  # Simulate surveys with different sample sizes
  #
  for(n in SAMPLE.SIZES)
    {
    cat(n,  " ", sep = "")
    for(replicates in 1:REPLICATES)
      {
      cat(".")
      rowIndex <- rowIndex + 1
      surveyResults[rowIndex, ] <- samplePopulation(name = file.name, population = y, size = n)
      }
    }
  }

#
# Bland & Altman style plots
#
quartz(height = 6, width = 6, pointsize = 9)
par(pty = "s")
for(sampleSize in unique(surveyResults$n))
  {
  y <- subset(surveyResults, n == sampleSize)
  ba.plot(y$estimate,
          y$true,
          main = paste("n = ", sampleSize, sep = ""),
          xlab = "True prevalence (%)",
          ylab = "Estimated - true prevalence (%)",
          xlim = c(0, 50),
          ylim = c(-15, 15))
  }
dev2bitmap(file = paste(analysis, ".BA.png", sep =""), type = "pngmono", height = 6, width = 6, pointsize = 9, res = 600)

#
# Relative precision
#
quartz(height = 6, width = 6, pointsize = 9)
par(pty = "s")
diff2uci <- surveyResults$uci - surveyResults$estimate
diff2lci <- surveyResults$estimate - surveyResults$lci
meanDiff <- (diff2uci + diff2lci) / 2
relativePrecision <- (meanDiff / surveyResults$estimate) * 100
boxplot(relativePrecision ~ surveyResults$n,
        xlab = "Sample Size",
        ylab = "Relative Precision (%)",
        outline = FALSE,
        frame.plot = FALSE)
abline(h = 30, lty = 3)
dev2bitmap(file = paste(analysis, ".BP.png", sep =""), type = "pngmono", height = 6, width = 6, pointsize = 9, res = 600)

#
# Save workspace for later use ...
#
save.image(file = paste(analysis, ".RData", sep =""))

