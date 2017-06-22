##  AdvImgMetricsSurvival.R - takes input of multiple advanced imaging metrics, uses
##    lasso to determine which metrics have the strongest relationship with survival.
##  Author: S. Hilz
##  Date: 2017.06.22

## Functions
buildDataStructure <- function(img, surv){ 
  # merges imaging with survival data; determines if conflicts exists
  colnames(surv)[grep('sfnumber',colnames(surv))] <- 'SF'
  merged <- merge(img, surv, by='SF')
  if(any(duplicated(merged$VIAL))==TRUE){
    print('WARNING: Duplicate entries with conflicting data detected in your input data. These entries will be dropped. Please manually resolve and re-read in data to include')
    duplicates <- merged$VIAL[duplicated(merged$VIAL)]
    print(merged[which(merged$VIAL %in% duplicates),])
  }
  merged <- merged[which(!merged$VIAL %in% duplicates),]
  return(merged)
}

buildAnalysisStructure <- function(data, dependent, independent){
  rownames(data) <- data$VIAL
  if (!dependent %in% c('osdaysspore','TTP','osmosfirstsx','osmossporesx','TTPmos')){
    print('ERROR: Response variable is not a valid option. Please check spelling and column headers in input data')
  } else {
    data <- data[,c(independent, dependent)]
    data <- data[complete.cases(data),]
    y <- as.vector(data[,dependent])
    x <- as.matrix(data[,independent])
    toReturn <- list('x'=x,'y'=y)
    return(toReturn)
  }
}

optimizeVariables <- function(data, dependent){
  rownames(data) <- data$VIAL
  if (!dependent %in% c('osdaysspore','TTP','osmosfirstsx','osmossporesx','TTPmos')){
    print('ERROR: Response variable is not a valid option. Please check spelling and column headers in input data')
  } else {
    y <- as.vector(data[,dependent])
    keep <- which(!is.na(y))
    y <- y[keep]
    x <- as.matrix(data[keep,8:33])
    toReturn <- c()
    for (i in 1:dim(x)[2]){
      toReturn <- rbind(toReturn, c(colnames(x)[i], sum(!is.na(x[,i]))))
    }
    toReturn <- as.data.frame(toReturn, rownames=FALSE)
    colnames(toReturn) <- c('metric','cases')
    toReturn$metric <- as.character(toReturn$metric) 
    toReturn$cases <- as.numeric(as.character(toReturn$cases))
    return(toReturn)
  }
}

lbs_fun <- function(fit, ...) {
  L <- length(fit$lambda)
  x <- log(fit$lambda[L])
  y <- fit$beta[, L]
  labs <- names(y)
  text(x, y, labels=labs, ...)
}

## Libraries
if("glmnet" %in% rownames(installed.packages()) == FALSE) {
  install.packages("dplyr", repos="http://cran.rstudio.com/")
}

library(glmnet)

## Input files
imagingMetricsFile <- 'SPORE_CPMG_FINAL_DATA_STEPHANIE.csv'
survivalFile <- 'SPORE_DATA_022715_Full_Data.csv'

#################### MAIN ####################

# bring in imaging data
imagingData <- read.table(imagingMetricsFile,sep=',',header=TRUE, stringsAsFactors=FALSE)

# bring in survival data
survivalData <- read.table(survivalFile,sep=',',header=TRUE, stringsAsFactors=FALSE)

# build joint data structure, which gives survival (at level of SF#) for each set of imaging metrics (at level of vial ID)
mergedData <- buildDataStructure(imagingData, survivalData)

# get stats on number of complete cases for all combs independent variable input based on y
# second argument is response varialbe (options are 'osdaysspore','TTP','osmosfirstsx','osmossporesx','TTPmos')
casesByMetric <- optimizeVariables(mergedData, 'TTP')

#pick independent variables by cutoff number
cutoff <- 80
independentVariables <- casesByMetric[which(casesByMetric$cases > cutoff),]$metric

# build compatible variable matrix, x, and response vector, y 
# second argument is response varialbe (options are 'osdaysspore','TTP','osmosfirstsx','osmossporesx','TTPmos')
# third argument is vector of independent varialbes (options are "X2HG", "IDH.IHC.SEQ", "X2HG...IDH", "Cho", "PC", "GPC", "ETH", "PE","GSH","TAU","H.TAU","NAA","ASP", "GLU", "GLN", "GLC", "MI", "SI", "GABA", "CRE", "LAC", "ALA", "GLY", "VAL", "THR", "ACE", "X2HG.1", "BET", "SUC") 
toFit <- buildAnalysisStructure(mergedData, 'TTP', independentVariables)

# perform lasso fit
fit = glmnet(toFit$x, toFit$y)

# plot fit
plot(fit, xvar="lambda")
lbs_fun(fit)

# do cross validation
cvfit = cv.glmnet(toFit$x, toFit$y)

# plot
plot(cvfit)

# check coefficients
coef(cvfit, s='lambda.1se')

