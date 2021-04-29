
#library(devtools)
#devtools::create("G:\\R packages\\predictiontest")

# Functions for Predictiontest package

# Function list:
# 1) Original exact/asymptotic prediction-test - DONE
      # printing function - DONE
# 2) Prediction bootstrap - DONE
# 3) Weights function - DONE
# 4) Prediction results - DONE
  #- Not done, error in gtvar vs group in group_by - dependent on specific naming
# 5) Plotting function
# 6) Add application
# 7) Make sure output from test works with summary()

# Package dependencies
library(tidyr)
library(dplyr)
library(robustbase)

# # Two different prediction test functions
# - one the original with exact and approximation
# - other the permutation version
# - prediction results function:
#   - given a data set and group or pre post designation and the variables and a set of up/down predictions
#  - return a vector of the results and weights and the test statistics
# - prediction plotting function
#   - 
# - O'Briens ols test?
# - Simulator for power calculations 


# Testing data sets ----------------
# data set 1 
# pre-post

id <- c(1:10,1:10)
time <- (c(rep(1,10), rep(2,10)))
v1 <- rnorm(20,0,1)
v2 <- rnorm(20,0,1)
v3 <- rnorm(20,0,1)
v4 <- rnorm(20,0,1)
v5 <- rnorm(20,0,1)
v6 <- rnorm(20,0,1)

data1 <- as.data.frame(cbind(id, time, v1, v2 ,v3, v4, v5,v6))


# data set 2 
# pre-post - different variable names and order

identifier <- c(1:10,1:10)
year <- (c(rep(1,10), rep(2,10)))
e1 <- rnorm(20,0,1)
e2 <- rnorm(20,0,1)
m3 <- rnorm(20,0,1)
endpoint <- rnorm(20,0,1)
check <- rnorm(20,0,1)
v6 <- rnorm(20,0,1)

data1 <- as.data.frame(cbind(identifier, v6, e1, e2 ,endpoint, m3, check,m3, year))


# data set 3 
# group

id <- 1:20
group <- c(rep(1,10), rep(2,10))
v1 <- rnorm(20,0,1)
v2 <- rnorm(20,0,1)
v3 <- rnorm(20,0,1)
v4 <- rnorm(20,0,1)
v5 <- rnorm(20,0,1)
v6 <- rnorm(20,0,1)

data3 <- as.data.frame(cbind(id, group, v1, v2 ,v3, v4, v5,v6))


# data set 4
# group - but with different variables and ordering

MRN <- 1:20
grp <- c(rep(1,10), rep(2,10))
bloodpress <- rnorm(20,0,1)
MoCA <- rnorm(20,0,1)
eGFR <- rnorm(20,0,1)
hyper <- rnorm(20,0,1)
x4 <- rnorm(20,0,1)
pmid4 <- rnorm(20,0,1)

data3 <- as.data.frame(cbind(MRN, bloodpress, MoCA, eGFR ,hyper, x4, grp, pmid4))





# Original prediction test ------------------------------------
# Inputs
  # Weights ( a vector)
  # Prediction results
  # Null phi default = 0.50
  # alpha, default = 0.05
  # exact  - true or false, default FALSE, and only available for m < 25

weights <- W
results <- prediction.results(data1, variables = variables,direction = "up", type = "prepost", gtvar = "time")[1]


prediction.test <- function(weights, results, nullphi = 0.50, alpha = 0.05, exact = TRUE){
  
  ntests <- length(weights)
  results <- unlist(results)
  teststat <- weights%*%results
  correct <- sum(results)
  
  if (exact == FALSE){
    type = 1
    z <- (teststat - nullphi*(sum(weights)) )/sqrt(nullphi*(1-nullphi)*sum(weights^2))
    pval <- 1-pnorm(z)
      if (pval < alpha)
        {
          decision <- 1
        } else if (pval >= alpha){
          decision <- 0
        }
  } else if( exact == TRUE & ntests < 25)
  {
    type = 2
    nperm <- 2^ntests
    perms <- as.matrix(expand.grid(rep(list(0:1), ntests)))
    values <- perms%*%as.matrix(weights)
    rank <- as.data.frame(cbind(values,rank(values)))
    pval <- 1-(rank[which(rank$V1 == as.numeric(teststat)),2]/nperm)
      if (pval < alpha)
        {
          decision <- 1
        } else if (pval >= alpha){
          decision <- 0
         }
    
  } else if (exact == TRUE & ntests >= 25)
  {
    type = 2
    stop("The exact test is only available for 25 or fewer endpoints due to the large number of permutations of the test statistic.")
  }
  
  x <- list(statistic = teststat, p.value = pval, null.value = nullphi, m = ntests, alpha = alpha, correct = correct, type = type )
  
  return(x)
  
}

pred1 <- prediction.test(weights,results= D50, exact = TRUE) 

print.prediction.test <- function(x){
  cat("\n\t\tResults for the Prediction Test")

  cat("\n\nOf the", paste(x$m), " endpoints of interest,", paste(x$correct), "were correctly predicted.")
  if (x$type == 1){
    cat("\nCalculated using the normal approximation:")
  } else if (x$type == 2){
    cat("\nCalculated using the exact distribution:")
  }
  cat("\n----------------------------------------------------------------")
  cat("\nTm:", paste(format(x$statistic,digits=4)), 
      "\t Null hypothesized", paste("\U03D5:"), paste(format(x$null.value)),
      "\t p.value:", paste(format(x$p.value, digits = 4)))
}

print.prediction.test( pred1)



# Permutation prediction test ------------------------------------
# Inputs
# Weights ( a vector)
# Prediction results
# Null phi (if one value assume to be all the same), default = 0.50
# alpha, default = 0.05

nullphi <- c(.6,.5,.9,.5,.7)

prediction.test(weights, D50, exact = TRUE)
prediction.bootstrap(weights,D50, nullphi = c(0.1,0.1,0.1,0.1,0.1))


prediction.bootstrap <- function(weights, results, nullphi = 0.50, alpha = 0.05, sims = 5000){
  
  teststat <- weights%*%results
  ntests <- length(weights)
  correct <- sum(results)
  
  if (length(nullphi) == 1 | length(nullphi) == ntests)
  {
    boots <- matrix(NA,sims,1)
    for (g in 1:sims)
    {
      boots[g,] <- ifelse((rbinom(ntests,1,nullphi))%*%weights >= teststat,1,0)
    }
    pval <- mean(boots)

  } else{
    stop("nullphi needs to be either a single value or specified for every endpoint")
  }

  x <- list(statistic = teststat, p.value = pval, null.value = nullphi, m = ntests, alpha = alpha, correct = correct )

return(x)
  
  
}

p1 <- prediction.bootstrap(weights, D50, nullphi <- c(0.5))


print.prediction.bootstrap <- function(x){
  cat("\n\t\tResults for the Prediction Test")
  cat("\n\nOf the", paste(x$m), " endpoints of interest,", paste(x$correct), "were correctly predicted.")
  
  cat("\n----------------------------------------------------------------")
  if (length(x$null.value) == 1){
    cat("\nTm:", paste(format(x$statistic,digits=4)),
        "\t Null hypothesized", paste("\U03D5:"), paste(format(x$null.value)),
        "\t p.value:", paste(format(x$p.value, digits = 4)))
  }else if(length(x$null.value == ntests)){
    cat("\nTm:", paste(format(x$statistic,digits=4)),
        "\t Null hypothesized", paste("\U03D5:"), "varied",
        "\t p.value:", paste(format(x$p.value, digits = 4)))
    
  }
  
}

print.prediction.bootstrap( prediction.bootstrap(weights, D50, nullphi <- c(0.5,.8)) )

print.prediction.bootstrap(p1)

# Weights function function ------------------------------------
# Inputs
# data frame
# list of names
# group is variable name of either group or time variable
# type = either a group based test or pre-post
# Only works for two groups

library(tidyr)
library(dplyr)





prediction.weights(data = data, variables = variables, id = "id", type = "prepost", gtvar = "time", corr = "pearson")


prediction.weights <- function(data, variables,id, timevar, type = "group", corr = "pearson"){

  if(type == "group")
  {
    dataset <- data
    endpoints <- subset(dataset, select = c(variables))
    samplec <- cor(endpoints, method = corr)
    weights <- 1/rowSums(samplec^2)

  }else if( type == "prepost"){
    
    dataset <- data
    gtvar <- dataset[, gtvar]
    
    if (is.null(gtvar) )
    {
      stop("Provide a numeric gtvar to indicate pre and post observations.")
      
    }
    if( length(unique(gtvar)) > 2 ){
      
      stop("Expecting two values for gtvar.")
    }
    mutatefunc <- function (column) (
      dataset %>%
        arrange(id) %>%
        group_by(id) %>%
        mutate( !!paste0('diff.',as.name(column)) :=  !!as.name(column)-first(!!as.name(column))  ) -> dataset
      
    )
    
    for ( i in 1:length(variables))
    {
    x <- variables[i]
    dataset <- mutatefunc(x)
    }
  
    samplec <- cor( dataset[,c(paste0( "diff.",variables)) ], method = corr)
    weights <- 1/rowSums(samplec^2)

  }
  outweight <- list(weights)
  return(outweight)
}



# Prediction results function ------------------------------------
# Inputs

# FOR TESTING
#Group example
id <- 1:20
group <- c(rep(1,10), rep(2,10))
v1 <- rnorm(20,0,1)
v2 <- rnorm(20,0,1)
v3 <- rnorm(20,0,1)
v4 <- rnorm(20,0,1)
v5 <- rnorm(20,0,1)
v6 <- rnorm(20,0,1)


data <- as.data.frame(cbind(id, group, v1, v2 ,v3, v4, v5,v6))

Group <- "group"
variables <- c("v1", "v2" ,"v3", "v4", "v5","v6")
colnames(data)
groupvar <- "group"




id <- c(1:10,1:10)
time <- (c(rep(1,10), rep(2,10)))
v1 <- rnorm(20,0,1)
v2 <- rnorm(20,0,1)
v3 <- rnorm(20,0,1)
v4 <- rnorm(20,0,1)
v5 <- rnorm(20,0,1)
v6 <- rnorm(20,0,1)

data1 <- as.data.frame(cbind(id, time, v1, v2 ,v3, v4, v5,v6))

timevar <- data$time
variables <- c("v1", "v2" ,"v3", "v4", "v5","v6")
colnames(data)


# Either all Up, down, different or supply a matrix/data frame with two columns 
# (the variable name, and the prediction)

# What group (or time) is up or down will depend on the order of the variable, change order to change reference
# default second factor level - first factor level (reasoning pre-post data) 

# 1 is up, 2 is down, 3 is diff
#for type (wilcoxon or normal) 1 is wilcoxon, 2 is normal (wilcoxon by default if not specified)
#Group variable named group, or ID variable named id
predictions <- matrix(c("v1", 1,
                        "v2", 2,
                        "v3", 3,
                        "v4", 3,
                        "v5", 3,
                        "v6", 2),byrow = TRUE, ncol = 2)



#testing:
bound = "wilcoxon"
type = "prepost"
direction <- "mixed"
gtvar <- "time"
bound = "normal"


dataz <- data
colnames(dataz) <- c("id", "grp", "e1", "e2", "e3", "e4", "e5", "e6")
varz <- c( "e1", "e2", "e3", "e4", "e5", "e6")
variables <- varz
gtvar <- "grp"

dataset <- dataz

prediction.results(data1, variables = varaibles,direction = "up", type = "time", 
                   gtvar = "time")

prediction.results <- function(dataset, direction, bound = "wilcoxon", variables, type = "group", 
                               gtvar,  phi_0 = 0.50, predictions, location = "mean"){
  
  if (type == "group"){
    levels <- unique(factor(dataset[,gtvar]))
    
    if (location == "mean")
    {
    dataset %>%
      group_by(group) %>%
        summarise_at(all_of(variables), mean, na.rm = TRUE) -> groupmeans
    } else if (location == "median")
    {
    dataset %>%
      group_by(group) %>%
      summarise_at(all_of(variables), median, na.rm = TRUE) -> groupmeans
    }
    groupmeans <- groupmeans[order(groupmeans$group),] # Prediction calculate as Grp1-Grp2
    
    results <- groupmeans[1,variables] - groupmeans[2,variables]
    differences <- results
    
    if (direction == "up"){
      for ( z in 1:length(variables)){
        results[c(variables[z])] <- ifelse(results[c(variables[z])] > 0,1,0 )
        predictions <- rep(phi_0, length(variables))
      }
    }else if (direction == "down"){
      for ( z in 1:length(variables)){
        results[c(variables[z])] <- ifelse(results[c(variables[z])] < 0,1,0 )
        predictions <- rep(phi_0, length(variables))
      }
    }else if (direction == "mixed"){
      
      for ( j in 1:dim(predictions)[1])
      {
      
      rules <- predictions[predictions[,1] ==names(results[j]),]
      if (as.numeric(rules[2]) == 1){
        
        results[j] <- ifelse(results[j] > 0, 1, 0)
        
      }else if (as.numeric(rules[2]) == 2){
        
        results[j] <- ifelse(results[j] < 0, 1, 0)
      
      }else if(as.numeric(rules[2]) == 3){
        
        if (bound == "wilcoxon")
        {
          split <- dataset[,c(gtvar,names(results[j]))]
          split1 <- split[split[,gtvar]==as.numeric(levels[1]),]
          split2 <- split[split[,gtvar]==as.numeric(levels[2]),]
          
          results[j] <- ifelse(wilcox.test(split1[,names(results[j])],split2[,names(results[j])])$p.value <
                                 phi_0, 1,0)
          
          
        } else if (bound == "normal")
        {
          split <- dataset[,c(gtvar,names(results[j]))]
          split1 <- split[split[,gtvar]==as.numeric(levels[1]),]
          split2 <- split[split[,gtvar]==as.numeric(levels[2]),]
          
          obdiff <- mean(split1[,names(results[j])]) - mean(split2[,names(results[j])])
          results[j] <- ifelse(obdiff >  qnorm(1-(phi_0/2)) | obdiff < qnorm(phi_0/2), 1, 0)   
        }
      }
      }
    } 
  }else if (type == "prepost"){
    reference <- factor(dataset$time)[1]
    
    post <- dataset[dataset$time != reference,]
    pre <- dataset[dataset$time == reference,]
    if (location == "mean")
    {
    results <- colMeans(as.data.frame(post[,variables] - pre[,variables]))
    } else if (location == "median")
    {
      results <- colMedians(as.matrix(post[,variables] - pre[,variables]))
    }
    differences <- results
    
    if (direction == "up"){
      for ( z in 1:length(variables)){
        results[c(variables[z])] <- ifelse(results[c(variables[z])] > 0,1,0 )
      }
    }else if (direction == "down"){
      for ( z in 1:length(variables)){
        results[c(variables[z])] <- ifelse(results[c(variables[z])] < 0,1,0 )
      }
    }else if (direction == "mixed"){
    
    for ( j in 1:dim(predictions)[1])
    {
      
      rules <- predictions[predictions[,1] ==names(results[j]),]
      if (as.numeric(rules[2]) == 1){
        results[j] <- ifelse(results[j] > 0, 1, 0)
        
      }else if (as.numeric(rules[2]) == 2){
        
        results[j] <- ifelse(results[j] < 0, 1, 0)
        
      }else if(as.numeric(rules[2]) == 3){
        
        if (bound == "wilcoxon")
        {
          split <- dataset[,c(gtvar,names(results[j]))]
          split1 <- split[split[,gtvar]==as.numeric(reference),]
          split2 <- split[split[,gtvar]!=as.numeric(reference),]
          
          results[j] <- ifelse(wilcox.test(split1[,names(results[j])],split2[,names(results[j])])$p.value < phi_0, 1,0)
          
          
        } else if (bound == "normal")
        {
          split <- dataset[,c(gtvar,names(results[j]))]
          split1 <- split[split[,gtvar]==as.numeric(reference),]
          split2 <- split[split[,gtvar]!=as.numeric(reference),]
          
          obdiff <- mean(split2[,names(results[j])]) - mean(split1[,names(results[j])])
          results[j] <- ifelse(obdiff >  qnorm(1-(phi_0/2)) | obdiff < qnorm(phi_0/2), 1, 0)   
        }
      }
    }
  }
}
  
  outresults <- list(results, differences)
  return(outresults)

}



# Prediction plotting function ------------------------------------

#INPUTS
# Data set forplot, with the endpoint, median difference and result of prediction
variables <- c("v1", "v2", "v3", "v4", "v5")

results <- prediction.results(data, variables = variables,direction = "up", type = "group", 
                   gtvar = "group")

unlist(results[1])


prediction.plot <- function(predictionresults = results){
  
  end <- colnames(results[[2]])
  diff <- (unlist(results[2]))
  resulted <- (unlist(results[1]))

  forplot <- as.data.frame(   cbind(end, diff  , resulted  )   )
  forplot$diff <- as.numeric(as.character(forplot$diff))
  forplot$resulted <- as.numeric(as.character(forplot$resulted))
  
  
  forplot1 <- ggplot(data = forplot, aes(end, diff, fill = factor(resulted)   )  ) +
    geom_bar(aes(x = end, y = diff),stat='identity') +
    scale_y_continuous(limits=c(-0.5,0.75)) +
    geom_bar(forplot, mapping = aes(end) ,alpha=0, size=1, color="black", stat='identity')+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(title = "Differences across endpoints", x = "", y = "Standardized median difference", fill = "Prediction\nresults") +
    scale_fill_manual(labels = c("Incorrect", "Correct"), values = c("White", "Black")) + 
    theme(plot.title = element_text(hjust = 0.5)) 
  
  
  +
    geom_segment(aes(x = 8.5, y = 0.46, xend = 9.5, yend = 0.46),size=1, linetype = 3)+
    geom_segment(aes(x = 8.5, y = -0.46, xend = 9.5, yend = -0.46),size=1,  linetype = 3)+
    geom_point(aes(x=9, y=0), shape = 1, fill = "white",size = 6, stroke = 1)
  
  class(forplot$end)
  
  output(forplot1)
  
}










