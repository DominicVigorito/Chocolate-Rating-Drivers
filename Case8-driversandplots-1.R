library(regclass)
library(multcompView)
library(stringr)
library(DALEX)
library(caret)
library(gbm)
library(randomForest)
library(pROC)


########################################################
############## Replace Rare Levels Code ################
########################################################

replace_rare_levels <- function(x, threshold = 20, newname = "Other")  {
  x <- factor(x)
  rare.levels <- names(which(sort(table(x)) <= threshold))
  if (length(rare.levels) == 0) { return(x) }
  levels(x)[which(levels(x) %in% rare.levels)] <- newname
  ST <- sort(table(x))
  if (ST[newname] <= threshold) {
    levels.to.combine <- which(levels(x) %in% c(newname, 
                                                names(ST)[2]))
    levels(x)[levels.to.combine] <- newname
    rare.levels <- c(rare.levels, names(ST)[2])
  }
  return(x)
}

########################################################
############## Examine Driver Y categor ################
########################################################

examine_driver_Ycat <- function(formula,data,sort=TRUE,inside=TRUE,equal=TRUE) { 
  require(regclass)
  require(multcompView)
  require(stringr)
  
  FORM <- as.character(formula)
  temp <- str_extract(FORM,"^.*==")
  temp <- temp[!is.na(temp)]
  y.label <- gsub(" ==","",temp)
  temp <- str_extract(FORM,'\\"[:alnum:]*')
  temp <- temp[!is.na(temp)]
  level <- substr(temp,2,nchar(temp))
  FORM <- as.formula(formula)
  
  variables <- as.character(attr(terms(FORM), "variables"))[-1]
  x.label <- variables[2]
  data <- data[,c(x.label,y.label)]
  if (head(class(data[,y.label]),1) != "ordered") {
    data[,y.label] <- factor(data[,y.label])
  }
  data[,y.label] <- factor(data[,y.label],ordered=TRUE,levels=c(level,setdiff(levels(data[,y.label]),level)))
  if(nlevels(data[,y.label])>2) { levels(data[,y.label]) <- c(levels(data[,y.label])[1],rep(paste("Not",level),nlevels(data[,y.label])-1)) }
  x <- data[,x.label]
  y <- data[,y.label]
  
  color <- FALSE
  labelat=c(); xlab=c(); ylab=c(); magnification=1
  complete.x <- which(!is.na(x))
  complete.y <- which(!is.na(y))
  complete.cases <- intersect(complete.x, complete.y)
  x <- x[complete.cases]
  y <- y[complete.cases]
  if (head(class(x),1) != "ordered") {
    x <- factor(x)
  }
  if (head(class(y),1) != "ordered") {
    y <- factor(y)
  }
  data[,1] <- x
  data[,2] <- y
  n <- length(x)
  nx.levels <- length(unique(x))
  ny.levels <- length(unique(y))
  if (nx.levels < 2 | ny.levels < 2) {
    stop(paste("Error:  need at least 2 levels to proceed.  x has", 
               nx.levels, "and y has", ny.levels))
  }
  if (nx.levels > 100) {
    stop(paste("Error:  function not designed for more than 100 levels of x"))
  }
  
  xlevel.names <- levels(x)
  if(length(labelat)>0) { xlevel.names[!(xlevel.names %in% labelat)] <- "" }
  
  ylevel.names <- levels(y)
  CONT.TAB <- table(x, y)
  CONT.TAB <- addmargins(CONT.TAB)
  rownames(CONT.TAB)[nx.levels + 1] <- "Total"
  colnames(CONT.TAB)[ny.levels + 1] <- "Total"
  O <- matrix(table(x, y), nrow = nx.levels, ncol = ny.levels)
  E <- (apply(O, 1, sum) %o% apply(O, 2, sum))/n
  plot(0, 0, col = "white", xlim = c(0, 1.3), ylim = c(-0.05, 
                                                       1), xlab=ifelse(equal==TRUE,"",x.label),cex.main=0.7,main = "", ylab = paste("Probability",y.label,"is",level), axes = FALSE)
  if(equal==FALSE) { axis(1, at = seq(0, 1, 0.05)) }
  axis(2, at = seq(0, 1, 0.05))
  COLORS <- grey(seq(0.1, 0.7, length = ny.levels))
  marginal.y <- apply(O, 2, sum)/n
  break.y <- c(0, cumsum(marginal.y))
  for (i in 1:ny.levels) {
    rect(1.01, break.y[i], 1.11, break.y[i + 1], col = COLORS[i])
    text(1.1, (break.y[i + 1] + break.y[i])/2, ylevel.names[i], 
         srt = 0, pos = 4, cex = 1)
  }
  marginal.x <- apply(O, 1, sum)/n
  
  
  if(equal==FALSE) { break.x <- c(0, cumsum(marginal.x)) } else { break.x <- seq(0,1,length=1+length(xlevel.names)) }
  
  for (i in 1:nx.levels) {
    marginal.y <- O[i, ]/sum(O[i, ])
    break.y <- c(0, cumsum(marginal.y))
    for (j in 1:ny.levels) {
      rect(break.x[i], break.y[j], break.x[i + 1], break.y[j + 
                                                             1], col = COLORS[j])
    }
    if(inside==TRUE) { 
      text((break.x[i + 1] + break.x[i])/2, 0.5, xlevel.names[i],cex=magnification,srt=90,col="white")
    } else {
      text((break.x[i + 1] + break.x[i])/2, -0.05, xlevel.names[i],cex=magnification) }
    
  }
  lines(c(0,1),rep(mean(data[,y.label]==level),2),lwd=2,col="red")
  if(equal==TRUE) { text(0.5,-0.05,x.label)   }
  FORM3 <- formula(paste(y.label,"=='",level,"'~", x.label,sep=""))
  SUMMARY <- aggregate(FORM3,data=data,FUN=mean)
  AOV <- aov(FORM3,data=data)
  TUKEY <- TukeyHSD(AOV)
  LETTERS <- multcompLetters4(AOV,TUKEY) 
  SUMMARY$letters <- LETTERS[[1]][1]$Letters[match(SUMMARY[,1],names( LETTERS[[1]][1]$Letters ) )]
  SUMMARY$n <- as.numeric(table(data[,x.label]))
  names(SUMMARY)[2] <- paste("Prob",level,sep="")
  if(sort==TRUE) { SUMMARY <- SUMMARY[order(SUMMARY[,2],decreasing=TRUE),] }
  rownames(SUMMARY) <- NULL
  print(SUMMARY)
  R2 <- summary(lm( formula, data ) )[8]$r.squared
  cat(paste("\nDriver Score:",round(R2,digits=3),"\n"))
  cat(paste("Scores range between 0 and 1.  Larger scores = stronger driver.\n"))
  cat(paste("Although context dependent, values above 0.02 or so are 'reasonably strong' drivers.\n"))
}

########################################################
############## Examine Driver Y numeric ################
########################################################

examine_driver_Ynumeric <- function(formula,data,sort=TRUE,lambda=NA) { 
  require(regclass)
  require(multcompView)
  require(stringr)
  
  FORM <- as.formula(formula)
  variables <- as.character(attr(terms(FORM), "variables"))[-1]
  x <- data[,variables[2]]
  y <- data[,variables[1]]
  complete.x <- which(!is.na(x))
  complete.y <- which(!is.na(y))
  complete.cases <- intersect(complete.x, complete.y)
  x <- x[complete.cases]
  y <- y[complete.cases]
  plot(y~x,xlab=variables[2],ylab=variables[1])
  M <- lm(FORM,data)
  if(class(x)[1] %in% c("integer","numeric")) { 
    #if(is.na(lambda)) { SS <- smooth.spline(x,y) } else { SS <- smooth.spline(x,y,lambda) }
    #lines(SS,col="red",lwd=3)
    abline(M,col="blue",lwd=3)
  }
  abline(h=mean(y),col="red")
  if(class(x)[1]  %in% c("integer","numeric") ) {
    print( summary(M) )
  } else { 
    
    SUMMARY <- aggregate(FORM,data=data,FUN=mean)
    AOV <- aov(FORM,data=data)
    TUKEY <- TukeyHSD(AOV)
    LETTERS <- multcompLetters4(AOV,TUKEY) 
    SUMMARY$letters <- LETTERS[[1]][1]$Letters[match(SUMMARY[,1],names( LETTERS[[1]][1]$Letters ) )]
    SUMMARY$n <- as.numeric(table(x))
    names(SUMMARY)[2] <- paste("Avg",variables[1],sep="")
    if(sort==TRUE) { SUMMARY <- SUMMARY[order(SUMMARY[,2],decreasing=TRUE),] }
    rownames(SUMMARY) <- NULL
    print(SUMMARY) }
  
  cat(paste("\nDriver Score:",round(summary(M)[8]$r.squared,digits=3),"\n"))
  cat(paste("Caution: misleading if driver is numerical and trend isn't linear.\n"))
  cat(paste("Scores range between 0 and 1.  Larger scores = stronger driver.\n"))
  cat(paste("Although context dependent, values above 0.02 or so are 'reasonably strong' drivers.\n"))
}


##################################################################################
#Option 14 (Chocolate Quality)
##################################################################################

CHOCOLATE <- read.csv("Case5-chocolate-1.csv",stringsAsFactors = TRUE)

#Remove identifier variables and other things that aren't interesting
CHOCOLATE$X <- NULL
CHOCOLATE$ref <- NULL
CHOCOLATE$review_date <- NULL
CHOCOLATE$beans<-NULL #only 1 value

#Combine rare levels
top <- 8 #going to write code to consider only top 10, rest get combined to "Other"


x <- combine_rare_levels(CHOCOLATE$company,threshold = sort(table(CHOCOLATE$company),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
CHOCOLATE$company <- x

x <- combine_rare_levels(CHOCOLATE$company_location,threshold = sort(table(CHOCOLATE$company_location),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
CHOCOLATE$company_location <- x

x <- combine_rare_levels(CHOCOLATE$country_of_bean_origin,threshold = sort(table(CHOCOLATE$country_of_bean_origin),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
CHOCOLATE$country_of_bean_origin <- x

x <- combine_rare_levels(CHOCOLATE$specific_bean_origin_or_bar_name,threshold = sort(table(CHOCOLATE$specific_bean_origin_or_bar_name),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
CHOCOLATE$specific_bean_origin_or_bar_name <- x

x <- combine_rare_levels(CHOCOLATE$first_taste,threshold = sort(table(CHOCOLATE$first_taste),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
CHOCOLATE$first_taste <- x

x <- combine_rare_levels(CHOCOLATE$second_taste,threshold = sort(table(CHOCOLATE$second_taste),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
levels(x)[1] <- "None"
CHOCOLATE$second_taste <- x

x <- combine_rare_levels(CHOCOLATE$third_taste,threshold = sort(table(CHOCOLATE$third_taste),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
levels(x)[1] <- "None"
CHOCOLATE$third_taste <- x

x <- combine_rare_levels(CHOCOLATE$fourth_taste,threshold = sort(table(CHOCOLATE$fourth_taste),dec=TRUE)[top])$values
levels(x)[which(levels(x) %in% c("Combined") )] <- "Other"
summary(x)
levels(x)[1] <- "None"
CHOCOLATE$fourth_taste <- x


#Combos
TREE <- rpart(rating ~ . , data=CHOCOLATE,cp=0,minbucket=25)
summarize_tree(TREE)
TREE$cptable
TREE <- rpart(rating ~ . , data=CHOCOLATE,cp=0.006,minbucket=25)
visualize_model(TREE)



#Make special datasets for examining partial dependence plots, interactions, and breakdown plots 
TRAIN.PREDICTORS <- CHOCOLATE
TRAIN.PREDICTORS$rating <- NULL
TRAIN.TARGET <- CHOCOLATE$rating

#Estimating generalization error (build on a 70% training set and evaluate on a 30% holdout sample
#Normally this is done with K-fold crossvalidation and caret, but this is faster!
#Drawback:  only one guess of the generalization error! 
DATA <- CHOCOLATE
y <- "rating"
set.seed(479); train.rows <- sample(1:nrow(DATA),0.7*nrow(DATA))
MODEL <- randomForest(formula(paste(y,"~.")),data=DATA[train.rows,])
#Estimated RMSE on holdout sample (one measure of generalization error)
postResample(predict(MODEL,newdata=DATA[-train.rows,]),DATA[-train.rows,y] )
#RMSE on holdout if predicted everything to be average value (naive model)
sd( DATA[-train.rows,y] - mean(DATA[train.rows,y]) )


#Building an explainer with a random forest
library(randomForest)
FOREST <- randomForest(rating~.,data=CHOCOLATE)
summarize_tree(FOREST)
model_explainer <- explain(model = FOREST, data = TRAIN.PREDICTORS,  y = TRAIN.TARGET)
model_vi <- model_parts(model_explainer, loss_function = loss_root_mean_square)
plot(model_vi)
model_vi


#Partial dependence plots 
numerics <- which( unlist(lapply(CHOCOLATE,function(x)tail(class(x),1))) %in% c("numeric","integer") ); names(CHOCOLATE)[numerics]
cats <- which( unlist(lapply(CHOCOLATE,function(x)tail(class(x),1))) %in% c("factor","character") ); names(CHOCOLATE)[cats]

#cat( paste('plot( model_profile(model_explainer,type = "accumulated"), variables = "',names(CHOCOLATE)[numerics],'")\n',sep="" ))
plot( model_profile(model_explainer,type = "accumulated"), variables = "cocoa_percent")
plot( model_profile(model_explainer,type = "accumulated"), variables = "counts_of_ingredients")

#cat(paste('plot( model_profile(model_explainer,type = "partial", variables = "',names(CHOCOLATE)[cats],'"))\n',sep=""))
plot( model_profile(model_explainer,type = "partial", variables = "company"))
plot( model_profile(model_explainer,type = "partial", variables = "company_location"))
plot( model_profile(model_explainer,type = "partial", variables = "country_of_bean_origin"))
plot( model_profile(model_explainer,type = "partial", variables = "specific_bean_origin_or_bar_name"))
plot( model_profile(model_explainer,type = "partial", variables = "cocoa_butter"))
plot( model_profile(model_explainer,type = "partial", variables = "vanilla"))
plot( model_profile(model_explainer,type = "partial", variables = "lecithin"))
plot( model_profile(model_explainer,type = "partial", variables = "salt"))
plot( model_profile(model_explainer,type = "partial", variables = "sugar"))
plot( model_profile(model_explainer,type = "partial", variables = "sweetener_without_sugar"))
plot( model_profile(model_explainer,type = "partial", variables = "first_taste"))
plot( model_profile(model_explainer,type = "partial", variables = "second_taste"))
plot( model_profile(model_explainer,type = "partial", variables = "third_taste"))
plot( model_profile(model_explainer,type = "partial", variables = "fourth_taste"))



#Examine interaction between num/cat or cat/cat variables (see variables in numerics and cats created above)
set.seed(2022); profile_group <- model_profile(explainer = model_explainer, 
                                               variables = c("cocoa_percent"), 
                                               groups = "first_taste", type = "partial")  
plot(profile_group)

#Breakdown plot showing how particular prediction is pieced together
#Choose two types of rows that tell a good story
which(CHOCOLATE$rating >= quantile(CHOCOLATE$rating,.99))  #Some of the rows with highest y
which(CHOCOLATE$rating >= quantile(CHOCOLATE$rating,.49) & CHOCOLATE$rating <= quantile(CHOCOLATE$rating,.51))  #Some of the rows with typical value of y
which(CHOCOLATE$rating <= quantile(CHOCOLATE$rating,.01))  #Some of the rows with highest y
specific.row <- CHOCOLATE[415,]  
plot( predict_parts(explainer = model_explainer, new_observation = specific.row, type = "break_down") )



#Single Drivers

#cat( paste('examine_driver_Ynumeric(rating ~',names(CHOCOLATE),', data=CHOCOLATE)\n' ))
examine_driver_Ynumeric(rating ~ company , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ company_location , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ country_of_bean_origin , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ specific_bean_origin_or_bar_name , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ cocoa_percent , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ counts_of_ingredients , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ cocoa_butter , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ vanilla , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ lecithin , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ salt , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ sugar , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ sweetener_without_sugar , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ first_taste , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ second_taste , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ third_taste , data=CHOCOLATE)
examine_driver_Ynumeric(rating ~ fourth_taste , data=CHOCOLATE)


