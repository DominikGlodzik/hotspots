library(MASS)
library(boot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(mosaic)
library(pscl)

# Axis rotate function
rotatedAxisElementText = function(angle,position='x'){
    angle     = angle[1]; 
    position  = position[1]
    positions = list(x=0,y=90,top=180,right=270)
    if(!position %in% names(positions))
        stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
    if(!is.numeric(angle))
        stop("'angle' must be numeric",call.=FALSE)
    rads  = (angle - positions[[ position ]])*pi/180
    hjust = 0.5*(1 - sin(rads))
    vjust = 0.5*(1 + cos(rads))
    element_text(angle=angle,vjust=vjust,hjust=hjust)
}


evaluateBgModel <- function(bgData, formula.nb, iterBoot = 100, propTrain = 0.7){
    # Performance metrics for SV null model
    # Estimates mean squared error (MSE) of null model predictions in bootstrapped datasets compared with random shuffling of predictions. 
    # 95 % confidence interval of bootstrapped MSE estimated by the quantile method.
    # p-value represents the proportion of permutation MSE values equal to or lower than the median MSE from prediction.
    # The p-value cannot be lower than 1 divided by the number of bootstrap iterations
    # Arguments:
    #   iterBoot: number of bootstrap iterations, default = 100
    #   propTrain: proportion of dataset used for model training, default 0.7  (test set = 1-propTrain)
    #   formula.nb: formula as input to negative binomial regression null model
    #   bgData: data frame of all data to use for model validation
    
    # Function
    val_function <- function(data, inds){
        temp_data <- data[inds,]
        nr <- nrow(temp_data)
        in_train <- sample(nr, propTrain*nr)
        train <- temp_data[in_train,]
        test <- temp_data[-in_train,]
        temp.fit <-  glm.nb(formula.nb, data = train)
        test$sv.count.pred <- exp(predict(temp.fit, test, type = "link", se.fit=TRUE)$fit)
        # MSE position 1 is observed vs predicted; position 2 is observed vs permutation of predicted
        MSE <- c(costFun(test$sv.count, test$sv.count.pred), costFun(test$sv.count, shuffle(test$sv.count.pred)))
        return(MSE)
    }
    
    # Defining MSE function
    costFun <- function(y, yhat) mean((y - yhat)^2) 
    
    # Run bootstrap
    cv.boot <- boot(data = bgData, statistic = val_function, R = iterBoot) # Bootstrapping test/train
    
    # Bootstrapping errors
    boot_errors <- data.frame(model = c('Null-model', 'Permutation'),
                              median = c(median(cv.boot$t[,1]), median(cv.boot$t[,2])),
                              CI2.5 = c(quantile(cv.boot$t[,1], probs = 0.025), quantile(cv.boot$t[,2], probs = 0.025)),
                              CI97.5 = c(quantile(cv.boot$t[,1], probs = 0.975), quantile(cv.boot$t[,2], probs = 0.975)))
    
    boot_errors_long <- boot_errors %>%
        melt(id.vars = c('model'),
             value.name = 'MSE')
    
    raw.mse <- rbind(data.frame(MSE = cv.boot$t[,1],
                                model = 'Null-model'),
                     data.frame(MSE = cv.boot$t[,2],
                                model = 'Permutation'))
    
    # p-value estimation
    median_prediction <- filter(boot_errors_long, model == 'Null-model', variable == 'median')$MSE
    permutation_all <- filter(raw.mse, model == 'Permutation')$MSE
    empirical_p <- sum(permutation_all <= median_prediction)/nrow(permutation_all) 
    
    boot_errors_long$p.val <- ifelse(length(empirical_p) == 0, rep(0, nrow(boot_errors_long)), empirical_p)
    
    # Generate test set for summary statistics
    
    # Final utput list
    set.seed(1)
    inTrainModel <- sample(nrow(bgData), nrow(bgData)*0.70)
    trainData <- bgData[inTrainModel,]
    testData <- bgData[-inTrainModel,]
    nb.fit.train <- glm.nb(formula.nb , data = trainData) # negative binomial fit
    testData$sv.count.predicted <- exp(predict(nb.fit.train, testData, type = "link", se.fit=TRUE)$fit)
    
    validationOut <- list()
    validationOut$MSEsummary <- boot_errors_long
    validationOut$rawBootMSE <- raw.mse
    validationOut$testData <- testData
    
    return(validationOut)
}


plotCompareModels <- function(inputList){
    # Generate summary plot for comparison of SV null model predictions versus permutation
    # Takes as input a list of 1) summary statistics for both models and 2) raw MSE values
    # Input argument corresponds to the output from bgModelEval
    
    summaryMSE <- inputList[[1]]
    rawMSE <- inputList[[2]]
    nBoot <- nrow(rawMSE)/2
    
    scale_col_model <- function(...){
        ggplot2:::manual_scale('colour', 
                               values = setNames(c('darkgray', "#E41A1C"),
                                                 c('Permutation', 'Null-model')), 
                               ...)
    }
    
    p <- ggplot() +
            geom_jitter(data = rawMSE, aes(model, MSE, col = model), size = 1, alpha = 0.2)+
            geom_point(data = filter(summaryMSE, variable == 'median'), aes(model, MSE, col = model), size = 3)+
            geom_path(data = filter(summaryMSE, variable != 'median'), aes(model, MSE, group = model, col = model), size = 1)+
            scale_col_model() +
            guides(col = F)+
            labs(y = 'Mean Squared Error, median (95 % CI) \n bootstraped null model vs permutation',
                 col = 'Null model')+
            geom_text(data = summaryMSE,
                      aes(x = 1.5, 
                          y = min(MSE)/2, 
                          label = ifelse(p.val < 1/nBoot, paste0('p < ', as.character(1/nBoot)), paste0('p = ', round(p.val,3)))), size = 5)+
            scale_y_continuous(limits = c(0,max(rawMSE$MSE)))+
            theme_bw()+
            theme(text = element_text(size = 12),
                  axis.text.x = rotatedAxisElementText(90, 'top'),
                  axis.title.x = element_blank())
    return(p)
}
