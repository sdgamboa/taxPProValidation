
args <- commandArgs(trailingOnly = TRUE)
physName <- args[[1]]

library(taxPPro)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(tidyr)

## Define funcions ####
extractPredValues <- function(testSet, predSet) {
    testValues <- testSet$trait1
    names(testValues) <- testSet[['assembly']]
    predValues <- predSet$trait1
    names(predValues) <- predSet[[1]]
    predValues <- predValues[names(predValues) %in% names(testValues)]
    predValues <- predValues[names(testValues)]
    predValues
}

calcMetrics <- function(testSet, predSet) {
    actual_values <- testSet$trait1
    predicted_values <- extractPredValues(testSet, predSet)
    
    mse <- mean((predicted_values - actual_values)^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(predicted_values - actual_values))
    
    ## use lm
    ss_total <- sum( (actual_values - mean(actual_values))^2 )
    ss_residual <- sum( (actual_values - predicted_values)^2 )
    r_squared <- 1 - (ss_residual / ss_total)
    
    mape <- mean(abs((actual_values - predicted_values) / actual_values)) * 100
    
    evaluation_results <- data.frame(
        # Metric = c("MSE", "RMSE", "MAE", "R_squared", "MAPE"),
        Metric = c("MSE", "RMSE", "R_squared"),
        # Value = c(mse, rmse, mae, r_squared, mape)
        Value = c(mse, rmse, r_squared)
    )
    
    evaluation_results
}

plot <- function(testSet, predSet) {
    testValues <- testSet$trait1
    names(testValues) <- testSet$assembly
    predValues <- extractPredValues(testSet, predSet)
    data.frame(actual = testValues, pred = predValues) |> 
        ggplot(aes(actual, pred)) +
        geom_point()
}

## Read files ####
testSetsF <- list.files('test_folds', full.names = TRUE)
testSetsF <- grep(paste0('\\b', physName, '_test'), testSetsF, value = TRUE)
testSets <- vector('list', length(testSetsF))
for (i in seq_along(testSets)) {
    ## the extension of the files is csv, but they're actually tsv
    testSets[[i]] <- readr::read_tsv(testSetsF[[i]], show_col_types = FALSE)
}

predSetsF <- list.files('predicted', full.names = TRUE)
predSetsF <- grep(paste0('\\b', physName, '_train'), predSetsF, value = TRUE)
predSets <- vector('list', length(predSetsF))
for (i in seq_along(predSets)) {
    predSets[[i]] <- readr::read_tsv(predSetsF[[i]], show_col_types = FALSE)
    
}

## Plot ####
plotList <- map2(testSets, predSets, plot)
p <- ggarrange(plotlist = plotList, ncol = 5, nrow = 2)
ggsave(filename = paste0(physName, '.png'), plot = p, width = 12, height = 6, units = 'in')

## Metrics ####
summaryDat <- map2(testSets, predSets, calcMetrics) |> 
    bind_rows() |> 
    group_by(Metric) |> 
    summarise(
        mean = round(mean(Value, na.rm = TRUE), 3),
        sd = round(sd(Value, na.rm = TRUE), 3)
    ) |> 
    ungroup() 


stats <- bind_rows(testSets) |> 
    summarise(
        mean_real_values = mean(trait1, na.rm = TRUE),
        min_real_values = min(trait1, na.rm = TRUE),
        max_real_values = max(trait1, na.rm = TRUE),
        sd_real_values = sd(trait1, na.rm = TRUE)
    ) |> 
    as.data.frame()

d1 <- as.data.frame(matrix(data = summaryDat$mean, nrow = 1))
colnames(d1) <- paste0(summaryDat$Metric, '_mean')

d2 <- as.data.frame(matrix(data = summaryDat$sd, nrow = 1))
colnames(d2) <- paste0(summaryDat$Metric, '_sd')
   
outputDat <- do.call('cbind', list(d1, d2, stats)) |> 
    mutate(phys = physName) |> 
    relocate(
        phys, R_squared_mean, R_squared_sd, MSE_mean, MSE_sd,
        RMSE_mean, RMSE_sd
    )

#histList <- map(predSets, ~ {
#    .x |> 
#        ggplot(aes(metadata_NSTI)) +
#        geom_histogram()
#})

write.table(outputDat, file = paste0(physName, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
