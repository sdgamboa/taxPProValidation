library(taxPPro)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)

## Define funcions ####
extractPredValues <- function(testSet, predSet) {
    testValues <- testSet$trait1
    names(testValues) <- testSet$taxa
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
    
    ss_total <- sum((actual_values - mean(actual_values))^2)
    ss_residual <- sum((actual_values - predicted_values)^2)
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
    names(testValues) <- testSet$taxa
    predValues <- extractPredValues(testSet, predSet)
    data.frame(actual = testValues, pred = predValues) |> 
        ggplot(aes(actual, pred)) +
        geom_point()
}

## Read files ####
testSetsF <- list.files('test_folds', full.names = TRUE)
testSetsF <- grep('growth_temperature', testSetsF, value = TRUE)
testSets <- vector('list', length(testSetsF))
for (i in seq_along(testSets)) {
    ## the extension of the files is csv, but they're actually tsv
    testSets[[i]] <- readr::read_tsv(testSetsF[[i]], show_col_types = FALSE)
}

predSetsF <- list.files('predicted', pattern = 'predict_traits.tab', recursive = TRUE, full.names = TRUE)
predSetsF <- grep('growth_temperature', predSetsF, value = TRUE)
predSets <- vector('list', length(predSetsF))
for (i in seq_along(predSets)) {
    predSets[[i]] <- readr::read_tsv(predSetsF[[i]], show_col_types = FALSE)
    
}

## Plot ####
plotList <- map2(testSets, predSets, plot)
p <- ggarrange(plotlist = plotList, ncol = 5, nrow = 2)
ggsave(filename = 'growth_temperature.png', plot = p, width = 15, height = 8, units = 'in')

## Metrics ####
summaryDat <- map2(testSets, predSets, calcMetrics) |> 
    bind_rows() |> 
    group_by(Metric) |> 
    summarise(
        mean = mean(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE)
    ) |> 
    ungroup()
write.table(summaryDat, file = 'growth_temperature.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
