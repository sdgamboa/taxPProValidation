library(taxPPro)
library(bugphyzz)
library(ape)
library(dplyr)
library(cvTools)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data
gn_tips <- ltp$gn_tips

physName <- 'growth temperature'

gt <- physiologies(physName)[[1]]

modifyNumeric <- function(x) {
    x |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown') |>
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max)) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_min))) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_max))) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
           .data$Confidence_in_curation, n = 1, with_ties = FALSE
        ) |>
        dplyr::mutate(
            Attribute_value = mean(.data$Attribute_value_min, .data$Attribute_value_max)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-.data$Attribute_value_min, -.data$Attribute_value_max) |>
        dplyr::distinct()
}

dat <- modifyNumeric(gt)
dat <- dat[dat$NCBI_ID %in% tip_data$taxid,]
valuesDat <- left_join(tip_data, dat, by = c('taxid' = 'NCBI_ID')) |>
    select(tip_label, Attribute_value) |>
    filter(!is.na(Attribute_value))
colnames(valuesDat) <- c('assembly', 'trait1')

## Create folds
set.seed(1234)
foldN <- cvFolds(n = nrow(valuesDat), K = 10)
testSets <- vector('list', 10)
trainSets <- vector('list', 10)
for (i in 1:10) {
    foldName <- paste0('Fold', i)
    testSets[[i]] <- valuesDat[foldN$subsets[foldN$which == i],]
    names(testSets)[i] <- foldName
    trainSets[[i]] <- valuesDat[foldN$subsets[foldN$which != i],]
    names(trainSets)[i] <- foldName
}

## Write test and train sets
for (i in 1:10) {
    testFileName <- paste0(gsub(' ', '_', physName), '_test_Fold_', i, '.csv')
    dir.create('test_folds')
    testFileName <- file.path('test_folds', testFileName)
    write.table(x = testSets[[i]], file = testFileName, quote = FALSE, row.names = FALSE, sep = '\t')
    
    trainFileName <- paste0(gsub(' ', '_', physName), '_train_Fold_', i, '.csv')
    dir.create('train_folds')
    trainFileName <- file.path('train_folds', trainFileName)
    write.table(x = trainSets[[i]], file = trainFileName, quote = FALSE, row.names = FALSE, sep = '\t')
    
}
