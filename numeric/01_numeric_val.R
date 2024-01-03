
args <- commandArgs(trailingOnly = TRUE)
# args <- c('growth_temperature', 'all')

library(taxPPro)
library(dplyr)
library(bugphyzz)
library(phytools)
library(castor)
library(cvTools)
library(purrr)
library(tidyr)
library(ggplot2)
library(logr)

phys_name <- args[[1]]
rank_arg <- args[[2]]

if (rank_arg == 'all') {
    rank_var <- c('genus', 'species', 'strain')
} else {
    rank_var <- rank_arg
}

logfile <- paste0(phys_name, '_', rank_arg, '_castor-ltp')
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

msg <- paste0('Running validation for: ', phys_name)
log_print(msg)

ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data

# label_data <- bind_rows(
#     distinct(select(as_tibble(tip_data), label = tip_label, taxid)),
#     distinct(select(as_tibble(node_data), label = node_label, taxid))
# )

## bugphyzz data ####
phys <- physiologies(gsub('_', ' ', phys_name))[[1]]
# cg <- physiologies('coding genes')[[1]]
# w <- physiologies('width')[[1]]

modifyNumeric <- function(x) {
    x |>
        dplyr::filter(!is.na(.data$NCBI_ID) | !.data$NCBI_ID == 'unkown') |>
        dplyr::filter(!is.na(.data$Attribute_value_min)) |>
        dplyr::filter(!is.na(.data$Attribute_value_max)) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_min))) |>
        dplyr::filter(!is.infinite(abs(.data$Attribute_value_max))) |>
        dplyr::mutate(Rank = taxizedb::taxid2rank(NCBI_ID, db = 'ncbi')) |> 
        dplyr::mutate(Taxon_name = taxizedb::taxid2name(NCBI_ID, db = 'ncbi')) |> 
        dplyr::mutate(NCBI_ID = addRankPrefix(NCBI_ID, Rank)) |> 
        dplyr::filter(!is.na(Rank) & !is.na(Taxon_name) & !is.na(NCBI_ID)) |>
        dplyr::group_by(.data$NCBI_ID) |>
        dplyr::slice_max(
           .data$Confidence_in_curation, n = 1, with_ties = FALSE
        ) |>
        dplyr::mutate(
            Attribute_value = mean(.data$Attribute_value_min, .data$Attribute_value_max)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-Attribute_value_min, -Attribute_value_max) |>
        dplyr::distinct()
}

dat <- modifyNumeric(phys) |> 
    filter(Rank %in% rank_var)

if (!nrow(dat)) {
    log_print('not enough data for propagation')
    quit(save = 'no')
    log_close()
}
    
dat$NCBI_ID[which(duplicated(dat$NCBI_ID))]

fdat <- dat |> 
    # filter(Rank %in% rank_var) |> 
    filter(NCBI_ID %in% unique(tip_data$NCBI_ID))

if (!nrow(fdat)) {
    log_print('not enough data for propagation')
    quit(save = 'no')
    log_close()
}
    
counts <- data.frame(
    physiology = gsub('_', ' ', phys_name),
    rank = rank_arg,
    ltp_and_bugphyzz = length(unique(fdat$NCBI_ID)),
    bugphyzz = length(unique(dat$NCBI_ID)),
    ltp = length(unique(tip_data$NCBI_ID))
)

tip_data_annotated <- left_join(tip_data, fdat, by = c('NCBI_ID')) |>
    select(tip_label, Attribute_value) # I don't need to filter NAs here

allTips <- tip_data_annotated$Attribute_value
names(allTips) <- tip_data_annotated$tip_label

values <- allTips[!is.na(allTips)]

## Create train and test sets ####
set.seed(1234)
foldN <- cvFolds(n = length(values), K = 10)
testSets <- vector('list', 10)
trainSets <- vector('list', 10)
for (i in 1:10) {
    foldName <- paste0('Fold', i)
    testSets[[i]] <- values[foldN$subsets[foldN$which == i]]
    names(testSets)[i] <- foldName
    trainSets[[i]] <- values[foldN$subsets[foldN$which != i]]
    names(trainSets)[i] <- foldName
}

## hidden-state-prediction ####
hsp <- vector('list', 10)
for (i in seq_along(trainSets)) {
    names(hsp)[i] <- names(trainSets)[i]

    x <- testSets[[i]]
    x[!is.na(x)] <- NA
    y <- trainSets[[i]]
    z = c(x, y)
    a = allTips[!names(allTips) %in% names(z)]
    input_vector <- c(z, a)
    input_vector <- input_vector[tree$tip.label]

    res <- hsp_squared_change_parsimony(
        ## slightly better results than hsp_independent_contrasts
        ## for growth temp
        tree = tree, tip_states = input_vector, weighted = TRUE,
        check_input = TRUE
    )

    # res <- hsp_independent_contrasts( # pic in ape::ace ?
    #     tree = tree, tip_states = input_vector, weighted = TRUE,
    #     check_input = TRUE
    # )

    statesDF <- data.frame(
        label = tree$tip.label,
        value = res$states[1:Ntip(tree)]
    )

    states <- statesDF$value
    names(states) <- statesDF$label

    # commonNames <- intersect(names(states), names(testSets[[i]]))
    # hsp[[i]] <- states[commonNames]
    hsp[[i]] <- states[names(testSets[[i]])]
}

## Calculate metrics
metrics <- map2(hsp, testSets, ~ {
    predicted_values <- .x
    actual_values <- .y

    mse <- mean((predicted_values - actual_values)^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(predicted_values - actual_values))

    ss_total <- sum((actual_values - mean(actual_values))^2)
    ss_residual <- sum((actual_values - predicted_values)^2)
    r_squared <- 1 - (ss_residual / ss_total)

    mape <- mean(abs((actual_values - predicted_values) / actual_values)) * 100

    evaluation_results <- data.frame(
        # Metric = c("MSE", "RMSE", "MAE", "R_squared", "MAPE"),
        # Metric = c("MSE", "RMSE", "r_squared"),
        Metric = c("r_squared"),
        Value = c(r_squared)
        # Value = c(mse, rmse, r_squared)
        # Value = c(mse, rmse, mae, r_squared, mape)
    )

    return(evaluation_results)
}) |>
    bind_rows() |>
    group_by(Metric) |>
    summarise(
        r_squared_mean = mean(Value, na.rm = TRUE),
        r_squared_sd = sd(Value, na.rm = TRUE)
    ) |> 
    select(-Metric) |> 
    as.data.frame()


output_tsv <- cbind(metrics, counts) |> 
    relocate(physiology) |> 
    mutate(method = 'castor-ltp')



# plotList <- map2(hsp, testSets, ~ {
#     df <- data.frame(pred = .x, actual = .y)
#     df |>
#         ggplot(aes(pred, actual)) +
#         geom_point()
# })
# 
# p <- ggpubr::ggarrange(plotlist = plotList, ncol = 5, nrow = 2)


write.table(
    x = output_tsv, file = paste0(logfile, '.tsv'), sep = '\t', quote = FALSE,
    row.names = FALSE
)

# ggsave('inst/scripts/castor_gt.png', width = 12, height = 4)
log_close()







