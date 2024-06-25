
args <- commandArgs(trailingOnly = TRUE)
# args <- list('growth_temperature', 'all')

library(taxPPro)
library(dplyr)
library(bugphyzz)
library(phytools)
library(castor)
library(cvTools)
library(purrr)
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

phys <- physiologies(gsub('_', ' ', phys_name))[[1]]

dat <- phys |> 
    filterData() |> 
    getDataReady() |> 
    filter(Rank %in% rank_var)

if (!nrow(dat)) {
    log_print('not enough data for propagation')
    quit(save = 'no')
    log_close()
}
    
fdat <- dat |> 
    filter(Rank %in% rank_var) |>
    filter(NCBI_ID %in% unique(tip_data$NCBI_ID))

if (!nrow(fdat)) {
    log_print('not enough data for propagation')
    quit(save = 'no')
    log_close()
}

tip_data_annotated <- left_join(tip_data, fdat, by = c('NCBI_ID')) |>
    select(tip_label, Attribute_value) # I don't need to filter NAs here

known_tips <- tip_data_annotated |> 
    filter(!is.na(Attribute_value)) |> 
    pull(tip_label)

nodes_with_taxid <- grep("^n", tree$node.label, value = TRUE, invert = TRUE)
nsti <- getNsti(tree, known_tips, nodes_with_taxid)

counts <- data.frame(
    physiology = gsub('_', ' ', phys_name),
    rank = rank_arg,
    ltp_bp_phys = length(unique(fdat$NCBI_ID)),
    bp_phys = length(unique(dat$NCBI_ID)),
    ltp = length(unique(tip_data$NCBI_ID)),
    nsti_mean = round(mean(nsti$nsti), 3),
    nsti_sd = round(sd(nsti$nsti), 3)
)

total_per <- floor(length(unique(fdat$NCBI_ID)) / Ntip(tree) * 100)
if (total_per < 1) {
    message("Not enough percentage data for propagation.")
    quit(save = "no")
}

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
        tree = tree, tip_states = input_vector, weighted = TRUE,
        check_input = TRUE
    )
    
    statesDF <- data.frame(
        label = tree$tip.label,
        value = res$states[1:Ntip(tree)]
    )

    states <- statesDF$value
    names(states) <- statesDF$label

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
    # r_squared <- 1 - (ss_residual / ss_total)
    r_squared <- summary(lm(formula = actual_values ~ predicted_values))$r.squared

    mape <- mean(abs((actual_values - predicted_values) / actual_values)) * 100

    evaluation_results <- data.frame(
        Metric = c("r_squared"),
        Value = c(r_squared)
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

write.table(
    x = output_tsv, file = paste0(logfile, '.tsv'), sep = '\t', quote = FALSE,
    row.names = FALSE
)

log_close()
