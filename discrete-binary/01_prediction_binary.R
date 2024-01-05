
args <- commandArgs(trailingOnly = TRUE)
suppressMessages({
    library(phytools)
    library(castor)
    library(cvTools)
    library(taxPPro)
    library(bugphyzz)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(mltools)
    library(logr)
})

seed <- 42342
phys_name <- args[[1]]
rank_arg <- args[[2]]

if (rank_arg == 'all') {
    rank_var <- c('genus', 'species', 'strain')
} else {
    rank_var <- rank_arg
}

logfile <- paste0(phys_name, '_', rank_arg, '_phytools-ltp')
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

msg <- paste0('Predicting ', phys_name, '. Rank: ', rank_arg)
log_print(msg, blank_after = TRUE)

msg <- paste0('Strating at ', Sys.time())
log_print(msg, blank_after = TRUE)

## Tree data (taxPPro) ####
ltp <- ltp()
tree <- ltp$tree
tip_data <- ltp$tip_data
node_data <- ltp$node_data

## bugphyzz data
phys <- physiologies(gsub('_', ' ', phys_name))[[1]]
dat <- getDataReady(filterData(phys)) |> 
    filter(!is.na(Evidence)) |>
    filter(Rank %in% rank_var) 

fdat <- dat |> 
    filter(NCBI_ID %in% unique(tip_data$NCBI_ID)) |> 
    filter(!duplicated(NCBI_ID))

if (!nrow(fdat)) {
    msg <- paste0(
        'Not enough data for propagation of: ', phys_name, '; rank: ', rank_arg,
        'Quitting.'
    )
    log_print(msg)
    quit(save = 'no')
}

remove_ids <- split(fdat, fdat$Attribute) |> 
    keep(~ nrow(.x) < 10) |> 
    map(~ pull(.x, NCBI_ID)) |> 
    unlist() |> 
    unique()
fdat <- filter(fdat, !NCBI_ID %in% remove_ids)

## Create folds ####
keep_this <- split(fdat, fdat$Attribute)

if (length(keep_this) < 2) {
    msg <- paste0(
        'Not enough data for propagation of: ', phys_name, '; rank: ', rank_arg,
        'Quitting. Number of attributes: ', length(keep_this)
    )
    log_print(msg)
    quit(save = 'no')
}

set.seed(seed)
cv_folds <- map(keep_this, ~ cvFolds(n = nrow(.x), K = 10))

test_folds <- vector('list', 10)
train_folds <- vector('list', 10)
for (i in 1:10) {
    fold_name <- paste0('Fold', i)
    
    ## Code for creating test folds
    test_dfs <- vector('list', length(cv_folds))
    for (j in seq_along(cv_folds)) {
        test_dfs[[j]] <- keep_this[[j]][cv_folds[[j]]$which == i,]
    }
    test_folds[[i]] <- bind_rows(test_dfs)
    names(test_folds)[i] <- fold_name
    
    ## Code fro creating train folds
    train_dfs <- vector('list', length(cv_folds))
    for (j in seq_along(cv_folds)) {
        train_dfs[[j]] <- keep_this[[j]][cv_folds[[j]]$which != i,]
    }
    train_folds[[i]] <- bind_rows(train_dfs)
    names(train_folds)[i] <- fold_name
}

## Check that none of the ids in the test sets are in the train sets
## Value here must be TRUE
all(map2_lgl(test_folds, train_folds, ~ !any(.x$NCBI_ID %in% .y$NCBI_ID)))

counts <- data.frame(
    Attribute = gsub("_", " ", phys_name),
    ltp_bp_phys = sum(unique(tip_data$NCBI_ID) %in% unique(dat$NCBI_ID)),
    bp_phys = length(unique(dat$NCBI_ID)),
    ltp = Ntip(tree)

    )

write.table(
    x = counts, file = paste0(gsub(' ', '_', phys_name), '_', rank_arg, '_phytools-ltp_counts', '.tsv'),
    row.names = FALSE, quote = FALSE, sep = '\t'
)

## Create input matrix for tree ####
known_priors <- vector('list', length(train_folds))
for (i in seq_along(train_folds)) {
    names(known_priors)[i] <- paste0('Fold', i)
    known_priors[[i]] <- tip_data |> 
        left_join(
            train_folds[[i]], by = 'NCBI_ID', relationship = 'many-to-many'
        ) |> 
        filter(!is.na(Attribute)) |> 
        select(tip_label, Attribute, Score) |> 
        pivot_wider(
            names_from = 'Attribute', values_from = 'Score', values_fill = 0
        ) |> 
        tibble::column_to_rownames(var = 'tip_label') |> 
        as.matrix()
}

input_mats <- vector('list', length(known_priors))
for (i in seq_along(known_priors)) {
   tips <- tree$tip.label[!tree$tip.label %in% rownames(known_priors[[i]])] 
   input_mat <- matrix(
       data = rep(1/ncol(known_priors[[i]]), length(tips) * ncol(known_priors[[i]])),
       ncol = ncol(known_priors[[i]]),
       dimnames = list(rownames = tips, colnames = colnames(known_priors[[i]]))
   )
   input_mat <- rbind(input_mat, known_priors[[i]])
   input_mat[tree$tip.label]
   names(input_mats)[i] <- paste0('Fold', i)
   input_mats[[i]] <- input_mat
}

## Run prediction/propagation ####
predictions <- vector('list', length(input_mats))
for (i in seq_along(input_mats)) {
    fold_name <- paste0('Fold', i)
    message('Running predictions for ', fold_name)
    tim <- system.time({
        fit <- fitMk(
            tree = tree, x = input_mats[[i]], model = 'ER', pi = 'equal'
        )
        names(predictions)[[i]] <- fold_name
        
        
        predictions[[i]] <- ancr(object = fit, tips = TRUE)$ace
    })
    log_print(tim, blank_after = FALSE)
}

## ASR only predictions
predictions_tips <- vector('list', length(predictions))
for (i in seq_along(predictions)) {
    x <- predictions[[i]]
    x <- x[1:Ntip(tree),]
    names(predictions_tips)[[i]] <- paste0('Fold', i)
    predictions_tips[[i]] <- x |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = 'tip_label') |> 
        pivot_longer(
        names_to = 'Attribute', values_to = 'Score', cols = 2:last_col()
        ) |> 
        left_join(tip_data, by = 'tip_label')
}

msg <- paste0('Writing files')
log_print(msg, blank_after = TRUE)

## Write test sets
for (i in seq_along(test_folds)) {
    fname <- paste0(
        gsub(' ', '_', phys_name), '_', rank_arg, '_phytools-ltp_test_Fold', i, '.csv'
    )
    write.csv(
        x = test_folds[[i]], file = fname, quote = TRUE, row.names = FALSE
    )
}

## Write predicted sets
for (i in seq_along(predictions_tips)) {
    fname <- paste0(
        gsub(' ', '_', phys_name), '_', rank_arg, '_phytools-ltp_predicted_Fold', i, '.csv'
    )
    write.csv(
        x = predictions_tips[[i]], file = fname, quote = TRUE, row.names = FALSE
    )
}

msg <- paste0('Script finihsed at ', Sys.time())
log_print(msg, blank_after = FALSE)
log_close()
