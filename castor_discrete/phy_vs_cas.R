library(phytools)
library(castor)
library(cvTools)
library(taxPPro)
library(bugphyzz)
library(dplyr)
library(tidyr)
library(purrr)
library(mltools)

seed <- 42342

phys_name <- 'aerophilicity'
aer <- physiologies(phys_name)[[1]]
dat <- getDataReady(filterData(aer))

ltp <- ltp()
tip_data <- ltp$tip_data
tree <- ltp$tree

mat <- tip_data |> 
    left_join(dat, by = 'NCBI_ID', relationship = 'many-to-many') |> 
    select(tip_label, Attribute, Score) |> 
    filter(!is.na(Attribute)) |> 
    pivot_wider(names_from = 'Attribute', values_from = 'Score') |> 
    tibble::column_to_rownames(var = 'tip_label') |> 
    as.matrix()

set.seed(seed)
cv_folds <- cvFolds(n = nrow(mat), K = 10)
test_folds <- vector('list', 10)
train_folds <- vector('list', 10)

for (i in 1:10) {
    fold_name <- paste0('Fold', i)
    names(test_folds)[i] <- fold_name
    test_folds[[i]] <- mat[cv_folds$which == i,]
    names(train_folds)[i] <- fold_name
    train_folds[[i]] <- mat[cv_folds$which != i,]
}

myFun <- function(m, td = tip_data, s, tr = tree) {
    if (isTRUE(!s %in% c('phy', 'cas')))
        stop('Use cas or phy for s', call. = FALSE)
    row_names <- td |> 
        filter(!tip_label %in% rownames(m)) |> 
        pull(tip_label)
    if (s == 'phy') {
        data <- rep(1/ncol(m), length(row_names) * ncol(m))
    } else if (s == 'cas') {
        data <- rep(NA, length(row_names) * ncol(m))
    }
    output <- matrix(
        data = data, ncol = ncol(m),
        dimnames = list(rownames = row_names, colnames = colnames(m))
    )
    output <- rbind(output, m)
    output <- output[tr$tip.label,]
    return(output)
}

predicted_phy <- vector('list', 10)
predicted_cas <- vector('list', 10)

for (i in 1:10) {
    fold_name <- paste0('Fold', i)
    message('Predicting states for ', fold_name, ' with phytools.' )
    system.time({
        input_phys <- myFun(m = train_folds[[i]], s = 'phy')
        fit <- fitMk(tree = tree, x = input_phys, model = 'ER', pi = 'equal')
        res_phy <- ancr(fit, tips = TRUE)$ace
        rownames(res_phy) <- c(tree$tip.label, tree$node.label)
        predicted_phy[[i]] <- res_phy
        names(predicted_phy)[i] <- fold_name
    })
    
    message('Predicting states for ', fold_name, ' with castor.' )
    system.time({
        input_cas <- myFun(m = train_folds[[i]], s = 'cas')
        res_cas <- hsp_mk_model(
            tree = tree, tip_priors = input_cas, tip_states = NULL,
            check_input = FALSE, rate_model = 'ER', root_prior = 'flat'
        )$likelihoods
        rownames(res_cas) <- c(tree$tip.label, tree$node.label)
        predicted_cas[[i]] <- res_cas
        names(predicted_cas)[i] <- fold_name
    })
}

getPN <- function(predf, testf) {
    predf <- predf[rownames(testf),]
    predf[predf > 0.5] <- 'P' ## becomes character
    predf[predf != 'P'] <- 'N'
    pred <- predf |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = 'label') |> 
        pivot_longer(
            names_to = 'attribute', values_to = 'predPN', cols = 2:last_col()
        ) |> 
        arrange(label, attribute)
        
    testf[testf > 0.5] <- 'P'
    testf[testf != 'P'] <- 'N'
    test <- testf |> 
        as.data.frame() |> 
        tibble::rownames_to_column(var = 'label') |> 
        pivot_longer(
            names_to = 'attribute', values_to = 'testPN', cols = 2:last_col()
        ) |> 
        arrange(label, attribute)
    
    df <- left_join(test, pred, by = c('label', 'attribute')) |> 
        mutate(TF = case_when(
            testPN == 'P' & predPN == 'P' ~ 'TP',
            testPN == 'P' & predPN == 'N' ~ 'FN',
            testPN == 'N' & predPN == 'N' ~ 'TN',
            testPN == 'N' & predPN == 'P' ~ 'FP'
        ))
    
    l3 <- split(df, df$attribute) |> 
        map(~ {count(.x, TF)}) |> 
        discard(~ nrow(.x) < 4) |> 
        map( ~ {
            mcc(
                TP = pull(filter(.x, TF == 'TP'), n), 
                FP = pull(filter(.x, TF == 'FP'), n), 
                TN = pull(filter(.x, TF == 'TN'), n),
                FN = pull(filter(.x, TF == 'FN'), n),
            )
        })
    data.frame(
        attribute = names(l3),
        mcc = unlist(l3, use.names = FALSE)
    )
}

hh <- getPN(predicted_phy[[1]], test_folds[[1]])

phyRes <- vector('list', 10)
casRes <- vector('list', 10)

for (i in 1:10) {
    fold_name <- paste0('Fold', i)
    names(phyRes)[i] <- fold_name
    phyRes[[i]] <- getPN(predicted_phy[[i]], test_folds[[i]])
    names(casRes)[i] <- fold_name
    casRes[[i]] <- getPN(predicted_cas[[i]], test_folds[[i]])
}

phyRes |> 
    bind_rows() |> 
    group_by(attribute) |> 
    summarise(
        mean = mean(mcc),
        sd = sd(mcc)
    ) |> 
    ungroup()

casRes |> 
    bind_rows() |> 
    group_by(attribute) |> 
    summarise(
        mean = mean(mcc),
        sd = sd(mcc)
    ) |> 
    ungroup()
