# args <- commandArgs(trailingOnly = TRUE)
args <- list('acetate producing', 'all')

suppressMessages({
    library(bugphyzz)
    library(taxPPro)
    library(data.tree)
    library(phytools)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(ggplot2)
    library(ape)
    library(bugphyzzExports)
    library(cvTools)
    library(ggplot2)
    library(forcats)
    library(BiocParallel)
})

n_cores <- parallel::detectCores()
if (n_cores <= 16) {
    n_cores <- 5
} else {
    n_cores <- 10
}

multicoreParam <- MulticoreParam(workers = n_cores)

data('tree_list')
ncbi_tree <- as.Node(tree_list)
ncbi_nodes <- ncbi_tree$Get(
    attribute = 'name', filterFun = function(node) grepl('^[gst]__', node$name),
    simplify = TRUE
) |> 
    unname()

getFolds10 <- function(dat, k_value = 10) {
    taxa <- unique(dat$NCBI_ID)
    k_value <- 10
    set.seed(1234)
    taxa_folds <- cvFolds(length(taxa), K = k_value)
    test_taxids <- vector('list', k_value)
    train_taxids <- vector('list', k_value)
    for(i in 1:k_value){
        fold_name <- paste0('Fold', i)
        names(test_taxids)[i] <- fold_name
        test_taxids[[i]] <- taxa[taxa_folds$subsets[taxa_folds$which == i]]
        names(train_taxids)[i] <- fold_name
        train_taxids[[i]] <- taxa[taxa_folds$subsets[taxa_folds$which != i]]
    }
    test_sets <- purrr::map(
        test_taxids, ~ dplyr::filter(dat, .data$NCBI_ID %in% .x)
    )
    train_sets <- purrr::map(
        train_taxids, ~ dplyr::filter(dat, .data$NCBI_ID %in% .x)
    )
    output <- list(
        test_sets = test_sets,
        train_sets = train_sets
    )
    return(output)
}

getNegatives <- function(dat) {
    ranks <- unique(dat$Rank)
    atgrp <- unique(dat$Attribute_group)
    attyp <- unique(dat$Attribute_type)
    ranks <- factor(ranks, levels = c('genus', 'species', 'strain'))
    ranks <- case_when(
        ranks == 'genus' ~ 'g__',
        ranks == 'species' ~ 's__',
        ranks == 'strain' ~ 't__'
    ) |> 
        {\(y) paste0('^(', paste0(y, collapse = '|'), ')')}()
    selected_nodes <- ncbi_nodes[which(grepl(ranks, ncbi_nodes))]
    l <- split(dat, dat$Attribute)
    if (length(l) == 2) {
        trues <- l[which(grepl('TRUE', names(l)))][[1]]
        falses <- l[which(grepl('FALSE', names(l)))][[1]]
        
        if (nrow(trues) > nrow(falses)) {
            ncbi_ids <- unique(dat$NCBI_ID)
            n <- length(ncbi_ids)
            selected_nodes <- selected_nodes[!selected_nodes %in% ncbi_ids]
            set.seed(20308)
            negatives_imputation <- sample(selected_nodes, n) 
            output <- data.frame(
                NCBI_ID = negatives_imputation,
                Attribute = sub('--TRUE$', '--FALSE', dat$Attribute),
                Score = 1
            ) |> 
                mutate(
                    taxid = sub('\\w__', '', NCBI_ID),
                    Taxon_name = taxizedb::taxid2name(taxid),
                    Rank = taxizedb::taxid2rank(taxid),
                    Attribute_group = atgrp,
                    Attribute_type = attyp
                )
            
        } else {
            return(NULL)
        }
        
    } else {
        ncbi_ids <- unique(dat$NCBI_ID)
        n <- length(ncbi_ids)
        selected_nodes <- selected_nodes[!selected_nodes %in% ncbi_ids]
        set.seed(20308)
        negatives_imputation <- sample(selected_nodes, n) 
        output <- data.frame(
            NCBI_ID = negatives_imputation,
            Attribute = sub('--TRUE$', '--FALSE', dat$Attribute),
            Score = 1
        ) |> 
            mutate(
                taxid = sub('\\w__', '', NCBI_ID),
                Taxon_name = taxizedb::taxid2name(taxid),
                Rank = taxizedb::taxid2rank(taxid),
                Attribute_group = atgrp,
                Attribute_type = attyp
            )
    }
    return(output)
}

phys_name <- args[[1]]
rank_var <- args[[2]]
if (rank_var == 'all') {
    rank_var <- c('genus', 'species', 'strain')
}
suppressMessages({
    bp_data <- physiologies(phys_name)[[1]]
})

attribute_type <- unique(bp_data$Attribute_type)
attribute_group <- unique(bp_data$Attribute_group)
if (attribute_type == 'range' && attribute_group %in% names(THRESHOLDS())) {
    res <- rangeToLogicalThr(bp_data, THRESHOLDS()[[attribute_group]])
    res$Attribute_type <- 'multistate-intersection'
    bp_data <- res
    attribute_type <- unique(bp_data$Attribute_type)
} else if (attribute_type == 'range' && !attribute_group %in% names(THRESHOLDS())) {
    quit(save = "no")
}

filtered_bp_data <- filterData(bp_data)

if (attribute_type == 'binary') {
    set_with_ids <- filtered_bp_data |> 
        getSetWithIDs() |>
        purrr::discard(~ all(is.na(.x))) |> 
        filter(Rank %in% rank_var)
    if (!nrow(set_with_ids)) {
        message('Not enough data for validation')
        quit(save = 'no')
    }
    folds <- getFolds10(set_with_ids)
    folds$test_sets <- map(folds$test_sets, ~ bind_rows(.x, getNegatives(.x)))
    set_without_ids <- getSetWithoutIDs(
        filtered_bp_data, set_with_ids = set_with_ids
    ) |>
        purrr::discard(~ all(is.na(.x)))
   phys_data_ready <- map(folds$train_sets, ~ bind_rows(.x, set_without_ids)) |> 
       map(~ completeBinaryData(.x)) |> 
       map(~ arrange(.x, NCBI_ID, Attribute))
   
} else if (attribute_type == 'multistate-intersection') {
    set_with_ids <- filtered_bp_data |>
        getSetWithIDs() |>
        purrr::discard(~ all(is.na(.x))) |> 
        filter(Rank %in% rank_var)
    if (!nrow(set_with_ids)) {
        message('Not enough data for validation')
        quit(save = 'no')
    }
    folds <- getFolds10(set_with_ids)
    set_without_ids <- getSetWithoutIDs(
        filtered_bp_data, set_with_ids = set_with_ids
    ) |>
        purrr::discard(~ all(is.na(.x)))
   phys_data_ready <- map(folds$train_sets, ~ bind_rows(.x, set_without_ids)) |> 
       map(~ complete(.x, NCBI_ID, Attribute, fill = list(Score = 0))) |> 
       map(~ arrange(.x, NCBI_ID, Attribute))
   
} else if (attr_type == 'multistate-union') {
    filtered_bp_data$Attribute_group_2 <- sub(
        '--(TRUE|FALSE)', '', filtered_bp_data$Attribute
    )
    l <- split(filtered_bp_data, factor(filtered_bp_data$Attribute_group_2))
    sets_with_ids <- vector('list', length(l))
    for (i in seq_along(sets_with_ids)) {
        names(sets_with_ids)[i] <- names(l)[i]
        suppressWarnings({
            res <- getSetWithIDs(l[[i]]) |>
                purrr::discard(~ all(is.na(.x))) |> 
                filter(Rank %in% rank_var)
        })
        if (length(res) > 0)
            sets_with_ids[[i]] <- res
    }
    select_names <- 
        names(head(sort(map_int(sets_with_ids, nrow), decreasing = TRUE), 3))
    sets_with_ids <- sets_with_ids[select_names] 
    folds <- purrr::map(sets_with_ids, getFolds10)
    for (i in seq_along(folds)) {
        folds[[i]]$test_sets <- purrr::map(folds[[i]]$test_sets, ~ bind_rows(.x, getNegatives(.x)))
        
    }
    test_sets <- purrr::map(folds, ~ .x[['test_sets']])
    train_sets <- purrr::map(folds, ~ .x[['train_sets']])
    l2 <- l[select_names]
    phys_data_ready <- vector('list', length(l2))
    for (i in seq_along(select_names)) {
        res <- map(train_sets[[i]], ~ getSetWithoutIDs(l2[[i]], .x)) |> 
            map(~ purrr::discard(.x, ~ all(is.na(.x))))
        datasets <- map2(.x = train_sets[[i]], .y = res, bind_rows)
        if (all(is.null(datasets)))
            next
        names(phys_data_ready)[i] <- names(l2)[i]
        phys_data_ready[[i]] <- map(datasets, completeBinaryData)
    }
    phys_data_ready <- discard(phys_data_ready, is.null)
    phys_data_ready <- list_flatten(phys_data_ready)
    test_sets <- list_flatten(test_sets)
    test_sets <- test_sets[names(phys_data_ready)]
    folds <- list(test_sets = test_sets) # convert to list just to make it compatible with export code (below)
}

suppressMessages({
    ltp <- ltp()
})
tree <- reorder(ltp$tree, 'postorder')
tip_data <- ltp$tip_data
node_data <- ltp$node_data

propagated <- bplapply(
    X = phys_data_ready,
    BPPARAM = multicoreParam,
    FUN = function(dat) {
        attribute_nms <- unique(dat$Attribute) |>
            {\(y) y[!is.na(y)]}()
        
        dat_n_tax <- length(unique(dat$NCBI_ID))
        node_list <- split(dat, factor(dat$NCBI_ID))
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
        ncbi_tree$Do(
            function(node) {
                taxPool(
                    node = node,
                    grp = attribute_group,
                    typ = attribute_type)
            },
            traversal = 'post-order'
        )
        ncbi_tree$Do(inh1, traversal = 'pre-order')
        
        new_dat <- ncbi_tree$Get(
            'attribute_tbl', filterFun = function(node) {
                grepl('^[gst]__', node$name)
            }
        ) |>
            discard(~ all(is.na(.x))) |>
            bind_rows() |>
            arrange(NCBI_ID, Attribute) |>
            filter(!NCBI_ID %in% dat$NCBI_ID) |>
            bind_rows(dat)
        
        new_taxids <- new_dat |> 
            pull(taxid) |> 
            {\(y) y[!is.na(y)]}()
        per <- mean(tip_data$taxid %in% new_taxids) * 100
        if (per < 1) {
            return(NULL)
        }
        
        tip_data_annotated <- left_join(
            x = tip_data,
            y = select(new_dat, taxid, Attribute, Score),
            by = 'taxid'
        )
        
        annotated_tips <- tip_data_annotated |>
            select(tip_label, Attribute, Score) |>
            filter(!is.na(Attribute)) |>
            pivot_wider(
                names_from = 'Attribute', values_from = 'Score', values_fill = 0
            ) |>
            tibble::column_to_rownames(var = 'tip_label') |>
            as.matrix()
        
        if (attribute_type %in% c('binary', 'multistate-union')) {
            no_annotated_tips <- tip_data |>
                filter(!tip_label %in% rownames(annotated_tips)) |>
                select(tip_label) |>
                mutate(
                    Attribute = factor(
                        attribute_nms[[1]], levels = attribute_nms
                    )
                ) |>
                complete(tip_label, Attribute) |>
                mutate(Score = ifelse(grepl('--FALSE$', Attribute), 1, 0)) |>
                pivot_wider(names_from = 'Attribute', values_from = 'Score') |>
                tibble::column_to_rownames(var = 'tip_label') |>
                as.matrix() |>
                {\(y) y[,colnames(annotated_tips)]}()
            
        } else if (attribute_type_var == 'multistate-intersection') {
            no_annotated_tip_names <- tip_data |>
                filter(!tip_label %in% rownames(annotated_tips)) |>
                pull(tip_label)
            fill_value <- 1 / length(attribute_nms)
            vct <- rep(
                fill_value,
                length(no_annotated_tip_names) * length(attribute_nms)
            )
            no_annotated_tips <- matrix(
                data = vct,
                nrow = length(no_annotated_tip_names),
                ncol = length(attribute_nms),
                dimnames = list(no_annotated_tip_names, attribute_nms)
            )
        }
        
        input_matrix <- rbind(annotated_tips, no_annotated_tips)
        input_matrix <- input_matrix[tree$tip.label,]
        
        fit <- fitMk(
            tree = tree, x = input_matrix, model = 'ER',
            pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
        )
        asr <- ancr(object = fit, tips = TRUE)
        res <- asr$ace
        node_rows <- length(tree$tip.label) + 1:tree$Nnode
        rownames(res)[node_rows] <- tree$node.label
        res <- res[tree$node.label,]
        res_df <- res |>
            as.data.frame() |>
            tibble::rownames_to_column(var = 'node_label') |>
            filter(!grepl('^n\\d+', node_label))
        
        node_data_annotated <- node_data |>
            filter(node_label %in% unique(res_df$node_label)) |>
            select(node_label, taxid, Taxon_name, Rank)
        
        nodes_annotated <- node_data_annotated |>
            left_join(res_df, by = 'node_label') |>
            mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
            mutate(
                NCBI_ID = case_when(
                    Rank == 'kingdom' ~ paste0('k__', taxid),
                    Rank == 'phylum' ~ paste0('p__', taxid),
                    Rank == 'class' ~ paste0('c__', taxid),
                    Rank == 'order' ~ paste0('o__', taxid),
                    Rank == 'family' ~ paste0('f__', taxid),
                    Rank == 'genus' ~ paste0('g__', taxid),
                    Rank == 'species' ~ paste0('s__', taxid),
                    Rank == 'strain' ~ paste0('t__', taxid)
                )
            ) |>
            filter(
                Rank %in% c(
                    'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                    'species', 'strain'
                )
            ) |>
            mutate(Evidence = 'asr') |>
            relocate(NCBI_ID, taxid, Taxon_name, Rank, Evidence) |>
            pivot_longer(
                cols = 7:last_col(), names_to = 'Attribute', values_to = 'Score'
            ) |>
            mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                Attribute_group = Attribute_group_var,
                Attribute_type = Attribute_type_var,
                Frequency = case_when(
                    Score == 1 ~ 'always',
                    Score > 0.9 ~ 'usually',
                    Score >= 0.5 ~ 'sometimes',
                    Score > 0 & Score < 0.5 ~ 'rarely',
                    Score == 0 ~ 'never'
                )
            )
        
        new_taxa_for_ncbi_tree <- nodes_annotated |>
            relocate(NCBI_ID, Rank, Attribute, Score, Evidence)
        new_taxa_for_ncbi_tree_list <- split(
            new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
        )
        
        ncbi_tree$Do(function(node) {
            cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
            cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
            if (cond1 && cond2) {
                node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
            }
        })
        
        ncbi_tree$Do(inh2, traversal = 'pre-order')
        
        result <- ncbi_tree$Get(
            attribute = 'attribute_tbl', simplify = FALSE,
            filterFun = function(node) {
                node$name != 'ArcBac' && !is.null(node$attribute_tbl)
            }
        ) |>
            bind_rows() |>
            discard(~ all(is.na(.x)))
        
        add_taxa_1 <- dat |>
            filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
            discard(~ all(is.na(.x)))
        add_taxa_2 <- new_taxa_for_ncbi_tree |>
            filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
            discard(~ all(is.na(.x)))
        final_result <- bind_rows(list(result, add_taxa_1, add_taxa_2))
        ncbi_tree$Do(cleanNode)
        return(final_result)
    }
)

if (all(map_lgl(propagated, is.null))) {
    message('Not enough data for ASR. Stopping after taxonomic pooling.')
    quit(save = 'no')
}

if (all(c('genus', 'species', 'strain') %in% rank_var)) {
    rank_var <- 'all'
}

for (i in seq_along(folds$test_sets)) {
    fold_n <- names(folds$test_sets)[i]
    fname <- paste0(phys_name, '_test_', rank_var, '_', fold_n, '.csv')
    fname <- gsub(" ", "_", fname)
    write.csv(
        x = folds$test_sets[[i]], file = fname, row.names = FALSE,
        quote = TRUE
    )
}

for (i in seq_along(propagated)) {
    fold_n <- names(propagated)[i]
    fname <- paste0(phys_name, '_propagated_', rank_var, '_', fold_n, '.csv')
    fname <- gsub(" ", "_", fname)
    write.csv(
        x = propagated[[i]], file = fname, row.names = FALSE,
        quote = TRUE
    )
}
