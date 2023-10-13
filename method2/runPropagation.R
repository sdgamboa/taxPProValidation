## Scrupt for performing propagation and mcc
## This script should take a single physiology name only
## The script should take a second argument indicating the rank
## (all, genus, species, strain)

args <- commandArgs(trailingOnly = TRUE)

## Setup ####
library(logr)
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
library(mltools)
library(ggplot2)
library(forcats)

library(BiocParallel)

multicoreParam <- MulticoreParam(workers = 11)

getFolds10 <- function(dat, set_seed = 9013) {
    taxids <- unique(dat$NCBI_ID)
    set.seed(set_seed)
    test_index <- caret::createFolds(
        taxids, k = 10, list = TRUE, returnTrain = FALSE
    )
    test_taxids <- lapply(test_index, function(x) taxids[x]) 
    test_sets <- lapply(test_taxids, function(x) dat[dat$NCBI_ID %in% x,])
    
    set.seed(set_seed)
    train_index <- caret::createFolds(
        taxids, k = 10, list = TRUE, returnTrain = TRUE 
    )
    train_taxids <- lapply(train_index, function(x) taxids[x]) 
    train_sets <- lapply(train_taxids, function(x) dat[dat$NCBI_ID %in% x,])
    
    output <- list(
        test_sets = test_sets, train_sets = train_sets
    )
    return(output)
}

phys_name <- args[[1]]
rank <- args[[2]]
if (rank == 'all')
    rank <- c('genus', 'species', 'strain')

# phys_name <- 'aerophilicity'
# rank <- c('genus', 'species', 'strain')

bp_data <- physiologies(phys_name)[[1]]

at <- unique(bp_data$Attribute_type)
dat_name <- bp_data$Attribute_group
if (at == 'range' && dat_name %in% names(THRESHOLDS())) {
    res <- rangeToLogicalThr(bp_data, THRESHOLDS()[[dat_name]])
    res$Attribute_type <- 'multistate-intersection'
    bp_data <- res
} else if (at == 'range' && !dat_name %in% names(THRESHOLDS())) {
    quit(save = "no")
}

filtered_bp_data <- filterData(bp_data)

attr_type <- unique(filtered_bp_data$Attribute_type)
if (attr_type == 'binary') {
    # set_with_ids <- getSetWithIDs(tbl) |>
    #     purrr::discard(~ all(is.na(.x)))
    # set_without_ids <- getSetWithoutIDs(tbl, set_with_ids) |>
    #     purrr::discard(~ all(is.na(.x)))
    # dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
    # output <- completeBinaryData(dataset)
} else if (attr_type == 'multistate-intersection') {
    set_with_ids <- getSetWithIDs(filtered_bp_data) |>
        purrr::discard(~ all(is.na(.x)))
    folds <- getFolds10(set_with_ids)
    set_without_ids <- getSetWithoutIDs(
        filtered_bp_data, set_with_ids = set_with_ids
    ) |>
        purrr::discard(~ all(is.na(.x)))
   phys_data_ready <- map(folds$train_sets, ~ bind_rows(.x, set_without_ids)) |> 
       map(~ complete(.x, NCBI_ID, Attribute, fill = list(Score = 0))) |> 
       map(~ arrange(.x, NCBI_ID, Attribute))
} else if (attr_type == 'multistate-union') {
    # tbl$Attribute_group_2 <- sub('--(TRUE|FALSE)', '', tbl$Attribute)
    # l <- split(tbl, factor(tbl$Attribute_group_2))
    # output <- vector('list', length(l))
    # for (i in seq_along(output)) {
    #     set_with_ids <- getSetWithIDs(l[[i]]) |>
    #         purrr::discard(~ all(is.na(.x)))
    #     set_without_ids <- getSetWithoutIDs(l[[i]], set_with_ids) |>
    #         purrr::discard(~ all(is.na(.x)))
    #     dataset <- dplyr::bind_rows(set_with_ids, set_without_ids)
    #     # output[[i]] <- completeBinaryData(dataset)
    #     if (is.null(dataset))
    #         next
    #     names(output)[i] <- names(l)[i]
    #     output[[i]] <- completeBinaryData(dataset)
    # }
}

# phys_data_ready <- list_flatten(phys_data_ready)

## Prepare tree data ####
data('tree_list')
ncbi_tree <- as.Node(tree_list)

ltp <- ltp()
tree <- reorder(ltp$tree, 'postorder')
tip_data <- ltp$tip_data

tx <- grep('_taxid$', colnames(tip_data), value = TRUE)
nodes <- flatten(map(tx, ~ split(tip_data, factor(tip_data[[.x]]))))
nodes <- map(nodes, ~ .x[['tip_label']])
node_names <- map_int(nodes, ~ getMRCATaxPPro(tree, .x))
node_names <- node_names[!is.na(node_names)]
nodes_df <- data.frame(
    node = unname(node_names),
    node_label = names(node_names)
) |>
    group_by(node) |>
    mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
    ungroup() |>
    distinct()

start_time <- Sys.time()
propagated <- bplapply(
    X = phys_data_ready,
    BPPARAM = multicoreParam,
    FUN = function(dat)
    {
        Attribute_group_var <- unique(dat$Attribute_group)
        Attribute_group_var <- Attribute_group_var[!is.na(Attribute_group_var)]
        Attribute_type_var <- unique(dat$Attribute_type)
        Attribute_type_var <- Attribute_type_var[!is.na(Attribute_type_var)]
        dat_n_tax <- length(unique(dat$NCBI_ID))
        node_list <- split(x = dat, f = factor(dat$NCBI_ID))
        ncbi_tree$Do(function(node) {
            if (node$name %in% names(node_list))
                node$attribute_tbl <- node_list[[node$name]]
        })
        ncbi_tree$Do(
            function(node) {
                taxPool(
                    node = node,
                    grp = Attribute_group_var,
                    typ = Attribute_type_var)
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
            bind_rows(dat) # After this chunk is run, new_data also includes dat
        
        if (all(!new_dat$taxid %in% tip_data$taxid))
            return(new_dat)
        
        tip_data_annotated <- left_join(
            tip_data,
            select(new_dat, taxid, Attribute, Score),
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
        
        pruned_tree <- ape::keep.tip(tree, tip = rownames(annotated_tips))
        pruned_tree <- reorder(pruned_tree, 'postorder')
        pruned_tip_data <- tip_data |>
            filter(tip_label %in% pruned_tree$tip.label)
        pruned_node_data <- data.frame(
            node = length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
        )
        
        tx <- grep('_taxid$', colnames(pruned_tip_data), value = TRUE)
        nodes <- tx |>
            map(~ split(pruned_tip_data, factor(pruned_tip_data[[.x]]))) |>
            flatten() |>
            map(~ .x[['tip_label']])
        node_names <- map_int(nodes, ~ getMRCATaxPPro(pruned_tree, .x))
        node_names <- node_names[!is.na(node_names)]
        nodes_df <- data.frame(
            node = unname(node_names),
            node_label = names(node_names)
        ) |>
            group_by(node) |>
            mutate(node_label = paste0(unique(node_label), collapse = '+')) |>
            ungroup() |>
            distinct() |>
            arrange(node)
        pruned_node_data <- left_join(pruned_node_data, nodes_df, by = 'node') |>
            mutate(
                node_label = ifelse(
                    is.na(node_label), paste0('n', as.character(node)), node_label
                )
            )
        pruned_tree$node.label <- pruned_node_data$node_label
        
        fit <- fitMk(
            tree = pruned_tree, x = annotated_tips, model = 'ER',
            pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
        )
        asr <- ancr(object = fit, tips = TRUE)
        res <- asr$ace
        node_rows <- length(pruned_tree$tip.label) + 1:pruned_tree$Nnode
        rownames(res)[node_rows] <- pruned_tree$node.label
        
        nodes_annotated <- res[which(grepl('^\\d+(\\+\\d+)*', rownames(res))),]
        new_taxa_from_nodes <- nodes_annotated |>
            as.data.frame() |>
            tibble::rownames_to_column(var = 'NCBI_ID') |>
            filter(grepl('^\\d+(\\+\\d+)*', NCBI_ID)) |>
            mutate(NCBI_ID = strsplit(NCBI_ID, '\\+')) |>
            tidyr::unnest(NCBI_ID) |>
            mutate(Rank = taxizedb::taxid2rank(NCBI_ID)) |>
            mutate(Rank = ifelse(Rank == 'superkingdom', 'kingdom', Rank)) |>
            mutate(
                NCBI_ID = case_when(
                    Rank == 'kingdom' ~ paste0('k__', NCBI_ID),
                    Rank == 'phylum' ~ paste0('p__', NCBI_ID),
                    Rank == 'class' ~ paste0('c__', NCBI_ID),
                    Rank == 'order' ~ paste0('o__', NCBI_ID),
                    Rank == 'family' ~ paste0('f__', NCBI_ID),
                    Rank == 'genus' ~ paste0('g__', NCBI_ID),
                    Rank == 'species' ~ paste0('s__', NCBI_ID),
                    Rank == 'strain' ~ paste0('t__', NCBI_ID)
                )
            ) |>
            filter(
                Rank %in% c(
                    'kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                    'species', 'strain'
                )
            ) |>
            mutate(Evidence = 'asr') |>
            relocate(NCBI_ID, Rank, Evidence) |>
            pivot_longer(
                cols = 4:last_col(), names_to = 'Attribute', values_to = 'Score'
            ) |>
            mutate(
                Attribute_source = NA,
                Confidence_in_curation = NA,
                Attribute_group = Attribute_group_var,
                Attribute_type = Attribute_type_var,
                taxid = sub('\\w__', '', NCBI_ID),
                Taxon_name = taxizedb::taxid2name(taxid, db = 'ncbi'),
                Frequency = case_when(
                    Score == 1 ~ 'always',
                    Score > 0.9 ~ 'usually',
                    Score >= 0.5 ~ 'sometimes',
                    Score > 0 & Score < 0.5 ~ 'rarely',
                    Score == 0 ~ 'never'
                )
            )
        new_taxa_for_ncbi_tree <- new_taxa_from_nodes |>
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
        
        ncbi_tree$Do(
            function(node_var) {
                taxPool(
                    node = node_var,
                    grp = Attribute_group_var,
                    typ = Attribute_type_var
                )
            },
            traversal = 'post-order'
        )
        
        ncbi_tree$Do(inh2, traversal = 'pre-order')
        
        result <- ncbi_tree$Get(
            attribute = 'attribute_tbl', simplify = FALSE,
            filterFun = function(node) {
                node$name != 'ArcBac' && !is.null(node$attribute_tbl)
            }
        ) |>
            bind_rows() |>
            discard(~ all(is.na(.x)))
        
        # min_thr <- 1 / length(unique(dat$Attribute))
        add_taxa_1 <- dat |>
            filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
            discard(~ all(is.na(.x)))
        add_taxa_2 <- new_taxa_for_ncbi_tree |>
            filter(!NCBI_ID %in% unique(result$NCBI_ID)) |>
            discard(~ all(is.na(.x)))
        final_result <- bind_rows(list(result, add_taxa_1, add_taxa_2))
        # filter(Score > min_thr)
        ncbi_tree$Do(cleanNode)
        return(final_result)

    }
)

## Export test sets
for (i in seq_along(folds$test_sets)) {
    fold_n <- names(folds$test_sets)[i]
    fname <- paste0('method2/',phys_name, '_test_', fold_n, '.csv')
    write.csv(
        x = folds$test_sets[[i]], file = fname, row.names = FALSE,
        quote = TRUE
    )
}

## Export propagated sets
for (i in seq_along(propagated)) {
    fold_n <- names(propagated)[i]
    fname <- paste0(phys_name, '_propagated_', fold_n, '.csv')
    write.csv(
        x = propagated[[i]], file = fname, row.names = FALSE,
        quote = TRUE
    )
}











