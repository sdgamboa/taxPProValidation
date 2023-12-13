
# args <- commandArgs(trailingOnly = TRUE)
args <- list('width', 'all')

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
    library(VennDiagram)
    library(logr)
})

logfile <- gsub(" ", "_", paste0("logFile_", args[[1]], "_", args[[2]]))
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

attribute_group <- args[[1]]
rank_arg <- args[[2]]

msg <- paste0('Propagating ', attribute_group, '. Rank: ', rank_arg)
log_print(msg, blank_after = FALSE)

msg <- paste0('Starting at ', Sys.time())
log_print(msg, blank_after = TRUE)

# Import tree data --------------------------------------------------------
data('tree_list')
ncbi_tree <- as.Node(tree_list)
ncbi_gst_nodes <- ncbi_tree$Get(
    attribute = 'name',
    filterFun = function(node) {
        grepl('^[gst]__', node$name)
    },
    simplify = TRUE
) |> 
    unname()
ltp <- ltp3()
tree <- ltp$tree

# Define functions --------------------------------------------------------
getFolds <- function(dat, k_value = 10, seed = 1234) {
    k_value <- k_value
    set.seed(seed)
    cv_folds <- cvTools::cvFolds(nrow(dat), K = k_value)
    test_sets <- vector('list', k_value)
    train_sets <- vector('list', k_value)
    for(i in 1:k_value){
        fold_name <- paste0('Fold', i)
        names(test_sets)[i] <- fold_name
        test_sets[[i]] <- dat[cv_folds$subsets[cv_folds$which == i],]
        names(train_sets)[i] <- fold_name
        train_sets[[i]] <- dat[cv_folds$subsets[cv_folds$which != i],]
    }
    output <- list(test_sets = test_sets, train_sets = train_sets)
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

# Import and prepare bugphyzz data ----------------------------------------
if (rank_arg == 'all') {
    rank_var <- c('genus', 'species', 'strain')
} else {
    rank_var <- rank_arg
} 
bp_data <- physiologies(attribute_group)[[1]]
attribute_type <- unique(bp_data$Attribute_type)

if (attribute_type == 'range' && attribute_group %in% names(THRESHOLDS())) {
    message('Converting ', attribute_group, ' from range to multistate-intersection.')
    res <- rangeToLogicalThr(bp_data, THRESHOLDS()[[attribute_group]])
    res$Attribute_type <- 'multistate-intersection'
    bp_data <- res
    attribute_type <- unique(bp_data$Attribute_type)
} else if (attribute_type == 'range' && !attribute_group %in% names(THRESHOLDS())) {
    # quit(save = "no")
}

## Remove rows that will not be used
## The Rank column is removed at this step.
## This Rank column is generated again later in the script
## The main reason for doing this is that some taxids have ranks that do not
## correspond to current annotations in the NCBI
filtered_bp_data <- filterData(bp_data)

## There will be some warnings here from taxizedb::taxid2rank
## These can be ignored since they will be dealt with in the bugphyzzExports
## repository.
if (attribute_type == 'binary') {
    ## section for binary data ####
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
    ## Section for multistate-intersection ####
    set_with_ids <- filtered_bp_data |>
        getSetWithIDs() |>
        purrr::discard(~ all(is.na(.x))) |> 
        filter(Rank %in% rank_var)
    if (!nrow(set_with_ids)) {
        ## Only data from sources are supposed to be used for propagation.
        message('Not enough data for validation')
        quit(save = 'no')
    }
    folds <- set_with_ids |> 
        {\(y) split(y, factor(y$Attribute))}() |> 
        discard(~ nrow(.x) < 10) |> 
        purrr::map(getFolds)
    testFolds <- folds |> # exported later in the script
        names()|> 
        purrr::map(~ pluck(folds, .x, 'test_sets')) |> 
        set_names(names(folds)) |> 
        pmap(bind_rows)
    trainFolds <- folds |> 
        names() |>
        purrr::map(~ pluck(folds, .x, 'train_sets')) |> 
        set_names(names(folds)) |> 
        pmap(bind_rows)
    set_without_ids <- filtered_bp_data |> 
        getSetWithoutIDs(set_with_ids) |>
        purrr::discard(~ all(is.na(.x)))
    phys_data_ready <- trainFolds |> 
        map(~ bind_rows(.x, set_without_ids)) |> 
        map(~ complete(.x, NCBI_ID, Attribute, fill = list(Score = 0))) |> 
        map(~ arrange(.x, NCBI_ID, Attribute))
    
} else if (attribute_type == 'multistate-union') {
    ## Section for multistate-union data ####
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

propagated <- vector('list', length(phys_data_ready))
for (i in seq_along(propagated)) {
    dat <- phys_data_ready[[i]]
    attribute_nms <- dat$Attribute |>
        {\(y) y[!is.na(y)]}() |>
        unique()
    
    node_list <- split(dat, factor(dat$NCBI_ID))
    
    tipVar <- ltp$tip_data$NCBI_ID |> 
        {\(y) y[!is.na(y)]}() |> 
        unique()
    
    setsForVenn <- list(
        ncbiSet = ncbi_gst_nodes,
        ltpSet = tipVar,
        sourceSet = names(node_list)
    )
    vennFileName <- paste0(attribute_group, '_venn_befTax_Fold', i, '.png')
    venn_plot <- venn.diagram(
        x = setsForVenn, disable.logging = TRUE, filename = vennFileName,
        imagetype = "png"
    )
    
    ncbi_tree$Do(function(node) {
        if (node$name %in% names(node_list))
            node$attribute_tbl <- node_list[[node$name]]
    })
    
    # ncbi_tree$Do(
    #     function(node) {
    #         taxPool(
    #             node = node,
    #             grp = attribute_group,
    #             typ = attribute_type)
    #     },
    #     traversal = 'post-order'
    # )
    # 
    # tax_dat <- ncbi_tree$Get(
    #     'attribute_tbl', filterFun = function(node) {
    #         grepl('^[gst]__', node$name)
    #     }
    # ) |>
    #     discard(~ all(is.na(.x))) |>
    #     bind_rows() |>
    #     arrange(NCBI_ID, Attribute) |>
    #     filter(Evidence %in% c('tax'))
    
    # dat <- bind_rows(dat, tax_dat) 
    
    setsForVenn <- list(
        ncbiSet = ncbi_gst_nodes,
        ltpSet = tipVar,
        sourceSet = unique(dat$NCBI_ID)
    )
    vennFileName <- paste0(attribute_group, '_venn_aftTax_Fold', i, '.png')
    venn_plot <- venn.diagram(
        x = setsForVenn, disable.logging = TRUE, filename = vennFileName,
        imagetype = "png"
    )
    
    perTips <- round(mean(tipVar %in% unique(dat$NCBI_ID)) * 100)
    msg <- paste0('Percentage of tips annotated: ', perTips, '%')
    log_print(msg, blank_after = TRUE)
    
    perAnn <- round(mean(unique(dat$NCBI_ID) %in% tipVar) * 100)
    msg <- paste0('Percentage of annotated taxa used: ', perAnn, '%')
    log_print(msg, blank_after = TRUE)
    
    if (perTips < 1) {
        message('Not enough annotations for ASR for ', attribute_group)
        propagated[[i]] <- perTips
        names(propagated)[i] <- names(phys_data_ready)[i]
    }
    
    tip_data_annotated <- left_join(
        x = ltp$tip_data,
        y = select(dat, taxid, Attribute, Score),
        by = 'taxid', # joining by taxid helps to join cases like 561 and g__561
        relationship = 'many-to-many'
    )
    
    annotated_tips <- tip_data_annotated |>
        select(tip_label, Attribute, Score) |>
        filter(!is.na(Attribute)) |> 
        pivot_wider(
            names_from = 'Attribute', values_from = 'Score', values_fill = 0
        ) |>
        tibble::column_to_rownames(var = 'tip_label') |> 
        as.matrix()
    
    annotated_gn_tips <- grep('g__', rownames(annotated_tips), value = TRUE)
    remove_gn_tips <- ltp$gn_tips[!ltp$gn_tips %in% annotated_gn_tips]
    
    tip_data_subset <- ltp$tip_data |>
        filter(!tip_label %in% remove_gn_tips)
    
    if (attribute_type %in% c('binary', 'multistate-union')) {
        no_annotated_tips <- tip_data_subset |>
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
        
    } else if (attribute_type == 'multistate-intersection') {
        no_annotated_tip_names <- tip_data_subset |>
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
    
    tree <- ape::drop.tip(phy = tree, tip = remove_gn_tips)
    input_matrix <- rbind(annotated_tips, no_annotated_tips)
    input_matrix <- input_matrix[tree$tip.label,]
    
    
    models <- c(ER = 'ER', ARD = 'ARD', SYM = 'SYM')
    fittedModels <- map(models, ~ {
        message('Fitting ', .x, ' for ', attribute_group)
        start_time <- Sys.time()
        message(start_time)
        if (.x == 'ER') {
            r <- fitMk(
                tree = tree, x = input_matrix, model = .x,
                pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE
            )
        } else {
            r <- fitMk.parallel(
                tree = tree, x = input_matrix, model = .x,
                pi = 'fitzjohn', lik.func = 'pruning', logscale = TRUE,
                ncores = 10 
            )
        }
        end_time <- Sys.time()
        elapsed_time <- difftime(end_time, start_time)
        message('fitting model took like ', elapsed_time)
        return(r)
    })

    bestModel <- fittedModels |>
        map_dbl(AIC) |>
        aic.w() |>
        sort(decreasing = TRUE) |>
        names() |>
        {\(y) y[1]}()

    asr <- ancr(fittedModels[[bestModel]], tips = TRUE)
    
    res <- asr$ace |> 
        as.data.frame() |> 
        mutate(label = c(tree$tip.label, tree$node.label)) |>
        filter(!grepl('^n\\d+$', label)) |> 
        filter(label != 'NA') |> 
        filter(!label %in% rownames(annotated_tips)) # this removes all genus annotations in the tree. They are already in the sources
    res_tips_df <- res |>
        as.data.frame() |>
        filter(!grepl('^\\d+(\\+\\d+)*$', label)) |> 
        rename(tip_label = label)
    res_nodes_df <- res |>
        as.data.frame() |>
        filter(grepl('^\\d+(\\+\\d+)*$', label)) |> 
        rename(node_label = label)
    
    ## Get annotations for tips and nodes
    new_tips_data <- ltp$tip_data |>
        filter(tip_label %in% unique(res_tips_df$tip_label)) |>
        select(tip_label, taxid, Taxon_name, Rank) |>
        group_by(taxid) |>
        slice_max(order_by = tip_label, n = 1) |>
        ungroup() |>
        left_join(res_tips_df, by = 'tip_label') |>
        mutate(
            NCBI_ID = case_when(
                Rank == 'species' ~ paste0('s__', taxid),
                Rank == 'strain' ~ paste0('t__', taxid)
            )
        ) |>
        filter(Rank %in% c('species', 'strain')) |>
        mutate(Evidence = 'asr') |>
        relocate(NCBI_ID, taxid, Taxon_name, Rank, Evidence) |>
        pivot_longer(
            cols = 7:last_col(), names_to = 'Attribute', values_to = 'Score'
        ) |>
        mutate(
            Attribute_source = NA,
            Confidence_in_curation = NA,
            Attribute_group = attribute_group,
            Attribute_type = attribute_type,
            Frequency = case_when(
                Score == 1 ~ 'always',
                Score > 0.9 ~ 'usually',
                Score >= 0.5 ~ 'sometimes',
                Score > 0 & Score < 0.5 ~ 'rarely',
                Score == 0 ~ 'never'
            )
        ) |>
        select(-tip_label)
    
    new_nodes_data <- ltp$node_data |>
        filter(node_label %in% unique(res_nodes_df$node_label)) |>
        select(node_label, taxid, Taxon_name, Rank) |>
        left_join(res_nodes_df, by = 'node_label') |>
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
            Attribute_group = attribute_group,
            Attribute_type = attribute_type,
            Frequency = case_when(
                Score == 1 ~ 'always',
                Score > 0.9 ~ 'usually',
                Score >= 0.5 ~ 'sometimes',
                Score > 0 & Score < 0.5 ~ 'rarely',
                Score == 0 ~ 'never'
            )
        ) |>
        select(-node_label) |>
        filter(!NCBI_ID %in% dat$NCBI_ID)
    
    new_taxa_for_ncbi_tree <- bind_rows(new_tips_data, new_nodes_data) |>
        relocate(NCBI_ID, Rank, Attribute, Score, Evidence)
    
    new_taxa_for_ncbi_tree_list <- split(
        new_taxa_for_ncbi_tree, factor(new_taxa_for_ncbi_tree$NCBI_ID)
    )
    
    ## Add annotations from ASR 
    ncbi_tree$Do(function(node) {
        cond1 <- node$name %in% names(new_taxa_for_ncbi_tree_list)
        cond2 <- is.null(node$attribute_tbl) || all(is.na(node$attribute_tbl))
        if (cond1 && cond2) {
            node$attribute_tbl <- new_taxa_for_ncbi_tree_list[[node$name]]
        }
    })
    
    ## Inheritance comes from either tax or asr
    ncbi_tree$Do(
        fun = function(nd) inh1(nd, adjF = 0.1, evidence_label = 'inh'),
        traversal = 'pre-order'
    )
    
    # min_thr <- 1 / length(attribute_nms)
    
    inh_dat <- ncbi_tree$Get(
        attribute = 'attribute_tbl', simplify = FALSE,
        filterFun = function(node) {
            grepl('^[gst]__', node$name) && !is.null(node$attribute_tbl)
            
        }
    ) |>
        bind_rows() |>
        discard(~ all(is.na(.x))) |> 
        filter(Evidence == 'inh') |> 
        filter(!NCBI_ID %in% unique(dat$NCBI_ID))
    
    result <- bind_rows(
        dat, # contains source and tax (tax comes preparation)
        new_taxa_for_ncbi_tree, # only contains asr
        inh_dat  # contains inheritance
    )
    ncbi_tree$Do(cleanNode)
    names(propagated)[i] <- names(phys_data_ready)[i]
    propagated[[i]] <- result
}

lgl_vct <- map_lgl(propagated, is.data.frame)
if (all(!lgl_vct)) {
    msg <- paste0(
        'Not enough data for ASR for ', attribute_group, '. Rank:', rank_var, '.',
        ' Quitting R script. Finishing at ', Sys.time(), '.'
    )
    log_print(msg, blank_after = TRUE)
    quit(save = 'no')
}

propagated <- propagated[which(lgl_vct)]

msg <- paste0(
    'Saving test sets for ', attribute_group, ' ', rank_arg, '.'
)
log_print(msg)
for (i in seq_along(testFolds)) {
    fold_n <- names(testFolds)[i]
    fname <- paste0(attribute_group, '_test_', rank_arg, '_', fold_n, '.csv')
    fname <- gsub(" ", "_", fname)
    write.csv(
        x = testFolds[[i]], file = fname, row.names = FALSE,
        quote = TRUE
    )
}

msg <- paste0(
    'Saving training sets for ', attribute_group, ' ', rank_arg, '.'
)
log_print(msg)
for (i in seq_along(propagated)) {
    fold_n <- names(propagated)[i]
    fname <- paste0(attribute_group, '_propagated_', rank_arg, '_', fold_n, '.csv')
    fname <- gsub(" ", "_", fname)
    write.csv(
        x = propagated[[i]], file = fname, row.names = FALSE,
        quote = TRUE
    )
}

msg <- paste0(
    'Propagation workflow finished for ', attribute_group, '. Rank: ', rank_arg, '.',
    ' Finished at ', Sys.time(), '.'
)
log_print(msg, blank_after = TRUE)
log_close()
