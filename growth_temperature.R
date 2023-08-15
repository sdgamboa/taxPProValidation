
library(logr)
library(bugphyzz) # BiocManager::install('waldronlab/bugphyzz', force = TRUE)
library(taxPPro) # BiocManager::install('sdgamboa/taxPPro', force = TRUE)
library(purrr)
library(rlang)
library(dplyr)
library(data.tree)
library(bugphyzzExports) # Install locally
library(BiocParallel)
library(tidyr)

## Start log
phys_name <- 'growth temperature'
logfile <- gsub(' ', '_', paste0(phys_name, "_log_file"))
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)

## Select number of cores
n_threads <- parallel::detectCores()
if (n_threads > 16) {
    n_threads <- round(n_threads * 0.6)
}
msg <- paste0('Using ', n_threads, ' cores.')
log_print(msg, blank_after = TRUE)

## Keep both TRUE and FALSE values for these physiologies
binaries <- c(
    "acetate producing",
    "animal pathogen",
    "antimicrobial sensitivity",
    "biofilm forming",
    "butyrate producing",
    "extreme environment",
    "health associated",
    "hydrogen gas producing",
    "lactate producing",
    "motility",
    "pathogenicity human",
    "plant pathogenicity",
    "sphingolipid producing",
    "spore formation"
)

## Import data
msg <- paste0('Creating data to validate ', phys_name, '.')
log_print(msg, blank_after = TRUE)
phys <- physiologies(phys_name, full_source = FALSE)
phys <- map(phys, ~ {
    attr_grp <- unique(.x$Attribute_group)
    if (attr_grp %in% binaries) {
        .x$Attribute <- paste0(.x$Attribute, '--', .x$Attribute_value)
        .x$Attribute_value <- TRUE
    }
    if (unique(.x$Attribute_type) == 'logical') {
        .x <- filter(.x, Attribuje_value == TRUE)
    }
    if ('Unit' %in% colnames(.x)) {
        ## fill empty units
        unit <- .x$Unit
        unit <- unique(unit[!is.na(unit)])
        .x$Unit <- unit
    }
    return(.x)
})

## Convert all attributes to categorical before propagation
categorical <- keep(phys, ~ unique(.x$Attribute_type) == 'logical')
range <- keep(phys, ~ unique(.x$Attribute_type == 'range'))
range <- range[which(names(range) %in% names(THRESHOLDS()))]
range_cat <- map2(range, names(range), ~ rangeToLogicalThr(.x, THRESHOLDS()[[.y]]))
categorical <- c(categorical, range_cat)

## Include <Attribute>--TRUE and <Attribute>--FALSE are added for binary attributes.
fname <- system.file('extdata', 'attributes.tsv', package = 'bugphyzz')
valid_attributes <- unique(read.table(fname, header = TRUE, sep = '\t')$attribute)
more_valid_attributes <- map(categorical, ~ {
    attr_grp <- unique(.x$Attribute_group)
    if (attr_grp %in% binaries) {
        binary_attr <- unique(.x$Attribute)
        return(binary_attr)
    }
}) |>
    discard(is.null) |>
    unlist(recursive = TRUE, use.names = FALSE)
valid_attributes <- unique(c(valid_attributes, more_valid_attributes))

## Keep only attributes with valid values
data <- map(categorical, ~ {
    attr_names <- unique(.x$Attribute)
    attr_grp <- unique(.x$Attribute_group)
    lgl <- sum(!attr_names %in% valid_attributes)
    if (lgl > 0) {
        invalid_values <- filter(.x, !Attribute %in% valid_attributes)
        invalid_values <- invalid_values |>
            select(Attribute_group, Attribute) |>
            unique() |>
            as_tibble()
    }
    output <- filter(.x, Attribute %in% valid_attributes)
    return(output)
}) |>
    discard(~ !nrow(.x))

data_discarded <- names(categorical)[which(!names(categorical) %in% names(data))]
if (length(data_discarded) > 0) {
    data_discarded <- paste0(data_discarded, collapse = ', ')
    msg <- paste0(
        "The following physiologies were discarded because they didn't have",
        " any valid Attribute: ", data_discarded, '.'
    )
}

## This chunk is to make sure that all NCBI_IDs have valid taxon names
data <- bplapply(data, BPPARAM = MulticoreParam(workers = n_threads), FUN = function(x) {
    set1 <- filter(x, !is.na(NCBI_ID))
    set2 <- filter(x, is.na(NCBI_ID))
    set1$Rank <- checkRank(set1$NCBI_ID)
    set1$Taxon_name <- checkTaxonName(set1$NCBI_ID)
    bind_rows(set1, set2)
})

## For those numeric attributes converted to categorical attributes,
## propagate them with their ranges and units.
data <- map(data, ~ {
    if ('Attribute_range' %in% colnames(.x)) {
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Attribute_range)
        .x$Attribute <- sub('\\)$', '', .x$Attribute)
        .x$Attribute <- paste0(.x$Attribute, ' ', .x$Unit, ')')
        .x$Attribute <- sub(' \\)', ')', .x$Attribute)
    }
    return(.x)
})

## Prepapre data for propagation
data_ready <- bplapply(data, BPPARAM = MulticoreParam(workers = n_threads), FUN = function(x) {
    tryCatch(
        error = function(e) e,
        {
            output <- x |>
                prepareDataForPropagation() |>
                mergeOriginalAndEarlyASR() |>
                group_by(NCBI_ID) |>
                mutate(Score = Score / sum(Score)) |>
                ungroup() |>
                distinct()
            return(output)
        }
    )
})

lgl <- map_lgl(data_ready, is_error)
if (any(lgl)) {
    discard_names <- names(data_ready)[which(lgl)]
    msg <- paste0(
        "The following physiologies will be discared for some error during",
        " preparation for propagation: ",
        paste0("'", paste0(discard_names, collapse = ', '), ".'")
    )

    data_ready <- discard(data_ready, is_error)
}

## taxids in data.tree structure
data('tree_list')
tree <- as.Node(tree_list)
tree_ids <- unname(tree$Get(function(node) node$name))

## I only want to compare  with taxids that can actually be found in the tree
myDF <- data_ready[[phys_name]]
taxids <- unique(myDF$NCBI_ID)
early_asr_taxids <- myDF |>
    
    filter(Evidence == 'asr') |> 
    pull(NCBI_ID)

holdout_options <- c(
    'holdout_all', 'holdout_gn', 'holdout_sp', 'holdout_st'
)
per <- 0.3

holdout_subsets <- vector('list', length(holdout_options))
for (i in seq_along(holdout_subsets)) {
    if (holdout_options[i] == 'holdout_all') {
        valid_taxids <- taxids[taxids %in% tree_ids & !taxids %in% early_asr_taxids]
        set.seed(1234)
        size_all <- round(length(valid_taxids) * per)
        holdout_taxids <- sample(x = valid_taxids, size = size_all, replace = FALSE)
        holdouts_df <- myDF |> 
            filter(NCBI_ID %in% holdout_taxids)
        if (!nrow(holdouts_df)) 
            next
        holdout_subsets[[i]] <- holdouts_df
        names(holdout_subsets)[i] <- 'holdout_all'
        datName <- gsub(' ', '_', paste0(phys_name, '_prediction_all'))
        data_ready[[datName]] <- myDF |> 
            filter(!NCBI_ID %in% holdout_taxids)
    }
    if (holdout_options[i] == 'holdout_gn') {
        regex_g <- '^g__'
        gn_taxids <- taxids[which(grepl(regex_g, taxids) & taxids %in% tree_ids & !taxids %in% early_asr_taxids)]
        set.seed(1234)
        size_gn <- round(length(gn_taxids) * per)
        holdout_gn_taxids <- sample(x = gn_taxids, size = size_gn, replace = FALSE)
        holdout_subsets[[i]] <- myDF |> 
            filter(NCBI_ID %in% holdout_gn_taxids)
        names(holdout_subsets)[i] <- 'holdout_gn'
        datName <- gsub(' ', '_', paste0(phys_name, '_prediction_gn'))
        data_ready[[datName]] <- myDF |> 
            filter(!NCBI_ID %in% holdout_gn_taxids)
    }
    if (holdout_options[i] == 'holdout_sp') {
        regex_s <- '^s__'
        sp_taxids <- taxids[which(grepl(regex_s, taxids) & taxids %in% tree_ids & !taxids %in% early_asr_taxids)]
        set.seed(1234)
        size_sp <- round(length(sp_taxids) * per)
        holdout_sp_taxids <- sample(x = sp_taxids, size = size_sp, replace = FALSE)
        holdout_subsets[[i]] <- myDF |> 
            filter(NCBI_ID %in% holdout_sp_taxids)
        names(holdout_subsets)[i] <- 'holdout_sp'
        datName <- gsub(' ', '_', paste0(phys_name, '_prediction_sp'))
        data_ready[['phys_prediction_sp']] <- myDF |> 
            filter(!NCBI_ID %in% holdout_sp_taxids)
    }
    if (holdout_options[i] == 'holdout_st') {
        regex_t <- '^t__'
        st_taxids <- taxids[which(grepl(regex_t, taxids) & taxids %in% tree_ids & !taxids %in% early_asr_taxids)]
        set.seed(1234)
        size_st <- round(length(st_taxids) * per)
        holdout_st_taxids <- sample(x = st_taxids, size = size_st, replace = FALSE)
        holdout_subsets[[i]] <- myDF |> 
            filter(NCBI_ID %in% holdout_st_taxids)
        names(holdout_subsets)[i] <- 'holdout_st'
        datName <- gsub(' ', '_', paste0(phys_name, '_prediction_st'))
        data_ready[[datName]] <- myDF |> 
            filter(!NCBI_ID %in% holdout_st_taxids)
    }
}
names(holdout_subsets) <- gsub(' ', '_', paste0(phys_name, '_', names(holdout_subsets)))

positions <- which(map_lgl(holdout_subsets, ~ !nrow(.x)))
names(holdout_subsets)[positions]



for (i in seq_along(holdout_subsets)) {
    write.table(
        x = holdout_subsets[[i]],
        file = paste0(names(holdout_subsets)[i], '.csv'),
        quote = TRUE, 
        sep = ',', 
        row.names = FALSE
    )
}

## Propagation with taxPPro
propagated <- bplapply(X = data_ready, BPPARAM = MulticoreParam(workers = n_threads), FUN = function(x) {

    input_tbl <- x |>
        select(NCBI_ID, Attribute, Score, Evidence) |>
        distinct() |>
        complete(NCBI_ID, Attribute, fill = list(Score = 0, Evidence = ''))
    l <- split(input_tbl, factor(input_tbl$NCBI_ID))
    tree$Do(function(node) {
        if (!is.null(l[[node$name]])) {
            node[['table']] <- l[[node$name]]
        }
    })
    tree$Do(asr, traversal = 'post-order')
    tree$Do(inh, traversal = 'pre-order')
    
    ## Get NCBI_IDs with propagation results
    data_tree_tbl <- tree$Get(function(node) node[['table']], simplify = FALSE) |>
        purrr::discard(~ all(is.na(.x))) |>
        dplyr::bind_rows() |>
        dplyr::relocate(NCBI_ID) |>
        dplyr::filter(Evidence %in% c('', 'asr', 'inh') | is.na(Evidence))
    
    ## Combine data from propagation and original annotations
    ## >>Taxon name, etc are missing here<<
    data_with_values <- bind_rows(data_tree_tbl, x)
    
    ## Add missing values for the tree (this is maybe just necessary for
    ## displaying stats
    all_node_names <- tree$Get(function(node) node$name, simplify = TRUE)
    missing_node_names <- all_node_names[which(!all_node_names %in% unique(data_with_values$NCBI_ID))]
    if (length(missing_node_names > 0)) {
        attrs <- unique(x$Attribute)
        empty_df <- data.frame(
            NCBI_ID = sort(rep(missing_node_names, length(attrs))),
            Attribute = rep(attrs, length(missing_node_names)),
            Score = 0,
            Evidence = NA
        )
        final_table <- bind_rows(data_with_values, empty_df)
    } else {
        final_table <- data_with_values
    }
    
    attr_grp <- unique(x$Attribute_group)
    attr_type <- unique(x$Attribute_type)
    
    final_table <- final_table |>
        filter(NCBI_ID != 'ArcBac') |>
        mutate(
            Rank = case_when(
                grepl('^d__', NCBI_ID) ~ 'domain',
                grepl('^p__', NCBI_ID) ~ 'phylum',
                grepl('^c__', NCBI_ID) ~ 'class',
                grepl('^o__', NCBI_ID) ~ 'order',
                grepl('^f__', NCBI_ID) ~ 'family',
                grepl('^g__', NCBI_ID) ~ 'genus',
                grepl('^s__', NCBI_ID) ~ 'species',
                grepl('^t__', NCBI_ID) ~ 'strain'
            )
        ) |>
        mutate(NCBI_ID = sub('^[dpcofgst]__', '', NCBI_ID)) |>
        mutate(Taxon_name = ifelse(is.na(Taxon_name), checkTaxonName(NCBI_ID), Taxon_name)) |>
        filter(!is.na(NCBI_ID) & !is.na(Taxon_name)) |>
        mutate(Frequency = taxPPro:::scores2Freq(Score)) |>
        mutate(Attribute_value = TRUE) |>
        mutate(Attribute_group = attr_grp) |>
        mutate(Attribute_type = attr_type)
    
    tree$Do(function(node) {
        node[['table']] <- NULL
    })
    
    return(final_table)
})

propagated$aerophilicity_prediction_gn |> 
    filter(NCBI_ID %in% sub('g__', '', unique(holdout_gn$NCBI_ID))) |> 
    count(Evidence)

for (i in seq_along(propagated)) {
    write.table(
        x = propagated[[i]],
        file = paste0(names(propagated)[i], '.csv'),
        quote = TRUE, 
        sep = ',', 
        row.names = FALSE
    )
}

log_close()
