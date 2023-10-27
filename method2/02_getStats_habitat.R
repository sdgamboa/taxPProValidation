library(purrr)
library(dplyr)
library(readr)
library(tidyr)
library(mltools)
library(ggplot2)
library(data.tree)
library(taxPPro)

listFiles <- function() {
    wd <- getwd()
    if (grepl('method2', wd)) {
        list.files(pattern = 'csv', full.names = TRUE)
    } else {
        list.files('method2', pattern = 'csv', full.names = TRUE)
    }
}

fileNames <- listFiles()

## name always with underscore instead of  space
physName <- 'habitat'

physFileNames <- sort(grep(physName, fileNames, value = TRUE))
tbls <- map(physFileNames, ~ read_csv(.x, show_col_types = FALSE))
names(tbls) <- sub('^.*/(.*)\\.csv', '\\1', physFileNames)
testSets <- tbls[grep('test', names(tbls))]
propSets <- tbls[grep('propagated', names(tbls))]








 
data('tree_list')
ncbi_tree <- as.Node(tree_list)
ncbi_nodes <- ncbi_tree$Get(
    attribute = 'name', filterFun = function(node) grepl('^[gst]__', node$name)
) |> 
    unname()
rm(ncbi_tree)
rm(tree_list)
gc()

attrs <- map(tbls, ~ unique(pull(.x, Attribute))) |> 
    unlist() |> 
    unique()
thr <- 1 / length(attrs)

attr_typ <- unique(propSets[[1]]$Attribute_type)
attr_typ <- attr_typ[!is.na(attr_typ)]
if (attr_typ == 'multisatate-union') {
    thr <- 0.5
}

lgl_vct <- map_lgl(testSets, ~ !is.null(.x))

testSets <- testSets[lgl_vct]
propSets <- propSets[lgl_vct]

testSets <- map(testSets, distinct)
propSets <- map(propSets, distinct)


getNegatives <- function(dat) {
    ranks <- unique(dat$Rank)
    ranks <- factor(ranks, levels = c('genus', 'species', 'strain'))
    ranks <- case_when(
        ranks == 'genus' ~ 'g__',
        ranks == 'species' ~ 's__',
        ranks == 'strain' ~ 't__'
    )
    ranks <- paste0('^(', paste0(ranks, collapse = '|'), ')')
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
            new_negatives <- sample(selected_nodes, n) 
            output <- data.frame(
                NCBI_ID = new_negatives,
                Attribute = unique(dat$Attribute),
                tScore = 0
            )
        } else {
            return(NULL)
        }
        
    } else {
        ncbi_ids <- unique(dat$NCBI_ID)
        n <- length(ncbi_ids)
        selected_nodes <- selected_nodes[!selected_nodes %in% ncbi_ids]
        set.seed(20308)
        new_negatives <- sample(selected_nodes, n) 
        output <- data.frame(
            NCBI_ID = new_negatives,
            Attribute = unique(dat$Attribute),
            Score = 0
        )
    }
    return(output)
}

# getNegatives(df)

testSetsPlus <- map(testSets, ~ bind_rows(.x, getNegatives(.x)))

sets <- map2(
    .x = testSetsPlus,
    .y = propSets,
    .f = ~ {
        x <- select(.x, NCBI_ID, Attribute, tScore = Score)
        y <- select(.y, NCBI_ID, Attribute, pScore = Score)
        left_join(x, y, by = c('NCBI_ID', 'Attribute')) |> 
            mutate_at(
                .vars = c('tScore', 'pScore'), .funs = function(x) {
                    ifelse(is.na(x), 0, x)
            })
    }
) |>
    # map(~ {
    #     complete(
    #         .x, NCBI_ID, Attribute, fill = list(tScore = 0, pScore = 0)
    #     )
    # }) |> 
    map(~ {
        mutate(
            .x,
            Attribute = sub('--(TRUE|FALSE)', '', Attribute),
            tPosNeg = ifelse(tScore > thr, 1, 0),
            pPosNeg = ifelse(pScore > thr, 1, 0),
            PosNeg = case_when(
                tPosNeg == 1 & pPosNeg == 1 ~ 'TP',
                tPosNeg == 1 & pPosNeg == 0 ~ 'FN',
                tPosNeg == 0 & pPosNeg == 0 ~ 'TN',
                tPosNeg == 0 & pPosNeg == 1 ~ 'FP'
                
            )
        )
    }) |> 
    {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}()

pn <- map(sets, ~ {
    dat <- .x |> 
        mutate(
            PosNeg = factor(
                PosNeg, levels = c('TP', 'FP', 'FN', 'TN')
            )
        ) |> 
        group_by(Attribute, PosNeg, .drop = FALSE) |> 
        summarise(
            n = n()
        ) |> 
        ungroup()
    dat
}) |> 
    bind_rows(.id = 'set_names') |> 
    mutate(
        set_names = sub(
            '_(all|genus|species|strain)_(.*)_(Fold[0-9]+)',
            " \\2 \\1 \\3",
            set_names
        )
    ) |> 
    separate(
        col = 'set_names', into = c('Attribute_group', 'Attribute0', 'Rank', 'Fold'),
        sep = " "
    ) |> 
    select(-Attribute0)

rank_order <- c('all', 'genus', 'species', 'strain')

(
    p1 <- pn |> 
        mutate(Attribute = paste0(Attribute_group, '_', Attribute)) |> 
        mutate(
            Rank = factor(Rank, levels = rank_order, ordered = TRUE)
        ) |> 
        ggplot(aes(Attribute, n)) +
        # geom_point(
        #     aes(color = PosNeg),
        #     shape = 21,
        #     position = position_jitterdodge(dodge.width = 0.9),
        #     alpha = 1
        # ) +
        geom_point(
            aes(color = PosNeg),
            shape = 21,
            position = 'jitter',
            # position = position_jitterdodge(dodge.width = 0.9),
            alpha = 1
        ) +
        # geom_boxplot(aes(color = PosNeg), alpha = 0) +
        facet_wrap(~Rank, scales = 'free_x') +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
        
)

mcc_res0 <- map(sets, ~ {
   dats <- split(.x, factor(.x$Attribute))
   map(dats, function(y) {
       mcc(preds = y$pPosNeg, actuals = y$tPosNeg)
   })
}) |> 
    list_flatten()

mcc_res <- data.frame(
    dat_name = gsub(' ', '_', names(mcc_res0)),
    mcc = as.double(mcc_res0)
) |> 
    mutate(
        dat_name = sub(
            '_(all|genus|species|strain)_(.*)_(Fold[0-9]+)_', 
            " \\2 \\1 \\3 ",
            dat_name 
        )
    ) |>  
    separate(
        col = 'dat_name', into = c('Attribute_group', 'Attribute0', 'Rank', 'Fold', 'Attribute'),
        sep = " "
    )


(
    p2 <- mcc_res |> 
        mutate(Attribute = paste0(Attribute_group, '_', Attribute)) |> 
        mutate(
            Rank = factor(Rank, levels = rank_order, ordered = TRUE)
        ) |> 
        ggplot(aes(Attribute, mcc)) +
        geom_boxplot(aes(group = Attribute, fill = Attribute)) +
        facet_wrap(~ Rank, scales = 'free_x') +
        # geom_point() +
        # scale_y_continuous(
        #     breaks = seq(0, 1, 0.1), limits = c(0.3, 1.1)
        # ) + 
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'none'
        )
)


# Export files ------------------------------------------------------------

wd <- getwd()
if (grepl('method2', wd)) {
    physName <- physName
} else {
    physName <- file.path('method2', physName)
}

pn_fname <- paste0(physName, '_', 'pn.tsv')
write_tsv(pn, pn_fname)

mcc_results_fname <- paste0(physName, '_', 'mcc.tsv')
write_tsv(mcc_res, mcc_results_fname)

p1_fname <- paste0(physName, '_', 'pn.png')
ggsave(p1_fname, p1, width = 8, height = 7, units = 'in')

p2_fname <- paste0(physName, '_', 'mcc.png')
ggsave(p2_fname, p2, width = 8, height = 7, units = 'in')

