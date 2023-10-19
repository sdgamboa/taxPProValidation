library(purrr)
library(dplyr)
library(readr)
library(tidyr)
library(mltools)
library(ggplot2)

listFiles <- function() {
    wd <- getwd()
    if (grepl('method2', wd)) {
        list.files(pattern = 'csv', full.names = TRUE)
    } else {
        list.files('method2', pattern = 'csv', full.names = TRUE)
    }
}

(fileNames <- listFiles())

physName <- 'optimal_ph'
# rankVar <- 'all'

physFileNames <- sort(grep(physName, fileNames, value = TRUE))
tbls <- map(physFileNames, ~ read_csv(.x, show_col_types = FALSE))
names(tbls) <- sub('^.*/(.*)\\.csv', '\\1', physFileNames)
testSets <- tbls[grep('test', names(tbls))]
propSets <- tbls[grep('propagated', names(tbls))]


attrs <- map(tbls, ~ unique(pull(.x, Attribute))) |> 
    unlist() |> 
    unique()

thr <- 1 / length(attrs)

sets <- map2(
    .x = testSets,
    .y = propSets,
    .f = ~ {
        x <- select(.x, NCBI_ID, Attribute, tScore = Score)
        y <- select(.y, NCBI_ID, Attribute, pScore = Score)
        left_join(x, y, by = c('NCBI_ID', 'Attribute'))
        # full_join(x, y, by = c('NCBI_ID', 'Attribute'))
    }
) |> 
    map(~ {
        complete(
            .x, NCBI_ID, Attribute, fill = list(tScore = 0, pScore = 0)
        )
    }) |> 
    map(~ {
        mutate(
            .x,
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
    # dat <- count(.x, Attribute, PosNeg)
    # colnames(dat)[1] <- c('PosNeg')
    dat
    # m <- matrix(as.character(table(.x$PosNeg)), nrow = 2, byrow = TRUE)
    # colnames(m) <- c('T', 'F')
    # rownames(m) <- c('T', 'F')
    # m
}) |> 
    bind_rows(.id = 'set_names') |> 
    mutate(
        set_names = sub('_(all|genus|species|strain)_(Fold[0-9]+)', " \\1 \\2", set_names)
    ) |> 
    separate(
        col = 'set_names', into = c('Attribute_group', 'Rank', 'Fold'),
        sep = " "
    )

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
        facet_wrap(~Rank) +
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
            '_(all|genus|species|strain)_(Fold[0-9]+)_', " \\1 \\2 ", dat_name 
        )
    ) |>  
    separate(
        col = 'dat_name', into = c('Attribute_group', 'Rank', 'Fold', 'Attribute'),
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
        facet_wrap(~Rank) +
        geom_point(position = 'jitter', shape = 2) +
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

