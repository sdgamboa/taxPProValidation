library(readr)
library(purrr)
library(dplyr)
library(tidyr)
library(mltools)

# args <- commandArgs(trailingOnly = TRUE)
# physName <- gsub('_', ' ', args[[1]])
physName <- 'spore formation'

links_tsv <- system.file(
    'extdata', 'links.tsv', package = 'bugphyzz', mustWork = TRUE
)
links <- read_tsv(links_tsv, show_col_types = FALSE)

attribute_type <- links |> 
    filter(physiology == physName) |> 
    pull(attribute_type)

if (!length(attribute_type)) {
    # quit(save = 'no')
    message('not here')
}

if (attribute_type %in% c('numeric', 'range')) {
    attribute_type <- 'multistate-intersection'
}

listFiles <- function(phys_name = NULL) {
    phys_name <- gsub(' ', '_', phys_name)
    wd <- getwd()
    if (grepl('method2', wd)) {
        physFileNames <- list.files(pattern = 'csv', full.names = TRUE)
    } else {
        physFileNames <- list.files('method2', pattern = 'csv', full.names = TRUE)
    }
    if (!is.null(phys_name))
        physFileNames <- sort(grep(phys_name, physFileNames, value = TRUE))
    return(physFileNames)
}


if (attribute_type == 'multistate-intersection') {
    fileNames <- listFiles(physName)
    tbls <- map(fileNames, ~ read.csv(.x))
    names(tbls) <- sub('^.*/(.*)\\.csv', '\\1', fileNames)
    testSets <- tbls[grep('test', names(tbls))]
    propSets <- tbls[grep('propagated', names(tbls))]
    
    attrs <- map(tbls, ~ unique(pull(.x, Attribute))) |> 
        unlist() |> 
        unique()
    thr <- 1 / length(attrs)
    
    foldData <- map2(.x = testSets, .y = propSets, .f = ~ {
        testDat <- select(.x, NCBI_ID, Attribute, tScore = Score) 
        predDat <- select(.y, NCBI_ID, Attribute, pScore = Score)
        left_join(testDat, predDat, by = c('NCBI_ID', 'Attribute')) |> 
            complete(NCBI_ID, Attribute, fill = list(tScore = 0, pScore = 0)) |> 
            mutate(
                tPosNeg = ifelse(tScore > thr, 1, 0),
                pPosNeg = ifelse(pScore > thr, 1, 0),
                PosNeg = case_when(
                    tPosNeg == 1 & pPosNeg == 1 ~ 'TP',
                    tPosNeg == 1 & pPosNeg == 0 ~ 'FN',
                    tPosNeg == 0 & pPosNeg == 0 ~ 'TN',
                    tPosNeg == 0 & pPosNeg == 1 ~ 'FP'
                ),
                PosNeg = factor(PosNeg, levels = c('TP', 'FP', 'FN', 'TN'))
            ) |> 
            group_by(Attribute, PosNeg, .drop = FALSE) |>
            summarise(n = n()) |>
            ungroup() |> 
            pivot_wider(
                names_from = 'PosNeg', values_from = 'n'
            ) |> 
            rowwise() |> 
            mutate(
                MCC = mltools::mcc(TP = TP, FP = FP, TN = TN, FN = FN)
            )
    }) |> 
        {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}() |> 
        {\(y) set_names(y, sub('^(.*)_(all|genus|species|strain)_(Fold.*)$', '\\1 \\2 \\3', names(y)))}() |> 
        bind_rows(.id = 'Fold') |> 
        separate(
            col = 'Fold', into = c('Physiology', 'Rank', 'Fold'), sep = ' '
        ) 
    
    foldDataSummary1 <- foldData |> 
        group_by(Physiology, Rank, Attribute) |> 
        summarise(
            meanTP = mean(TP),
            sdTP = sd(TP),
            meanFP = mean(FP),
            sdFP = sd(FP),
            meanFN = mean(FN),
            sdFN = sd(FN),
            meanTN = mean(TN),
            sdTN = sd(TN),
            meanMCC = mean(MCC),
            sdMCC = sd(MCC)
        ) |> 
        ungroup()
    
} else if (attribute_type == 'binary') {
    fileNames <- listFiles(physName)
    tbls <- map(fileNames, ~ read.csv(.x))
    names(tbls) <- sub('^.*/(.*)\\.csv', '\\1', fileNames)
    testSets <- tbls[grep('test', names(tbls))]
    propSets <- tbls[grep('propagated', names(tbls))]
    
    thr <- 0.5
    
    testSets <- map(testSets, ~ {
        chr <- .x$Attribute |> 
            {\(y) y[!is.na(y)]}() |> 
            {\(y) sub('--.*$', '', y)}() |> 
            {\(y) gsub(' ', '_', y)}() |> 
            unique() |> 
            paste0(c('--TRUE', '--FALSE'))
        .x |> 
            mutate(Attribute = gsub(' ', '_', Attribute)) |> 
            mutate(Attribute = factor(Attribute, levels = chr)) |> 
            complete(NCBI_ID, Attribute, fill = list(Score = 0)) |> 
            group_by(NCBI_ID) |> 
            slice_max(Score, n = 1, with_ties = FALSE) |> 
            ungroup() |> 
            mutate(
                Score = ifelse(grepl('--FALSE', Attribute), 0, 1),
                Attribute = sub('--(TRUE|FALSE)$', '', Attribute)
            )
    })
    
    
    propSets <- map(propSets, ~ {
        chr <- .x$Attribute |> 
            {\(y) y[!is.na(y)]}() |> 
            {\(y) sub('--.*$', '', y)}() |> 
            {\(y) gsub(' ', '_', y)}() |> 
            unique() |> 
            paste0(c('--TRUE', '--FALSE'))
        .x |> 
            mutate(Attribute = gsub(' ', '_', Attribute)) |> 
            mutate(Attribute = factor(Attribute, levels = chr)) |> 
            # complete(NCBI_ID, Attribute, fill = list(Score = 0)) |> 
            group_by(NCBI_ID) |> 
            slice_max(Score, n = 1, with_ties = FALSE) |> 
            ungroup() |> 
            mutate(
                Score = ifelse(grepl('--FALSE', Attribute), 0, 1),
                Attribute = sub('--(TRUE|FALSE)$', '', Attribute)
            )
    })
    
    foldData <- map2(.x = testSets, .y = propSets, .f = ~ {
        testDat <- select(.x, NCBI_ID, Attribute, tScore = Score) 
        predDat <- select(.y, NCBI_ID, Attribute, pScore = Score)
        left_join(testDat, predDat, by = c('NCBI_ID', 'Attribute')) |> 
            complete(NCBI_ID, Attribute, fill = list(tScore = 0, pScore = 0)) |> 
            mutate(
                # tPosNeg = ifelse(tScore > thr, 1, 0),
                tPosNeg = tScore,
                # pPosNeg = ifelse(pscore > thr, 1, 0),
                pPosNeg = pScore,
                PosNeg = case_when(
                    tPosNeg == 1 & pPosNeg == 1 ~ 'TP',
                    tPosNeg == 1 & pPosNeg == 0 ~ 'FN',
                    tPosNeg == 0 & pPosNeg == 0 ~ 'TN',
                    tPosNeg == 0 & pPosNeg == 1 ~ 'FP'
                ),
                PosNeg = factor(PosNeg, levels = c('TP', 'FP', 'FN', 'TN'))
            ) |> 
            group_by(Attribute, PosNeg, .drop = FALSE) |>
            summarise(n = n()) |>
            ungroup() |> 
            pivot_wider(
                names_from = 'PosNeg', values_from = 'n'
            ) |> 
            rowwise() |> 
            mutate(
                MCC = mltools::mcc(TP = TP, FP = FP, TN = TN, FN = FN)
            )
    }) |> 
        {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}() |> 
        {\(y) set_names(y, sub('^(.*)_(all|genus|species|strain)_(Fold.*)$', '\\1 \\2 \\3', names(y)))}() |> 
        bind_rows(.id = 'Fold') |> 
        separate(
            col = 'Fold', into = c('Physiology', 'Rank', 'Fold'), sep = ' '
        ) 
    
    foldDataSummary1 <- foldData |> 
        group_by(Physiology, Rank, Attribute) |> 
        summarise(
            meanTP = mean(TP),
            sdTP = sd(TP),
            meanFP = mean(FP),
            sdFP = sd(FP),
            meanFN = mean(FN),
            sdFN = sd(FN),
            meanTN = mean(TN),
            sdTN = sd(TN),
            meanMCC = mean(MCC),
            sdMCC = sd(MCC)
        ) |> 
        ungroup()
    
} else if (attribute_type == 'multistate-union') {
    fileNames <- listFiles(physName)
    tbls <- map(fileNames, ~ read.csv(.x))
    names(tbls) <- sub('^.*/(.*)\\.csv', '\\1', fileNames)
    testSets <- tbls[grep('test', names(tbls))]
    propSets <- tbls[grep('propagated', names(tbls))]
    
    thr <- 0.5
    
    testSets <- map(testSets, ~ {
        chr <- .x$Attribute |> 
            {\(y) y[!is.na(y)]}() |> 
            {\(y) sub('--.*$', '', y)}() |> 
            {\(y) gsub(' ', '_', y)}() |> 
            unique() |> 
            paste0(c('--TRUE', '--FALSE'))
        .x |> 
            mutate(Attribute = gsub(' ', '_', Attribute)) |> 
            mutate(Attribute = factor(Attribute, levels = chr)) |> 
            complete(NCBI_ID, Attribute, fill = list(Score = 0)) |> 
            group_by(NCBI_ID) |> 
            slice_max(Score, n = 1, with_ties = FALSE) |> 
            ungroup() |> 
            mutate(
                Score = ifelse(grepl('--FALSE', Attribute), 0, 1),
                Attribute = sub('--(TRUE|FALSE)$', '', Attribute)
            )
    })
    
    
    propSets <- map(propSets, ~ {
        chr <- .x$Attribute |> 
            {\(y) y[!is.na(y)]}() |> 
            {\(y) sub('--.*$', '', y)}() |> 
            {\(y) gsub(' ', '_', y)}() |> 
            unique() |> 
            paste0(c('--TRUE', '--FALSE'))
        .x |> 
            mutate(Attribute = gsub(' ', '_', Attribute)) |> 
            mutate(Attribute = factor(Attribute, levels = chr)) |> 
            # complete(NCBI_ID, Attribute, fill = list(Score = 0)) |> 
            group_by(NCBI_ID) |> 
            slice_max(Score, n = 1, with_ties = FALSE) |> 
            ungroup() |> 
            mutate(
                Score = ifelse(grepl('--FALSE', Attribute), 0, 1),
                Attribute = sub('--(TRUE|FALSE)$', '', Attribute)
            )
    })
    
    foldData <- map2(.x = testSets, .y = propSets, .f = ~ {
        testDat <- select(.x, NCBI_ID, Attribute, tScore = Score) 
        predDat <- select(.y, NCBI_ID, Attribute, pScore = Score)
        left_join(testDat, predDat, by = c('NCBI_ID', 'Attribute')) |> 
            complete(NCBI_ID, Attribute, fill = list(tScore = 0, pScore = 0)) |> 
            mutate(
                # tPosNeg = ifelse(tScore > thr, 1, 0),
                tPosNeg = tScore,
                # pPosNeg = ifelse(pscore > thr, 1, 0),
                pPosNeg = pScore,
                PosNeg = case_when(
                    tPosNeg == 1 & pPosNeg == 1 ~ 'TP',
                    tPosNeg == 1 & pPosNeg == 0 ~ 'FN',
                    tPosNeg == 0 & pPosNeg == 0 ~ 'TN',
                    tPosNeg == 0 & pPosNeg == 1 ~ 'FP'
                ),
                PosNeg = factor(PosNeg, levels = c('TP', 'FP', 'FN', 'TN'))
            ) |> 
            group_by(Attribute, PosNeg, .drop = FALSE) |>
            summarise(n = n()) |>
            ungroup() |> 
            pivot_wider(
                names_from = 'PosNeg', values_from = 'n'
            ) |> 
            rowwise() |> 
            mutate(
                MCC = mltools::mcc(TP = TP, FP = FP, TN = TN, FN = FN)
            )
    }) |> 
        {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}() |> 
        {\(y) set_names(y, sub('^(.*)_(all|genus|species|strain)_(.*)_(Fold.*)$', '\\1 \\3 \\2 \\4', names(y)))}() |> 
        bind_rows(.id = 'Fold') |> 
        separate(
            col = 'Fold', into = c('Physiology', 'Attribute', 'Rank', 'Fold'), sep = ' '
        ) 
    
    foldDataSummary1 <- foldData |> 
        group_by(Physiology, Rank, Attribute) |> 
        summarise(
            meanTP = mean(TP),
            sdTP = sd(TP),
            meanFP = mean(FP),
            sdFP = sd(FP),
            meanFN = mean(FN),
            sdFN = sd(FN),
            meanTN = mean(TN),
            sdTN = sd(TN),
            meanMCC = mean(MCC),
            sdMCC = sd(MCC)
        ) |> 
        ungroup()
}

## Export tables ####
outputFileName1 <- gsub(' ', '_', paste0(physName, '_fold_data.tsv'))
write_tsv(x = foldData, file = outputFileName1)
outputFileName2 <- gsub(' ', '_', paste0(physName, '_mcc_summary.tsv'))
write_tsv(x = foldDataSummary1, file = outputFileName2)

