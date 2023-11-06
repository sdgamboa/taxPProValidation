library(purrr)
library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
physName <- args[[1]]

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

fileNames <- listFiles(physName)
fileNames <- grep('propagated', fileNames, value = TRUE)
sets <- map(fileNames, ~ read_csv(.x, show_col_types = FALSE))

attribute_type <- sets[[1]]$Attribute_type |> 
    {\(y) y[!is.na(y)]}() |> 
    unique()

names(sets) <- sub('^.*/(.*)\\.csv', '\\1', fileNames)
countSummary <- function(dat, thr) {
    evCodes <- c('exp', 'tas', 'nas', 'igc')
    dat |> 
        dplyr::mutate(
            Attribute = sub('--(TRUE|FALSE)', '', Attribute) ## This won't affect the attributes of type multistate-intersection
        ) |> 
        dplyr::mutate(
            Evidence = ifelse(Evidence %in% evCodes, 'Source', Evidence)
        ) |> 
        dplyr::filter(
            (Evidence %in% c('asr', 'tax', 'inh', 'inh2') & Score >= thr) |
                Evidence == 'Source'
        ) |> 
        dplyr::count(Attribute, Evidence) |> 
        tidyr::drop_na() |> 
        tidyr::pivot_wider(
            names_from = 'Evidence', values_from = 'n', values_fill = 0
        ) 
}

if (attribute_type == 'multistate-intersection') {
    thr <- map(sets, ~ pull(.x, Attribute)) |> 
        unlist() |> 
        {\(y) y[!is.na(y)]}() |> 
        unique() |> 
        {\(y) 1 / length(y)}()
    count_summary <- map(sets, countSummary, thr = thr) |> 
        {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}() |> 
        {\(y) set_names(y, sub('^(.*)_(all|genus|species|strain)_(Fold.*)$', '\\1 \\2 \\3', names(y)))}() |> 
        bind_rows(.id = 'Fold') |> 
        separate(
            col = 'Fold', into = c('Physiology', 'Rank', 'Fold'), sep = ' '
        ) |>  
        group_by(Physiology, Rank, Attribute) |> 
        summarise(
            meanSource = mean(Source),
            sdSource = sd(Source),
            meanTaxpool = mean(tax),
            sdTaxpool = sd(tax),
            meanInh1 = mean(inh),
            sdInh1 = sd(inh),
            meanAsr = mean(asr),
            sdAsr = sd(asr),
            meanInh2 = mean(inh2),
            sdInh2 = sd(inh2)
        ) |> 
        ungroup()
} else if (attribute_type == 'binary') {
    count_summary <- map(sets, countSummary, thr = 0.5) |> 
        {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}() |> 
        {\(y) set_names(y, sub('^(.*)_(all|genus|species|strain)_(Fold.*)$', '\\1 \\2 \\3', names(y)))}() |> 
        bind_rows(.id = 'Fold') |> 
        separate(
            col = 'Fold', into = c('Physiology', 'Rank', 'Fold'), sep = ' '
        ) |> 
        group_by(Physiology, Rank, Attribute) |> 
        summarise(
            meanSource = mean(Source),
            sdSource = sd(Source),
            meanTaxpool = mean(tax),
            sdTaxpool = sd(tax),
            meanInh1 = mean(inh),
            sdInh1 = sd(inh),
            meanAsr = mean(asr),
            sdAsr = sd(asr),
            meanInh2 = mean(inh2),
            sdInh2 = sd(inh2)
        ) |> 
        ungroup()
} else if (attribute_type == 'multistate-union') {
    count_summary <- map(sets, countSummary, thr = 0.5) |>
        {\(y) set_names(y, sub('_(test|propagated)', '', names(y)))}() |> 
        {\(y) set_names(y, sub('^(.*)_(all|genus|species|strain)_(.*)_(Fold.*)$', '\\1 \\3 \\2 \\4', names(y)))}() |> 
        bind_rows(.id = 'Fold') |> 
        separate(
            col = 'Fold', into = c('Physiology', 'Attribute', 'Rank', 'Fold'), sep = ' '
        ) |> 
        group_by(Physiology, Rank, Attribute) |> 
        summarise(
            meanSource = mean(Source),
            sdSource = sd(Source),
            meanTaxpool = mean(tax),
            sdTaxpool = sd(tax),
            meanInh1 = mean(inh),
            sdInh1 = sd(inh),
            meanAsr = mean(asr),
            sdAsr = sd(asr),
            meanInh2 = mean(inh2),
            sdInh2 = sd(inh2)
        ) |> 
        ungroup()
}
outputFileName <- paste0(physName, '_count_summary.tsv')
write_tsv(x = count_summary, file = outputFileName)
