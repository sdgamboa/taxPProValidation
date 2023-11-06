library(purrr)
library(readr)
library(dplyr)
library(tidyr)

physName <- 'aerophilicity'


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

names(sets) <- sub('^.*/(.*)\\.csv', '\\1', fileNames)
countSummary <- function(dat) {
    evCodes <- c('exp', 'tas', 'nas', 'igc')
    dat |> 
        dplyr::mutate(
            Evidence = ifelse(Evidence %in% evCodes, 'Source', Evidence)
        ) |> 
        dplyr::filter(
            (Evidence %in% c('asr', 'tax', 'inh', 'inh2') & Score >= 0.25) |
                Evidence == 'Source'
        ) |> 
        dplyr::count(Attribute, Evidence) |> 
        tidyr::drop_na() |> 
        tidyr::pivot_wider(
            names_from = 'Evidence', values_from = 'n', values_fill = 0
        ) 
}
count_summary <- map(sets, countSummary) |> 
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








