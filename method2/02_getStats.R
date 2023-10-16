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

physName <- 'aerophilicity'
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

result <- map(sets, ~ {
   dats <- split(.x, factor(.x$Attribute)) 
   map(dats, function(y) {
       mcc(preds = y$pPosNeg, actuals = y$tPosNeg)
   })
}) |> 
    list_flatten() |> 
    as.data.frame() |> 
    t() |> 
    as.data.frame() |> 
    tibble::rownames_to_column(var = 'fold') |> 
    as_tibble() |> 
    set_names(c('fold', 'mcc')) |> 
    mutate(
        fold = sub('_Fold[0-9]+_', " Fold ", fold)
    ) |> 
    separate(
        col = 'fold', into = c('Attribute_group', 'Fold', 'Attribute'),
        sep = " "
    )

result |> 
    mutate(Attribute = paste0(Attribute_group, '_', Attribute)) |> 
    ggplot(aes(Attribute, mcc)) +
    geom_boxplot(aes(group = Attribute)) +
    scale_y_continuous(
        breaks = seq(0, 1, 0.1), limits = c(0.3, 1.1)
    ) + 
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
    


