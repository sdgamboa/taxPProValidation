library(vroom)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(rlang)
source('function.R')


phys_name <- 'aerophilicity'

holdouts <- importHoldouts(phys_name)
predictions <- importPredictions(phys_name)

df1 <- holdouts[[1]]
df2 <- predictions[[1]]


getTestSet <- function(df) {
    df |> 
        dplyr::select(.data$NCBI_ID, .data$Attribute, .data$Score) |> 
        dplyr::arrange(.data$NCBI_ID) |> 
        dplyr::distinct() |> 
        tidyr::pivot_wider(
            names_from = 'Attribute', values_from = 'Score', values_fill = 0
        ) |> 
        pivot_longer(
            cols = 2:last_col(), names_to = 'Attribute', values_to = 'Score'
        ) |> 
        rename(tScore = Score) |> 
        mutate(tPosNeg = ifelse(tScore > 0, 1, 0))
}




test_set <- df1 |> 
    select(NCBI_ID, Attribute, Score) |> 
    arrange(NCBI_ID) |> 
    unique() |> 
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score',
        values_fill = 0
    ) |> 
    # mutate(aerotolerant = 0) |>
    pivot_longer(
        cols = 2:last_col(), names_to = 'Attribute', values_to = 'Score'
    ) |> 
    rename(tScore = Score) |> 
    mutate(tPosNeg = ifelse(tScore > 0, 1, 0))

predicted_set <- df2 |> 
    select(NCBI_ID, Attribute, Score) |> 
    arrange(NCBI_ID) |> 
    unique() |> 
    rename(pScore = Score) |> 
    mutate(NCBI_ID = as.character(NCBI_ID)) |> 
    mutate(Attribute = sub('^.*:', '', Attribute)) |> 
    mutate(Attribute = sub('_', ' ', Attribute))

res <- left_join(test_set, predicted_set, by = c('NCBI_ID', 'Attribute')) |> 
    mutate(pScore = ifelse(is.na(pScore), 0, pScore))

splitted_res <- split(res, res$Attribute)

roc_objs <- map(splitted_res, ~  {
    tryCatch(
        error = function(e) e,
        roc(.x$tPosNeg, as.double(.x$pScore))
    )
}) |> 
    discard(is_error)

roc_plot_all <- ggroc(roc_objs, legacy.axes = TRUE) +
    labs(x = '1 - Specificity', y = 'Specificity') +
    scale_color_discrete(name = 'Attribute') + 
    theme_bw() 
