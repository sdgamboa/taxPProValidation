importHoldouts <- function(phys = 'aerophilicity') {
    regex <- paste0('^', phys, '_holdout.*csv')
    fnames <- list.files(pattern = regex)
    purrr::map(
        fnames, ~ vroom::vroom(.x, delim = ',', show_col_types = FALSE)
    ) |>
        purrr::set_names(sub('\\.csv$', '', fnames)) |> 
        purrr::map(~ {
            .x$NCBI_ID <- sub('^\\w__', '', .x$NCBI_ID)
            .x
        })
}

importPredictions <- function(phys = 'aerophilicity') {
    regex <- paste0('^', phys, '_prediction.*csv')
    fnames <- list.files(pattern = regex)
    purrr::map(
        fnames, ~ vroom::vroom(
            .x, delim = ',', show_col_types = FALSE, guess_max = 500000,
            col_types = vroom::cols(
                NCBI_ID = vroom::col_character(),
                Attribute = vroom::col_character(),
                Score = vroom::col_double(),
                Evidence = vroom::col_character(),
                Attribute_source = vroom::col_character(),
                Confidence_in_curation= vroom::col_character()
            )
        )
    ) |> 
        purrr::set_names(sub('\\.csv$', '', fnames)) |> 
        purrr::map(~ {
            .x$NCBI_ID <- sub('^\\w__', '', .x$NCBI_ID)
            .x
        })
}

getTestSet <- function(df) {
    df |> 
        dplyr::select(NCBI_ID, Attribute, Score) |> 
        dplyr::arrange(NCBI_ID) |> 
        dplyr::distinct() |> 
        tidyr::pivot_wider(
            names_from = 'Attribute', values_from = 'Score', values_fill = 0
        ) |> 
        tidyr::pivot_longer(
            cols = 2:last_col(), names_to = 'Attribute', values_to = 'Score'
        ) |> 
        dplyr::rename(tScore = Score) |> 
        dplyr::mutate(tPosNeg = ifelse(tScore > 0, 1, 0))
}

getPredictedSet <- function(df) {
    df |> 
        dplyr::select(NCBI_ID, Attribute, Score) |> 
        dplyr::arrange(NCBI_ID) |> 
        dplyr::distinct() |> 
        dplyr::rename(pScore = Score) |> 
        dplyr::mutate(NCBI_ID = as.character(NCBI_ID)) |> 
        dplyr::mutate(Attribute = sub('^.*:', '', Attribute)) |> 
        dplyr::mutate(Attribute = sub('_', ' ', Attribute))
}

doRoc <- function(df1, df2) {
    testSet <- getTestSet(df1)
    predictedSet <- getPredictedSet(df2)
    res <- dplyr::left_join(
        testSet, predictedSet, by = c('NCBI_ID', 'Attribute')
    ) |> 
        dplyr::mutate(pScore = ifelse(is.na(pScore), 0, pScore))
    splitted_res <- split(res, res$Attribute)
    roc_objs <- purrr::map(splitted_res, ~  {
        tryCatch(
            error = function(e) e,
            pROC::roc(.x$tPosNeg, as.double(.x$pScore))
        )
    }) |> 
        purrr::discard(is_error)
    return(roc_objs)
}

getAucTable <- function(rocRes, physName) {
    aucs <- purrr:::map_dbl(rocRes, ~ pROC::auc(.x))
    df <- data.frame(attribute = names(aucs), attribute_auc = aucs)
    rownames(df) <- NULL
    df <- df |> 
        dplyr::mutate(attribute_group = physName) |> 
        dplyr::relocate(attribute_group) |> 
        dplyr::mutate(attribute_group_auc = mean(attribute_auc)) |> 
        dplyr::mutate(
            attribute_auc = round(attribute_auc, digits = 2),
            attribute_group_auc = round(attribute_group_auc, digits = 2)
        )
    return(df)
}

plotAucRoc <- function(rocRes) {
    pROC::ggroc(rocRes, legacy.axes = TRUE) +
        ggplot2::labs(x = '1 - Specificity', y = 'Specificity') +
        ggplot2::scale_color_discrete(name = 'Attribute') + 
        ggplot2::theme_bw() 
}
