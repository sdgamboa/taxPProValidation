importHoldouts <- function(phys = 'aerophilicity') {
    regex <- paste0(phys, '.*holdout.*csv')
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
    regex <- paste0(phys, '.*prediction.*csv')
    fnames <- list.files(pattern = regex)
    purrr::map(
        fnames, ~ vroom::vroom(
            .x, delim = ',', show_col_types = FALSE, guess_max = 500000,
            col_types = vroom::cols(
                NCBI_ID = vroom::col_character(),
                Attribute = vroom::col_character(),
                Score = vroom::col_double(),
                Evidence = vroom::col_character(),
            )
        )
    ) |> 
        purrr::set_names(sub('\\.csv$', '', fnames)) |> 
        purrr::map(~ {
            .x$NCBI_ID <- sub('^\\w__', '', .x$NCBI_ID)
            .x
        })
}




