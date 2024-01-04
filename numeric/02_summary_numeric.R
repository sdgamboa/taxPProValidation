library(purrr)
library(dplyr)
library(readr)

fnames <- list.files(path = '.', pattern = '(all|genus|species|strain).*tsv', full.names = TRUE)

dat <- map(fnames, ~ read_tsv(.x, show_col_types = FALSE) ) |> 
    bind_rows() |> 
    relocate(physiology, rank) |>
    mutate(r_squared_mean = round(r_squared_mean, 2)) |>
    mutate(r_squared_sd = round(r_squared_sd, 2)) |>
    rename(r2_mean = r_squared_mean, r2_sd = r_squared_sd)

write.table(
    x = dat, file = 'numeric_castor-ltp_summary.tsv', 
    sep = '\t', row.names = FALSE, quote = FALSE
)

