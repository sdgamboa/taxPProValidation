library(purrr)
library(readr)
library(dplyr)

fileNames <- list.files(
    'method2/validation_tables/', pattern = "mcc", full.names = TRUE
)

dat <- map(fileNames, ~ read_tsv(.x, show_col_types = FALSE)) |> 
    bind_rows() |> 
    arrange(-meanMCC)
write_tsv(x = dat, file = 'method2/validation_tables/merged_data.tsv')
