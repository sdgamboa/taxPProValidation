library(purrr)
library(dplyr)
library(readr)

fnames <- list.files(path = '.', pattern = '(all|genus|species|strain).*tsv', full.names = TRUE)

dat <- map(fnames, ~ read_tsv(.x, show_col_types = FALSE) ) |> 
    bind_rows() |> 
    relocate(physiology, rank)

write.table(
    x = dat, file = 'numeric_castor-ltp_summary.tsv', 
    sep = '\t', row.names = FALSE, quote = FALSE
)

