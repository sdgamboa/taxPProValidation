library(purrr)
library(dplyr)
library(readr)

fnames <- list.files(path = 'numeric', pattern = '(all|genus|species|strain).*tsv', full.names = TRUE)

dat <- map(fnames, ~ read_tsv(.x, show_col_types = FALSE) ) |> 
    bind_rows() |> 
    relocate(physiology, rank)

write.table(
    x = dat, file = 'numeric/numeric_castor-ltp.tsv', 
    sep = '\t', row.names = FALSE, quote = FALSE
)
