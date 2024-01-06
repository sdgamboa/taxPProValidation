library(purrr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

dir_path <- file.path('.')
fnames <- list.files(dir_path, '(counts|mcc).*tsv', full.names = TRUE, recursive = FALSE)

l <- map(fnames, ~ read_tsv(.x, show_col_types = FALSE))
names(l) <- sub("(^.*/)*(.*)\\.tsv$", "\\2", fnames)


counts <- l[grep('counts', names(l))] |> 
    bind_rows(.id = 'var') |> 
    mutate(var = sub('_counts$', '', var)) |> 
    mutate(physiology = gsub("_", " ", sub("^(.*)_(all|genus|species|strain).*$", "\\1", var))) |> 
    mutate(rank = sub('^.*(all|genus|species|strain).*$', "\\1", var)) |> 
    mutate(method = sub('^.*(all|genus|species|strain)_(.*)$', '\\2', var)) |> 
    select(-var) |> 
    {\(y) set_names(y, tolower(names(y)))}()

mcc <- l[grep('mcc', names(l))] |> 
    bind_rows() |> 
    {\(y) set_names(y, tolower(names(y)))}()

summary <- left_join(
    mcc, counts, by = c('physiology', 'attribute', 'rank', 'method')
) |> 
    rename(
        mean_mcc = mean,
        sd_mcc = sd
    ) |> 
    mutate(mean_mcc = round(mean_mcc, 2)) |> 
    mutate(sd_mcc = round(sd_mcc, 2)) 

write.table(
    x = summary, file = 'discrete_multistate_summary.tsv', 
    row.names = FALSE, sep = '\t', quote = FALSE
)



