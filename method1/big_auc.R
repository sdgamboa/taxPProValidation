library(purrr)
library(readr)
library(dplyr)
library(forcats)
library(tidyr)
library(ggplot2)

fnames <- list.files(pattern = 'auc_table.tsv')
auc_table <- map(fnames, ~ read_tsv(.x, show_col_types = FALSE)) |> 
    bind_rows() |> 
    discard(~ all(is.na(.x)))

df <- auc_table |> 
    rename(physiology_rank = rank, auc = attribute_group_auc) |> 
    mutate(physiology = sub('^(.*)_(all|gn|sp|st)$', '\\1', physiology_rank)) |> 
    mutate(rank = sub('^(.*)_(all|gn|sp|st)$', '\\2', physiology_rank)) |> 
    distinct() |> 
    group_by(physiology, rank) |> 
    mutate(avg = mean(auc)) |> 
    ungroup() |> 
    arrange(-avg) |> 
    mutate() |> 
    mutate(physiology_rank = fct_inorder(physiology_rank))

(
p <- df |> 
    select(physiology, rank, auc) |> 
    distinct() |> 
    ggplot(aes(physiology, auc)) +
    geom_point(shape = 4, color = 'blue') +
    geom_hline(yintercept = 0.5, linetype = 2, color = 'red', linewidth = 0.1) +
    labs(
        x = 'Attribute group', y = 'Attribute group AUC (avg of individual attributes)'
    ) +
    facet_wrap(~ rank, nrow = 1) +
    theme_bw() +
    coord_flip()
)

ggsave(
    filename = 'auc_curves.png', plot = p,
    width = 9, height = 6, units = 'in', bg = 'white'
)
