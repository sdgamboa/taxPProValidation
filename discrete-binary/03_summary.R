library(purrr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

dir_path <- file.path('discrete-binary')
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
    # mutate(
    #     ltp_per = round(inltp / totalltp * 100),
    #     dat_per = round(indat / totaldat * 100)
    # ) |> 
    select(-indat, -notinltp, -indat, -notindat) |> 
    rename(
        mean_mcc = mean,
        sd_mcc = sd,
        ltp_and_bugphyzz = inltp,
        ltp = totalltp,
        bugphyzz = totaldat
    ) |> 
    mutate(mean_mcc = round(mean_mcc, 2)) |> 
    mutate(sd_mcc = round(sd_mcc, 2)) 

# summary |> 
#     mutate(
#         ltp_per = round(inltp / totalltp * 100),
#         dat_per = round(inltp / totaldat * 100)
#     ) |> 
#     ggplot(aes(ltp_per, mean)) +
#     geom_point() +
#     geom_text_repel(aes(label = attribute)) +
#     labs(
#         y = 'mean MCC', x = 'LTP tips annotated (%)'
#     ) +
#     scale_color_viridis_c(
#         option = 'C', name = 'bugphyzz annotations used (%)'
#     ) +
#     # facet_wrap(~ rank, ncol = 2, nrow = 2) +
#     facet_grid(rank ~ physiology) +
#     theme_bw() +
#     theme(legend.position = 'bottom') 

write.table(
    x = summary, file = 'discrete_multistate_summary.tsv', 
    row.names = FALSE, sep = '\t', quote = FALSE
)



