library(purrr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

dirpath <- file.path('method2', 'validation_tables')

mccFileNames <- list.files(
    path = dirpath, 
    pattern = "mcc", 
    full.names = TRUE
)

countFileNames <- list.files(
    path = dirpath,
    pattern = 'count',
    full.names = TRUE
)

mcc <- map(mccFileNames, ~ {
    read.table(.x, header = TRUE, sep = '\t')
})
names(mcc) <- sub('^.*/(.*)_mcc_.*$', '\\1', mccFileNames)

counts <- map(countFileNames, ~ {
    read.table(.x, header = TRUE, sep = '\t')
})
names(counts) <- sub('^.*/(.*)_count_.*$', '\\1', countFileNames)
counts <- counts[names(mcc)]

l <- map2(mcc, counts, ~ {
    left_join(.x, .y)
})

dat <- l |> 
    bind_rows() |> 
    arrange(Physiology, Rank, Attribute, -meanMCC) |> 
    relocate(Physiology, Rank, Attribute, meanMCC, sdMCC) |> 
    arrange(-meanMCC) |> 
    drop_na()
write_tsv(x = dat, file = 'method2/validation_tables/merged_data.tsv')



dat |> 
    ggplot(aes(x = meanMCC)) +
    geom_histogram()


dat |> 
    arrange(Physiology, Rank, Attribute, -meanMCC) |> 
    relocate(Physiology, Rank, Attribute, meanMCC) |> 
    View()

dat |> 
    mutate(
        Physiology = case_when(
            Physiology == 'COGEM_pathogenicity_rating' ~ 'COGEN_rating',
            Physiology == 'mutation_rate_per_site_per_generation' ~ 'mut_per_year',
            Physiology == 'mutation_rate_per_site_per_year' ~ 'mut_per_generation',
            TRUE ~ Physiology
            
        )
    ) |> 
    ggplot(aes(meanSource, meanMCC)) +
    geom_point(shape = 1) +
    geom_hline(yintercept = 0.4, color = 'red') +
    geom_text_repel(aes(label = Attribute), size = 2) +
    facet_grid(Rank ~ Physiology) +
    labs(
        x = 'Number of annotations from sources (Mean of 10 Folds)',
        y = 'Matthews Correlation Coefficient (Mean of 10 Folds)'
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

ggsave(
    filename = 'mcc_vs_source.png', height = 10, width = 25, dpi = 300
)

