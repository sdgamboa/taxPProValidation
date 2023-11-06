library(purrr)
library(readr)
library(dplyr)
library(tidyr)

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





