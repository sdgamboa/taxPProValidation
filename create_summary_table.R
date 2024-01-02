library(dplyr)

fnames <- c(
    list.files(
        path = 'discrete-multistate-intersection',
        pattern = 'summary.*tsv',
        full.names = TRUE
    ),
    list.files(
        path = 'numeric',
        pattern = 'summary.*tsv',
        full.names = TRUE
    )
)

dat <- lapply(fnames, function(x) {
    read.table(x, sep = '\t', header = TRUE, row.names = NULL)
}) |> 
    bind_rows() |> 
    relocate(
        method, rank, physiology, attribute,
        mcc_mean = mean_mcc, mcc_sd = sd_mcc,
        r_squared_mean, r_squared_sd,
    ) |> 
    mutate(
        attribute = ifelse(is.na(attribute), physiology, attribute)
    )

write.table(
    x = dat, file = 'validation_summary.tsv', quote = FALSE, sep = '\t'
)
