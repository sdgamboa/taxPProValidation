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
    ),
    list.files(
        path = 'discrete-binary',
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
        r2_mean, r2_sd,
    ) |> 
    mutate(
        attribute = ifelse(is.na(attribute), physiology, attribute)
    ) |> 
    mutate(
        ltp_bp_per = round(ltp_bp / ltp * 100),
        ltp_bp_phys_per = round(ltp_bp_phys / ltp * 100)
    )


write.table(
    x = dat, file = 'validation_summary.tsv', quote = FALSE, sep = '\t',
    row.names = FALSE
)
