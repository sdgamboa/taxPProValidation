library(logr)
library(vroom)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(rlang)
library(ggpubr)
source('function.R')

args <- commandArgs(trailingOnly = TRUE)
args <- 'gram stain'
phys_name <- gsub(" ", "_", args[[1]])

logfile <- gsub(' ', '_', paste0(phys_name, "_auc_log_file"))
lf <- log_open(logfile, logdir = FALSE, compact = TRUE, show_notes = FALSE)
msg <- paste0(
    'Doing AUC-ROC for ', phys_name
)
log_print(msg, blank_after = TRUE)
holdouts <- importHoldouts(phys_name)
holdouts <- map(holdouts, ~ modify_at(.x, .at = 'Attribute', function(x) gsub(' ', '_', x)))
predictions <- importPredictions(phys_name)
predictions <- map(predictions, ~ modify_at(.x, .at = 'Attribute', function(x) gsub(' ', '_', x)))
roc_res <- map2(holdouts, predictions, ~ doRoc(.x, .y)) |> 
    discard(~ !length(.x))
names(roc_res) <- sub("_(holdout|prediction)", "", names(roc_res))
auc_table <- map(roc_res, ~ getAucTable(.x, phys_name)) |> 
    bind_rows(.id = 'rank')

table_fname <- paste0(phys_name, '_auc_table.tsv')
write.table(
    x = auc_table, file = table_fname, sep = "\t", quote = TRUE, 
    row.names = FALSE
)

roc_plots <- map(roc_res, ~ plotAucRoc(.x))
p <- ggarrange(
    plotlist = roc_plots,
    labels = names(roc_plots)
)
plot_fname <- paste0(phys_name, '_auc_plot.png')
png(
    filename = plot_fname, bg = 'white', 
    width = 11.3, height = 7.5, units = 'in', res = 100
)
p
dev.off()

log_close()

