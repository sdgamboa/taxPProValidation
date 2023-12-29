library(castor)
library(taxPPro)
library(bugphyzz)
library(dplyr)
library(tidyr)

ltp <- ltp()
tip_data <- ltp$tip_data
tree <- ltp$tree
a <- physiologies('aerophilicity')[[1]]
a <- getDataReady(filterData(a))

tip_data_annotated <- left_join(
    tip_data, a, by = 'NCBI_ID', relationship = 'many-to-many'
)

states <- tip_data_annotated |> 
    filter(!is.na(Attribute)) |> 
    select(tip_label, Attribute, Score) |> 
    pivot_wider(names_from = 'Attribute', values_from = 'Score') |> 
    tibble::column_to_rownames(var = 'tip_label') |> 
    as.matrix()


no_annotated_tips <- tip_data$tip_label[!tip_data$tip_label %in% rownames(states)]

unknown_states <- matrix(
    data = rep(0.25, length(no_annotated_tips) * ncol(states)), 
    nrow = length(no_annotated_tips)
)
rownames(unknown_states) <- no_annotated_tips
colnames(unknown_states) <- colnames(states)

tip_priors <- rbind(states, unknown_states)
tip_priors <- tip_priors[tree$tip.label,]

all(rownames(tip_priors) == tree$tip.label)

result <- hsp_mk_model(
    tree = tree, tip_priors = tip_priors, tip_states = NULL,
    check_input = TRUE
)


