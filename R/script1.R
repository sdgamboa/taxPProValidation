library(bugphyzz)
library(dplyr)
library(taxPPro)
library(bugphyzzExports)

phys <- physiologies('aerophilicity')
df <- phys[[1]]

fname <- system.file(
    'extdata', 'attributes.tsv', package = 'bugphyzz', mustWork = TRUE
)

attributes <- readr::read_tsv(fname, show_col_types = FALSE)
valid_attributes <- attributes |> 
    filter(attribute_group == 'aerophilicity') |> 
    pull(attribute) |> 
    unique()

df <- df |> 
    filter(Attribute_value == TRUE) |> 
    filter(Attribute %in% valid_attributes)

set1 <- filter(df, !is.na(NCBI_ID))
set2 <- filter(df, is.na(NCBI_ID))
set1$Rank <- checkRank(set1$NCBI_ID)
set1$Taxon_name <- checkTaxonName(set1$NCBI_ID)
df2 <- bind_rows(set1, set2)
