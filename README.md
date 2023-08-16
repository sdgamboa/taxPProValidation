# taxPProValidation

A repo for testing some validation strategies for bugphyzz

aerophilicity   X
growth temperature  X
disease association X   Note: Too many attributes. Only use the table of AUC values.
biosafety level X   Note: Some biosafety data (sp) couldn't be used for auc beacause there were too few holdouts (like 5).
gram stain X Note: Very low auc values.
shape   X   Note: AUC was good for some and not so good for others (ranks).
motility    X   Note: good auc values.
arrangement X   Note: Bad auc values. Possibly too many individual attributes. Maybe I'll just use the tables for these ones.
coding genes X Note: Very bad performance for AUC values.
genome size X Note: Very bad performance for AUC values.
extreme environment X   Note: Very bad performance for AUC values.
spore formation X Note: very bad.
animal pathogen
plant pathogenicity
optimal ph
width
COGEM pathogenicity rating
mutation rate per site per year
antimicrobial sensitivity
pathogenicity human
hemolysis
length
biofilm forming
growth medium
mutation rate per site per generation
spore shape
health associated
sphingolipid producing
acetate producing
butyrate producing
lactate producing
hydrogen gas producing



These physiologies throw an error beacuse:
1) too few data.
2) all test are 1 and all predicted are 0.

pathogenicity_human
sphingolipid_producing
acetate_producing
butyrate_producing
lactate_producing





