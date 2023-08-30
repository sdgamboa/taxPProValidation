library(mltools)

actual <- rep(c(1, 0), times = c(20, 380))
preds <- rep(c(1, 0, 1, 0), times = c(15, 5, 5, 375))
mcc(preds = preds, actuals = actual)


conf_matrix <- matrix(c(15, 5, 5, 375), nrow = 2)
mcc(confusionM = conf_matrix)
