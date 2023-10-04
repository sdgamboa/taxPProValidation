library(phytools)
data('anoletree')
t <- anoletree
rm(anoletree)
tipStates <- as.factor(getStates(t, 'tips'))
mat <- to.matrix(tipStates, seq = levels(tipStates))

fit <- fitMk(
    tree = t, x = mat, model = "ARD", pi = "fitzjohn", logscale=TRUE,
    lik.func = "pruning", CI = TRUE  
)
ace <- ancr(fit, tips = TRUE, CI = TRUE)
