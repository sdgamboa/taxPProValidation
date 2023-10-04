library(phytools)
data("primate.tree")
data("primate.data")
tree <- primate.tree
data <- primate.data
rm(primate.tree)
rm(primate.data)
activity <- data$Activity_pattern
names(activity) <- rownames(data)
# n <- round(nrow(m) * 0.9)
# labs <- sample(tree$tip.label, n)
# m[labs,] <- rep(1/ncol(m), ncol(m))
m <- to.matrix(activity, levels(activity))
m[,] <- rep(1/ncol(m), ncol(m))

labs <- rownames(m)[which(grepl('Galago', rownames(m)))]


for (i in seq_along(labs)) {
    m[labs[i],] <- c(0, 0, 1)
}

# m[labs,] <- rep(1, ncol(m))
fit <- fitMk(
    tree = tree, x = m, model = 'ARD', pi = "fitzjohn", logscale = TRUE, 
    lik.func = "pruning"
)
ace <- ancr(fit, tips = TRUE)
res <- ace$ace
plot(ace, args.plotTree = list(direction = "upwards"))
tips <- sapply(labs, function(x, y) which(y == x), y = tree$tip.label)
add.arrow(tree, tips, arrl = 3, offset = 2, lwd = 2, col = palette()[4])
