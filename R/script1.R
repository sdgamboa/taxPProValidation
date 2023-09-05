library(phytools)

data("primate.tree")
data("primate.data")
tree <- primate.tree
data <- primate.data
rm(primate.tree)
rm(primate.data)
activity <- data$Activity_pattern
names(activity) <- rownames(data)
labs <- sample(tree$tip.label, 10)

m <- to.matrix(activity, levels(activity))
m[labs,] <- rep(1/ncol(m), ncol(m))
m['Lagothrix_lagotricha',] <- c(1, 0, 0)
m2 <- matrix(c(1, 0, 0), nrow = 1)
rownames(m2) <- '91'
colnames(m2) <- colnames(m)
mat <- rbind(m, m2)

# m <- m[-which(rownames(m) == 'Lagothrix_lagotricha'),]
fit <- fitMk(tree, m)
ace <- ancr(fit, tips = TRUE, internal = TRUE)

res <- ace$ace
plot(ace, args.plotTree = list(direction = "upwards"))
tips <- sapply(labs, function(x, y) which(y == x), y = tree$tip.label)
add.arrow(tree, tips, arrl = 3, offset = 2, lwd = 2, col = palette()[4])
