library(phangorn)

# https://www.mail-archive.com/r-sig-phylo@r-project.org/msg05854.htmhttps://www.mail-archive.com/r-sig-phylo@r-project.org/msg05854.htmll


# define ambiguous characters
contrast <- matrix(data = c(
    1,0,0,0,0,0,
    0,1,0,0,0,0,
    0,0,1,0,0,0,
    0,0,0,1,0,0,
    0,0,0,0,1,0,
    1,1,0,0,0,0,
    1,0,1,0,0,0,
    1,0,0,1,0,0,
    0,1,1,0,0,0,
    0,1,0,0,1,0,
    0,0,1,1,0,0,
    0,0,1,0,1,0,
    0,0,0,0,0,1),
    ncol = 6, byrow = TRUE)
dimnames(contrast) <- list(c("0", 
                             "1","2","3","4","01","02","03","12","14","23","24","-"),
                           c("0", "1", "2", "3", "4", "-"))
contrast

traits <- matrix(
    data = c("01", "24", "4", "1", "-", "03", "0", "01", "23", "01"),
    ncol = 1,
    dimnames = list(
        paste0("t",1:10),
        "trait"
    )
)

data <- phyDat(traits,
               type="USER",
               levels=c("0", 
                        "1","2","3","4","01","02","03","12","14","23","24","-"),
               contrast=contrast)
data

tree <- rtree(10)

# maximum parsimony ancestral reconstruction
anc.pars <- ancestral.pars(tree, data, type="MPR")

cols <- c("#33a02c", "#a6cee3", "#1f78b4", "#000000", "#000000", "#000000")
plotAnc(tree, anc.pars, 1, col=cols, cex.pie=.75)
