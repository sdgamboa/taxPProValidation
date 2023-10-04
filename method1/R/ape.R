
## Meaning of confidence intervals in ape::ace (or in ASR in general?)
## CI a logical specifying whether to return the 95% confidence intervals of
## the ancestral state estimates (for continuous characters) or the likelihood
## of the different states (for discrete ones'

library(phytools)
data(anoletree)
x<-getStates(anoletree,"tips")
tree<-anoletree
rm(anoletree)
tree
plotTree(tree,type="fan",fsize=0.8,ftype="i")
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)


fitER<-ace(x,tree,model="ER",type="discrete", CI = TRUE)
fitER
res <- round(fitER$lik.anc, 2)
