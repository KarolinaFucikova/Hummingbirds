# this script explores character evolution in the original hummingbird data set
# i.e., without the added taxa for which we have no molecular data.

library("corHMM")
# the package above is for estimating rate heterogeneity in trait evolution
# https://doi.org/10.1111/2041-210X.13534

library("phytools")
library("parallel")
library("geiger")
library("ape")

# the tree first needs to be cleaned as for the Pagel's lambda calculations

HumData <- read.csv(file="data/HumData2.csv", header=TRUE)

humtree <- "data/HumtreeCons.tre"
Hum<-read.nexus(humtree)
#Hum<-reroot(Hum,103)
plot(Hum)

#check if names match
name.check(Hum, HumData, data.names = HumData$Species)

#fixing common mistake - space instead of underscore
#replace any whitespace with underscore
#dataset$Species<-gsub(" ", "_", HumData$Species)

#if there's nothing to fix but rather taxa should be dropped
#to exclude taxa that are extra in either data set
#first find species in tree that ARE in the dataset
matches <- match(HumData$Species, Hum$tip.label, nomatch = 0)

#Remove species in dataset not in the tree
data <-subset(HumData, matches !=0)
not_in_data<-setdiff(Hum$tip.label, data$Species)

#Remove species in tree not in the dataset
cleantree<-drop.tip(Hum, not_in_data)
plot(cleantree)

#check for polytomies
is.binary(cleantree)

#The multi2di command will break polytomies in random order.
Humtree2 <- multi2di(cleantree, random=TRUE)
#in case there are singles
Humtree2<-collapse.singles(Humtree2,root.edge=TRUE)

#check for polytomies again
is.binary(Humtree2)
plot(Humtree2)

#check if names match again
name.check(Humtree2, data, data.names = data$Species)

# save the relevant subset of data
data_feeding <- data[c(1,17,18)]

# save the clinging column and names separately to use it later
clinging <- data$Clinging
# zeroes are not allowed in the model so have to convert to clinger vs. legitimate feeder
clinging<-sapply(clinging,function(x) if(x==0) "legit" else "clinger")
species<-data$Species
names(clinging)<-species

# following procedure from http://www.phytools.org/Cordoba2017/ex/6/Discrete-char-models.html
# first, fitting the simple equal rates (ER) model
#fitER<-fitDiscrete(Humtree2,clinging,model="ER")
#fitER
# fitting the ARD model
#fitARD<-fitDiscrete(Humtree2,clinging,model="ARD")
#fitARD
# using AIC to select the best model
#aicc<-setNames(
#  c(fitER$opt$aicc,fitARD$opt$aicc),
#  c("ER","ARD"))
#aic.w(aicc)
# the SYM model here is the same as ER (?); since there are only two states

# using fitMk instead of geiger's fitDiscrete
fitER<-fitMk(Humtree2,clinging,model="ER")
fitER
fitARD<-fitMk(Humtree2,clinging,model="ARD")
fitARD

# using AIC to select the best model
AIC<-setNames(sapply(list(fitER,
                          fitARD),AIC),c("ER","ARD"))
AIC
aic.w(AIC)
# AIC has slightly higher (worse) value for ARD and weight of evidence is also in favor of ER.

# alternative: likelihood ratio test
#k0<-length(fitER$rates)
#k1<-length(fitARD$rates)
#LR<--2*(logLik(fitER)-logLik(fitARD))
#LR
#P_chisq<-pchisq(LR,df=k1-k0,lower.tail=FALSE)
#P_chisq

# next, trying stochastic mapping simulations in phytools
# goal is to plot a tree with ancestral state probabilities at nodes/branches
# following http://blog.phytools.org/search?q=stochastic+mapping
# and http://blog.phytools.org/2019/07/stochastic-character-mapping-with.html
# the line below can be changed to implement a different model instead of "ER"
trees<-make.simmap(Humtree2,clinging,nsim=100,model="ER")
trees
#fmode<-setNames(clinging,
#                names(clinging))
cols<-setNames(c("black","yellow"), c("legit","clinger"))
plot(summary(trees), colors=cols, fsize=0.3, ftype="i",cex=c(0.2,0.1))
legend(x="bottomleft",legend=c("hover", "cling"),pt.bg=cols,pt.cex=2.4,pch=21)

#fan-type tree can be made if preferred aesthetically
#plot(summary(trees), colors=cols, type="fan")
#legend(x="bottomleft",legend=c("legit", "clinger"),pt.bg=cols,pt.cex=2.4,pch=21)

#density map
object<-densityMap(trees,states=names(cols),plot=FALSE)
object<-setMap(object,c("black","yellow"))
plot(object,lwd=c(2,6),ftype="off")

# testing more complex models according to https://rdrr.io/cran/corHMM/f/inst/doc/corHMMv2.1-vignette.pdf
require(ape)
require(expm)
require(corHMM)

phy <- Humtree2
# here we'll use the data_feeding set
plot(phy, show.tip.label = FALSE)
# the following will plot feeding styles onto the tips
tiplabels(pch = 16, col = data_feeding[,2], cex = 0.5)
# or could plot the actual feeding styles
plot(phy, show.tip.label = FALSE)
tiplabels(pch = 16, col = data_feeding[,3], cex = 0.5)

# trying to do the simplest model fitting
# only going to do it for the clinging character
data_clinging <- data_feeding[c(1,2)]
MK_2state <- corHMM(phy = phy, data = data_clinging, rate.cat = 1)
MK_2state
# nifty visualization of the model
# 1 is legitimate feeder, 2 is clinger
plotMKmodel(MK_2state)

# trying models with multiple rate categories/ rate heterogeneity
MK_2cat <- corHMM(phy = phy, data = data_clinging, rate.cat = 2, model = "SYM", get.tip.states = TRUE)
MK_2cat
plotMKmodel(MK_2cat)

#MK_3cat <- corHMM(phy = phy, data = data_clinging, rate.cat = 3, model = "SYM", get.tip.states = TRUE)
#MK_3cat
#plotMKmodel(MK_3cat)

#MK_4cat <- corHMM(phy = phy, data = data_clinging, rate.cat = 4, model = "SYM", get.tip.states = TRUE)
#MK_4cat
#plotMKmodel(MK_4cat)

# AIC to see which model (number of rate categories) fits best
# the simplest model has the best (lowest AIC), the 4 rate model has the highest AICc but 2-4 are very close.
# simple one rate simmap:

phy = MK_2state$phy
datasim = MK_2state$data
model = MK_2state$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
# run get simmap (can be plotted using phytools)
# simulating trees 100x
simmap <- makeSimmap(tree=phy, data=datasim, model=model, rate.cat=1, nSim=100, nCores=1)
cols1<-setNames(c("black", "yellow"),
                c("0", "1"))
plot(summary(simmap), colors = cols1, cex = c(0.2, 0.1), fsize=0.3, ftype="i")
legend(x="bottomleft",legend=c("hover", "cling"),pt.bg=cols1,pt.cex=2.4,pch=21)
#density map
object1<-densityMap(simmap,states=names(cols1),plot=FALSE)
object1<-setMap(object,c("black","yellow"))
plot(object1,lwd=c(2,6),ftype="off")


# 2cat version:
phy = MK_2cat$phy
datasim = MK_2cat$data
model = MK_2cat$solution
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)
cols2<-setNames(c("grey50", "black", "yellow", "#CCCC00"),
               c("(0,R1)","(0,R2)","(1,R1)","(1,R2)"))
# run get simmap (can be plotted using phytools)
simmap <- makeSimmap(tree=phy, data=datasim, model=model, rate.cat=2, nSim=100, nCores=1)

plot(summary(simmap), colors = cols2, cex=c(0.2, 0.1), fsize=0.3, ftype="i")
legend(x="bottomleft",legend=c("hoverR1","hoverR2","clingR1","clingR2"),pt.bg=cols2,pt.cex=2.4,pch=21)
# density can be plotted here but doesn't make much sense because there are two ways of
# being a clinger and two ways of being a legitimate feeder (rate 1, rate 2)
object2<-densityMap(simmap,states=names(cols2),plot=FALSE)
object2<-setMap(object,c("grey50", "black", "yellow", "#CCCC00"))
plot(object2,lwd=c(2,6),ftype="off")
