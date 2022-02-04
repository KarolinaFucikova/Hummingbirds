# calculating Pagel's lambda
library("phytools")
library("parallel")
library("viridis")
library("geiger")
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


#because there are different subsets for the different values, separate the columns as objects
BodyMass <- data$WtSpMean
Wing <- data$WingSpMean
Bill <- data$ECulmSpMean
Tarsus <- data$TarSpMean
Hallux <- data$HClawSpMean
#Midtoe <- data$LogMidToeClaw
# leaving the Midtoe out for now

#data_with_names <- data$Species
#names(data_with_names) <- data_with_names$species_names
#simple, without accounting for within-species variance
names(BodyMass) <- data$Species
phylosig(Humtree2, BodyMass, method="lambda", test=TRUE)

names(Wing) <- data$Species
phylosig(Humtree2, Wing, method="lambda", test=TRUE)

names(Bill) <- data$Species
phylosig(Humtree2, Bill, method="lambda", test=TRUE)

names(Tarsus) <- data$Species
phylosig(Humtree2, Tarsus, method="lambda", test=TRUE)

names(Hallux) <- data$Species
phylosig(Humtree2, Hallux, method="lambda", test=TRUE)

# names(Midtoe) <- data$Species
# phylosig(Humtree2, Midtoe, method="lambda", test=TRUE)


# the following code conducts the same test but also accounts for variance (se) in measurments
# first we must replace zeroes with a small number, 0.000001
data1<-data
data1[data1 == 0] <- 0.000001 

# next, assign the standard error values to new vectors for easier handling
bmse <-data1$WtSpSE
wingse<-data1$WingSpSE
billse<-data1$ECulmSpSE
tarse<-data1$TarSpSE
halse<-data1$HClawSpSE

# conduct the test for each variable respectively
names(BodyMass) <- data$Species
names(bmse)<-data$Species
phylosig(Humtree2, BodyMass, method="lambda", test=TRUE, se=bmse)

names(Wing) <- data$Species
names(wingse)<-data$Species
phylosig(Humtree2, Wing, method="lambda", test=TRUE, se=wingse)

names(Bill) <- data$Species
names(billse)<-data$Species
phylosig(Humtree2, Bill, method="lambda", test=TRUE, se=billse)

names(Tarsus) <- data$Species
names(tarse)<-data$Species
phylosig(Humtree2, Tarsus, method="lambda", test=TRUE, se=tarse)

names(Hallux) <- data$Species
names(halse)<-data$Species
phylosig(Humtree2, Hallux, method="lambda", test=TRUE, se=halse)


#testing if different from Brownian Motion model, gotta remove missing data
BMClean <- na.omit(data1[,c(1,2,4)])
x<-BMClean
not_in_bmclean<-setdiff(Humtree2$tip.label, x$Species)
#Remove species in tree not in the dataset
prunedtree<-drop.tip(Humtree2, not_in_bmclean)
plot(prunedtree)
#The multi2di command will break polytomies in random order.
prunedtree <- multi2di(prunedtree, random=TRUE)
#in case there are singles
prunedtree<-collapse.singles(prunedtree,root.edge=TRUE)
#check for polytomies again
is.binary(prunedtree)
plot(prunedtree)
# fit lambda model
y=x$WtSpMean
z=x$WtSpSE
names(y) <- x$Species
names(z) <- x$Species
lambda<-phylosig(prunedtree, y, method="lambda", se=z)
lambda
# get likelihood for Brownian Motion model
prunedtree$mapped.edge<-matrix(prunedtree$edge.length, nrow(prunedtree$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(prunedtree,y)$logL1
# conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P

#testing if different from Brownian Motion model, gotta remove missing data
WinClean <- na.omit(data1[,c(1,5,7)])
x<-WinClean
not_in_winclean<-setdiff(Humtree2$tip.label, x$Species)
#Remove species in tree not in the dataset
prunedtree<-drop.tip(Humtree2, not_in_winclean)
plot(prunedtree)
#The multi2di command will break polytomies in random order.
prunedtree <- multi2di(prunedtree, random=TRUE)
#in case there are singles
prunedtree<-collapse.singles(prunedtree,root.edge=TRUE)
#check for polytomies again
is.binary(prunedtree)
plot(prunedtree)
# fit lambda model
y=x$WingSpMean
z=x$WingSpSE
names(y) <- x$Species
names(z) <- x$Species
lambda<-phylosig(prunedtree, y, method="lambda", se=z)
lambda
# get likelihood for Brownian Motion model
prunedtree$mapped.edge<-matrix(prunedtree$edge.length, nrow(prunedtree$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(prunedtree,y)$logL1
# conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P

#testing if different from Brownian Motion model, gotta remove missing data
TarsClean <- na.omit(data1[,c(1,11,13)])
x<-TarsClean
not_in_tarsclean<-setdiff(Humtree2$tip.label, x$Species)
#Remove species in tree not in the dataset
prunedtree<-drop.tip(Humtree2, not_in_tarsclean)
plot(prunedtree)
#The multi2di command will break polytomies in random order.
prunedtree <- multi2di(prunedtree, random=TRUE)
#in case there are singles
prunedtree<-collapse.singles(prunedtree,root.edge=TRUE)
#check for polytomies again
is.binary(prunedtree)
plot(prunedtree)
# fit lambda model
y=x$TarSpMean
z=x$TarSpSE
names(y) <- x$Species
names(z) <- x$Species
lambda<-phylosig(prunedtree, y, method="lambda", se=z)
lambda
# get likelihood for Brownian Motion model
prunedtree$mapped.edge<-matrix(prunedtree$edge.length, nrow(prunedtree$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(prunedtree,y)$logL1
# conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P

#testing if different from Brownian Motion model, gotta remove missing data
BillClean <- na.omit(data1[,c(1,8,10)])
x<-BillClean
not_in_billclean<-setdiff(Humtree2$tip.label, x$Species)
#Remove species in tree not in the dataset
prunedtree<-drop.tip(Humtree2, not_in_billclean)
plot(prunedtree)
#The multi2di command will break polytomies in random order.
prunedtree <- multi2di(prunedtree, random=TRUE)
#in case there are singles
prunedtree<-collapse.singles(prunedtree,root.edge=TRUE)
#check for polytomies again
is.binary(prunedtree)
plot(prunedtree)
# it is not necessary to plot the trees every time
# but it may be a useful visual check for any issues with the phylogeny

# fit lambda model
y=x$ECulmSpMean
z=x$ECulmSpSE
names(y) <- x$Species
names(z) <- x$Species
lambda<-phylosig(prunedtree, y, method="lambda", se=z)
lambda
# get likelihood for Brownian Motion model
prunedtree$mapped.edge<-matrix(prunedtree$edge.length, nrow(prunedtree$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(prunedtree,y)$logL1
# conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P


#testing if different from Brownian Motion model, gotta remove missing data
HalClean <- na.omit(data1[,c(1,14,16)])
x<-HalClean
not_in_halclean<-setdiff(Humtree2$tip.label, x$Species)
#Remove species in tree not in the dataset
prunedtree<-drop.tip(Humtree2, not_in_halclean)
plot(prunedtree)
#The multi2di command will break polytomies in random order.
prunedtree <- multi2di(prunedtree, random=TRUE)
#in case there are singles
prunedtree<-collapse.singles(prunedtree,root.edge=TRUE)
#check for polytomies again
is.binary(prunedtree)
plot(prunedtree)
# fit lambda model
y=x$HClawSpMean
z=x$HClawSpSE
names(y) <- x$Species
names(z) <- x$Species
lambda<-phylosig(prunedtree, y, method="lambda", se=z)
lambda
# get likelihood for Brownian Motion model
prunedtree$mapped.edge<-matrix(prunedtree$edge.length, nrow(prunedtree$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(prunedtree,y)$logL1
# conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P


#if there are no missing data it's much easier:
#testing if lambda is significantly different from 1, rather than from 0
x<-Bill
# the above line can be changed to test the different traits
# see lines 50-54 for the names of traits
# fit lambda model
lambda<-phylosig(Humtree2, x, method="lambda")
lambda
# get likelihood for Brownian Motion model
Humtree2$mapped.edge<-matrix(Humtree2$edge.length, nrow(Humtree2$edge),1,dimnames=list(NULL,"1"))
bm.logL<-brownie.lite(Humtree2,x)$logL1
# conduct hypothesis test using chi-square
LR<--2*(bm.logL-lambda$logL)
P<-pchisq(LR,df=1,lower.tail=F)
P
