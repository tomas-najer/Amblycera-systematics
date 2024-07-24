library(ape)
#library(phyloch)
library(strap)
library(phylotate)
library(MCMCtreeR)
library(phytools)
library(geiger)
library(ggtree)
library(caper)
library(ggplotify)
library(treeio)
library("readxl")
library("ape")
library("ggplot2")
library("phytools")
library("geiger")
library(ggplot2)
library(gggenes)
library(ggfittext)
library(ggrepel)
library(dplyr)
library("ggtree")
library(ggplotify)

setwd('/home/tom/Documents/Amblycera_systematics/ancestral_hosts')
louse_tree <- read.nexus("amblycera_dated.nex")

plot.phylo(louse_tree,cex=0.4,no.margin=T)
nodelabels()
louse_tree <- extract.clade(phy=louse_tree, 149)

louse_tree$root.time <- 68.3596 #root age
plot.phylo(louse_tree,cex=0.4,no.margin=T)
x <- read.csv("tree_labels.csv")
louse_tree <- rename_taxa(louse_tree, x, label, label2)
plot.phylo(louse_tree,cex=0.4,no.margin=T)

geoscalePhylo(tree=louse_tree, units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, lwd=3, width=2,)

hosts <- read.csv('hosts.csv',header=T)
hosts
arch <- setNames(as.factor(hosts[,2]),hosts[,1])
arch
to.matrix(arch,levels(arch))
name.check(louse_tree,arch)
louse_tree$tip.label
structure_matrix<-to.matrix(arch[louse_tree$tip.label],levels(arch))
structure_matrix
geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,81), lwd=3, width=2,quat.rm = T)
plot.phylo(louse_tree,cex=0.5)
tiplabels(pie=structure_matrix,piecol=c("orange","blue"),cex=0.15)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2.2)

ace(arch,louse_tree,model='ER',type='discrete')
fitER <- ace(arch,louse_tree,model='ER',type='discrete')
ace(arch,louse_tree,model='ARD',type='discrete')
fitARD <- ace(arch,louse_tree,model='ARD',type='discrete')
ace(arch,louse_tree,model='SYM',type='discrete')
fitSYM <- ace(arch,louse_tree,model='SYM',type='discrete')
fitUSR <- ace(arch, louse_tree, model=matrix(c(0,1,0,0),2),type='discrete')
model_geiger_ER <- fitDiscrete(louse_tree,arch,type='discrete',model="ER")
model_geiger_ARD <- fitDiscrete(louse_tree,arch,type='discrete',model="ARD")
model_geiger_SYM <- fitDiscrete(louse_tree,arch,type='discrete',model="SYM")
model_geiger_USR <- fitDiscrete(louse_tree,arch,type='discrete',model=matrix(c(0,1,0,0)))

ERvSYM <- abs(model_geiger_ER$opt$aicc - model_geiger_SYM$opt$aicc)
ERvARD <- abs(model_geiger_ER$opt$aicc - model_geiger_ARD$opt$aicc)
ERvUSR <- abs(model_geiger_ER$opt$aicc - model_geiger_USR$opt$aicc)
SYMvARD <- abs(model_geiger_SYM$opt$aicc - model_geiger_ARD$opt$aicc)

comp <- data.frame(model = c("ER","SYM","ARD","USR"),AIC = c(model_geiger_ER$opt$aic,model_geiger_SYM$opt$aic,model_geiger_ARD$opt$aic,model_geiger_USR$opt$aic),
                   AICc = c(model_geiger_ER$opt$aicc,model_geiger_SYM$opt$aicc,model_geiger_ARD$opt$aicc,model_geiger_USR$opt$aicc))
comp

1-pchisq(2*abs(fitER$loglik - fitARD$loglik), 1)
1-pchisq(2*abs(fitER$loglik - fitUSR$loglik), 1)
1-pchisq(2*abs(fitARD$loglik - fitUSR$loglik), 1)

round(fitER$lik.anc,3)
round(fitARD$lik.anc,3)
round(fitUSR$lik.anc,3)

geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,81), lwd=3, width=2,quat.rm = T)
tiplabels(pie=structure_matrix,piecol=c("orange","blue"),cex=0.15)
nodelabels(node=1:louse_tree$Nnode+Ntip(louse_tree),
           pie=fitER$lik.anc,piecol=c("orange","blue"),cex=0.3)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2)

?make.simmap
mtrees<-make.simmap(louse_tree,arch,model="ER",nsim=1000)
summary(mtrees)
pd<-summary(mtrees,plot=T)

plot(mtrees[[10]],cols,fsize=0.6,ftype="i")
geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,81), lwd=3, width=2,quat.rm = T)
plot.phylo(louse_tree,cex=0.4,no.margin=T)
tiplabels(pie=structure_matrix,piecol=c("orange","blue"),cex=0.15)
nodelabels(pie=pd$ace,piecol=c("orange","blue"),cex=0.3)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2)

# try it again with coalescent tree
coalescent <- read.nexus("coalescent.nex")
plot.phylo(coalescent,cex=0.4,no.margin=T)
coalescent <- drop.tip(coalescent,"Risp")
x <- read.csv("tree_labels.csv")
coalescent <- rename_taxa(coalescent, x, label, label2)
plot.phylo(coalescent,cex=0.4,no.margin=T)
nodelabels()
coalescent <- extract.clade(phy=coalescent, 149)
plot.phylo(coalescent,cex=0.4,no.margin=T)
coalescent_matrix<-to.matrix(arch[coalescent$tip.label],levels(arch))
tiplabels(pie=coalescent_matrix,piecol=c("orange","blue"),cex=0.15)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2.2)

ace(arch,coalescent,model='ER',type='discrete')
fitER <- ace(arch,coalescent,model='ER',type='discrete')
ace(arch,coalescent,model='ARD',type='discrete')
fitARD <- ace(arch,coalescent,model='ARD',type='discrete')
ace(arch,coalescent,model='SYM',type='discrete')
fitSYM <- ace(arch,coalescent,model='SYM',type='discrete')
fitUSR <- ace(arch, coalescent, model=matrix(c(0,1,0,0),2),type='discrete')
model_geiger_ER <- fitDiscrete(coalescent,arch,type='discrete',model="ER")
model_geiger_ARD <- fitDiscrete(coalescent,arch,type='discrete',model="ARD")
model_geiger_SYM <- fitDiscrete(coalescent,arch,type='discrete',model="SYM")
model_geiger_USR <- fitDiscrete(coalescent,arch,type='discrete',model=matrix(c(0,1,0,0)))

ERvSYM <- abs(model_geiger_ER$opt$aicc - model_geiger_SYM$opt$aicc)
ERvARD <- abs(model_geiger_ER$opt$aicc - model_geiger_ARD$opt$aicc)
ERvUSR <- abs(model_geiger_ER$opt$aicc - model_geiger_USR$opt$aicc)
SYMvARD <- abs(model_geiger_SYM$opt$aicc - model_geiger_ARD$opt$aicc)

comp <- data.frame(model = c("ER","SYM","ARD","USR"),AIC = c(model_geiger_ER$opt$aic,model_geiger_SYM$opt$aic,model_geiger_ARD$opt$aic,model_geiger_USR$opt$aic),
                   AICc = c(model_geiger_ER$opt$aicc,model_geiger_SYM$opt$aicc,model_geiger_ARD$opt$aicc,model_geiger_USR$opt$aicc))
comp

#AIC weights
aic.all<-setNames(c(AIC(fitER),AIC(fitSYM),AIC(fitARD),AIC(fitUSR)),
                  c("ER","SYM","ARD","USR"))
aic.w(aic.all)

plot.phylo(coalescent,cex=0.4,no.margin=T)
tiplabels(pie=coalescent_matrix,piecol=c("orange","blue"),cex=0.15)
nodelabels(node=1:coalescent$Nnode+Ntip(coalescent),
           pie=fitARD$lik.anc,piecol=c("orange","blue"),cex=0.3)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2)

coaltrees<-make.simmap(coalescent,arch,model="ARD",nsim=1000)
summary(coaltrees)
coal_pd<-summary(coaltrees,plot=T)

plot.phylo(coalescent,cex=0.4,no.margin=T)
tiplabels(pie=coalescent_matrix,piecol=c("orange","blue"),cex=0.15)
nodelabels(pie=coal_pd$ace,piecol=c("orange","blue"),cex=0.3)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2)

# make nice tree
?plot.phylo
?tiplabels
plot.phylo(coalescent,type="fan",cex=0.5,label.offset=0.3,no.margin=F)
tiplabels(pie=coalescent_matrix,piecol=c("orange","#00BFFF"),cex=0.2)
nodelabels(pie=coal_pd$ace,piecol=c("orange","#00BFFF"),cex=0.2)
legend("bottom",title="Host",levels(arch),pch=21,pt.bg=c("orange","#00BFFF"),
       pt.cex=2.2)

geoscalePhylo(tree=ladderize(coalescent,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, lwd=3, width=2,quat.rm = T)
MCMC.tree.plot(louse_tree,ladderize.tree=F)

# try it with mcmc tree
phy <- readMCMCtree("FigTree_Amblycera.tre",from.file=T)
phy
phy$apePhy
plot.phylo(phy$apePhy,cex=0.5)
d <- read.csv("tree_labels.csv")
d
phy$apePhy <- rename_taxa(phy$apePhy, d, label, label2)
plot.phylo(phy$apePhy,cex=0.5)
mcmc_tree <- phy$apePhy
mcmc_tree$tip.label
mcmc_matrix<-to.matrix(arch[mcmc_tree$tip.label],levels(arch))

?MCMC.tree.plot
MCMC.tree.plot(phy,analysis.type="MCMCtree",cex.tips=0.4, 
               time.correction=100,plot.type="phylogram",lwd.bar=20, 
               scale.res=c("Period","Epoch"),node.method="full.length", 
               all.nodes=c(161,164),col.age="#00BFFF40",no.margin=TRUE,label.offset=1)
nodelabels()
?tiplabels
mcmc_matrix
tiplabels(pie=mcmc_matrix,piecol=c("orange","blue"),cex=0.15)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2.2)

ace(arch,mcmc_tree,model='ER',type='discrete')
fitER <- ace(arch,mcmc_tree,model='ER',type='discrete')
ace(arch,mcmc_tree,model='ARD',type='discrete')
fitARD <- ace(arch,mcmc_tree,model='ARD',type='discrete')
ace(arch,mcmc_tree,model='SYM',type='discrete')
fitSYM <- ace(arch,mcmc_tree,model='SYM',type='discrete')
fitUSR <- ace(arch, mcmc_tree, model=matrix(c(0,1,0,0),2),type='discrete')
model_geiger_ER <- fitDiscrete(mcmc_tree,arch,type='discrete',model="ER")
model_geiger_ARD <- fitDiscrete(mcmc_tree,arch,type='discrete',model="ARD")
model_geiger_SYM <- fitDiscrete(mcmc_tree,arch,type='discrete',model="SYM")
model_geiger_USR <- fitDiscrete(mcmc_tree,arch,type='discrete',model=matrix(c(0,1,0,0)))

ERvSYM <- abs(model_geiger_ER$opt$aicc - model_geiger_SYM$opt$aicc)
ERvARD <- abs(model_geiger_ER$opt$aicc - model_geiger_ARD$opt$aicc)
ERvUSR <- abs(model_geiger_ER$opt$aicc - model_geiger_USR$opt$aicc)
SYMvARD <- abs(model_geiger_SYM$opt$aicc - model_geiger_ARD$opt$aicc)

comp <- data.frame(model = c("ER","SYM","ARD","USR"),AIC = c(model_geiger_ER$opt$aic,model_geiger_SYM$opt$aic,model_geiger_ARD$opt$aic,model_geiger_USR$opt$aic),
                   AICc = c(model_geiger_ER$opt$aicc,model_geiger_SYM$opt$aicc,model_geiger_ARD$opt$aicc,model_geiger_USR$opt$aicc))
comp

aic.all<-setNames(c(AIC(fitER),AIC(fitSYM),AIC(fitARD),AIC(fitUSR)),
                  c("ER","SYM","ARD","USR"))
#AIC weights
aic.w(aic.all)
aic.w(model_geiger_ER$opt$aic)
aic.w(model_geiger_SYM$opt$aic)


1-pchisq(2*abs(fitER$loglik - fitARD$loglik), 1)
1-pchisq(2*abs(fitER$loglik - fitUSR$loglik), 1)
1-pchisq(2*abs(fitARD$loglik - fitUSR$loglik), 1)

round(fitER$lik.anc,3)
round(fitARD$lik.anc,3)
round(fitUSR$lik.anc,3)

MCMC.tree.plot(phy,analysis.type="MCMCtree",cex.tips=0.4, 
               time.correction=100,plot.type="phylogram",lwd.bar=2, 
               scale.res=c("Period","Epoch"),node.method="none", 
               col.age="navy",no.margin=TRUE,label.offset=1)
nodelabels(node=1:mcmc_tree$Nnode+Ntip(mcmc_tree),
           pie=fitER$lik.anc,piecol=c("orange","blue"),cex=0.3)
nodelabels()
tiplabels(pie=mcmc_matrix,piecol=c("orange","blue"),cex=0.15)
legend("topleft",levels(arch),pch=21,pt.bg=c("orange","blue"),
       pt.cex=2.2)

mcmctrees<-make.simmap(mcmc_tree,arch,model="ER",nsim=1000)
summary(mcmctrees)
mcmc_pd<-summary(mcmctrees,plot=T)

?MCMC.tree.plot

#MCMC.tree.plot(phy,analysis.type="MCMCtree",cex.tips=0.4, 
#               time.correction=100,plot.type="phylogram",lwd.bar=2, 
#               scale.res=c("Period","Epoch"),node.method="full.length", 
#               all.nodes=c(161,164),col.age="#ff000040",no.margin=TRUE,label.offset=1)

MCMC.tree.plot(phy,analysis.type="MCMCtree",cex.tips=0.5, 
               time.correction=100,plot.type="phylogram",lwd.bar=2, 
               scale.res=c("Period","Epoch"),node.method="none", 
               all.nodes=c(161,164),col.age="#ff000040",no.margin=TRUE,label.offset=1)

#phy
#nodelabels()
nodelabels(pie=mcmc_pd$ace,piecol=c("orange","#00BFFF"),cex=0.3)
#nodelabels()
tiplabels(pie=mcmc_matrix,piecol=c("orange","#00BFFF"),cex=0.15)
#?legend
legend("left",inset=0.05,title="Host",levels(arch),pch=21,pt.bg=c("orange","#00BFFF"),
       pt.cex=2.2)

# mcmctree with mammals
full_mcmc <- readMCMCtree("FigTree.tre",from.file=T)
full_mcmc
full_mcmc$apePhy
plot.phylo(full_mcmc$apePhy,cex=0.5)
nodelabels()
full_mcmc$apePhy <- drop.tip(full_mcmc$apePhy,"Risp")
MCMC.tree.plot(full_mcmc,analysis.type="MCMCtree",cex.tips=0.4, 
               time.correction=100,plot.type="phylogram",lwd.bar=20, 
               scale.res=c("Period","Epoch"),node.method="full.length", 
               all.nodes=c(194,197),col.age="#00BFFF40",no.margin=TRUE,label.offset=1)
