
##############################################
#### Random regression in BLUPF90 ############
#####      Mark Cheng       ##################
##############################################
default_lib_path <- .libPaths()
.libPaths(unique(c("/mnt/custom_lib/", default_lib_path)))

plot_dir<- "/mnt/gnuog/GxE/"

#install.packages("plot3D")
#install.packages("GGally")
#install.packages("ggplot2")
library(ggplot2)
library(GGally)
library(readxl)
library(tidyverse)
library(data.table)
library(orthopolynom) # install.packages("orthopolynom")

rm(list=ls())

# GxE for YLD, using RMAAS as covariate

# All data
data<- read.table("/mnt/GWS_training/ALLIND_data_recode.txt", header = T, sep=" ", na="-999")

# data from submarket K
test<- data[data$GROW_YEAR< 2021 & data$recode %in% genoK$V1,]
table(test$GROW_YEAR)
table(test$MG)
plot(test$RMAAS,test$YLD_BE)
test<- test[,c("recode", "mean", "YR_LOCSET", "RMAAS", "YLD_BE")]

# add Legendre polynomials (0,1,2) (BLUPF90 manual, random regression model)
leg2coef <- legendre.polynomials(n=2, normalized=TRUE)

#standardize the RMAAS
min<- min(test$RMAAS, na.rm = T)
max<- max(test$RMAAS, na.rm = T)
test$RMAAS_sd<- -1+2*(test$RMAAS-min)/(max-min)

M <- matrix(c(rep(1, nrow(test)),test$RMAAS_sd, test$RMAAS_sd^2), ncol=3)
A<- matrix(c(0.7071,0,-0.7906,0,1.2247,0,0,0,2.3717), nrow=3, byrow = T)
pol<- as.data.frame(M %*% A)
names(pol)<- c("poly0","poly1","poly2")

test<- cbind(test,pol)
write.table(test, file="/mnt/gnuog/GxE/Reaction_norm/trainK_RN.dat", quote = F, sep=" ", na="-999", col.names = T, row.names = F)

val<- data[data$GROW_YEAR == 2021 & data$recode %in% genoK$V1,]
val<- val[!is.na(val$recode),]
table(val$GROW_YEAR)
table(val$MG)
val<- val[,c("recode", "mean", "YR_LOCSET", "RMAAS", "YLD_BE")]
write.table(val, file="/mnt/gnuog/GxE/Reaction_norm/valK_RN.dat", quote = F, sep=" ", na="-999", col.names = T, row.names = F)


# run gibbs sampling
setwd("/mnt/gnuog/GxE/Reaction_norm/SubK")

# renumber parameter file
system("nohup /mnt/blupf90/renumf90 Random_Regre.par > renum.log &")

# run renumbered parameter file
system("ulimit -s unlimited
          export OMP_STACKSIZE=64M
          niter=100000 ;   nbin=10000 ;  nsp=100
          { echo renf90.par ; echo $niter; echo $nbin; echo $nsp; } | nohup /mnt/blupf90/gibbs3f90 > gibbsf90.log &
         ")




# prediction results

# liner qudratic Legendre polynomials
# check convergence
post<- fread("/mnt/gnuog/GxE/Reaction_norm/SubK/postgibbs_samples")
ggplot(post, 
       aes(x=V1, y=V16)) +
        geom_line(color="red", aes(alpha=0.5)) +
        ggtitle("Trace plot") +
        xlab("Iteration") +
        ylab("Genetic parameter")

train<- read.table("/mnt/gnuog/GxE/Reaction_norm/trainK_RN.dat", header = T, sep=" ", na="-999")
G1<- matrix(c(9.3138, 4.2813 , 4.4904,  4.2813, 3.2723, 1.6134, 4.4904 , 1.6134 , 2.3689), nrow = 3, byrow = T)
G2<- matrix(c(1.7284 ,0.14104E-01 ,1.8829, 0.14104E-01,1.8034, 0.15899, 1.8829 , 0.15899 , 2.5488), nrow = 3, byrow = T)
R= 12.321

RM<- train[sample(1:nrow(train),1000),]
RM<- RM[!is.na(RM$RMAAS),]
hist(RM$RMAAS)

Z<- as.matrix(RM[,c("poly0", "poly1","poly2")])
Vg<- diag(Z %*% G1 %*% t(Z))
Vp<- diag(Z %*% G2 %*% t(Z))
h2<- Vg/(Vg+Vp+R)

RM<- cbind(RM, h2)
plot(RM$RMAAS, RM$h2)


library(plot3D)

plot_dir<- "/mnt/gnuog/GxE/"

RM<- RM[order(RM$RMAAS),]
RM<- RM[!duplicated(RM$RMAAS),]
Z<- as.matrix(RM[,c("poly0", "poly1","poly2")])
rg<- cov2cor(Z %*% G1 %*% t(Z))

x<- RM$RMAAS
y<- RM$RMAAS

pdf(paste0(plot_dir, 
           "3D plot YLD on RMAAS SubK.pdf"), 
    width=8, height=6)

print(persp3D(x,-y,rg,theta=30, phi=50, axes=TRUE,scale=2, box=TRUE, nticks=5, 
              
              ticktype="detailed",xlab="RMAAS", ylab="RMAAS", zlab="rg", 
              
              main="genetic correlation of YLD_BE among RMAAS"))

dev.off()




# prediction results SubK
val<- read.table("/mnt/gnuog/GxE/Reaction_norm/valK_RN.dat", header = T, sep=" ", na="-999")
val<- val[!is.na(val$YLD_BE),]

sol<- fread("/mnt/gnuog/GxE/Reaction_norm/SubK/final_solutions.txt") # final_solutions.txt file can be obtained running binsol_to_textsol.f90
ebv1<- sol%>%filter(effect ==6)
ebv2<- sol%>%filter(effect ==7)
ebv3<- sol%>%filter(effect ==8)

ped<- read.table("/mnt/gnuog/GxE/Reaction_norm/SubK/renadd06.ped", header=F)
ped<- ped[,c(1,10)]
names(ped)<- c("level","recode")
ebv1<- left_join(ebv1, ped, by="level")
ebv1<- ebv1[,c(6,4)]
names(ebv1)<- c("recode","bv1")

ebv2<- left_join(ebv2, ped, by="level")
ebv2<- ebv2[,c(6,4)]
names(ebv2)<- c("recode","bv2")

ebv3<- left_join(ebv3, ped, by="level")
ebv3<- ebv3[,c(6,4)]
names(ebv3)<- c("recode","bv3")

val<- left_join(val, ebv1, by="recode")
val<- left_join(val, ebv2, by="recode")
val<- left_join(val, ebv3, by="recode")

# add Legendre polynomials (0,1,2) (BLUPF90 manual, random regression model)
leg2coef <- legendre.polynomials(n=2, normalized=TRUE)

#standardize the RMAAS
min<- min(val$RMAAS, na.rm = T)
max<- max(val$RMAAS, na.rm = T)
val$RMAAS_sd<- -1+2*(val$RMAAS-min)/(max-min)

M <- matrix(c(rep(1, nrow(val)),val$RMAAS_sd, val$RMAAS_sd^2), ncol=3)
A<- matrix(c(0.7071,0,-0.7906,0,1.2247,0,0,0,2.3717), nrow=3, byrow = T)
pol<- as.data.frame(M %*% A)
names(pol)<- c("poly0","poly1","poly2")
val<- cbind(val,pol)

val$ebv<- val$bv1*val$poly0+val$bv2*val$poly1+val$bv3*val$poly2

val2<- val[val$RMAAS> 2 & val$RMAAS< 3,]
ggplot(val2, aes(RMAAS, ebv))+ geom_line(aes(group=recode))

cor(val$YLD_BE, val$ebv)
mean<- val%>%group_by(recode)%>%summarise(meanY=mean(YLD_BE), meanBV=mean(ebv))
cor(mean$meanY, mean$meanBV)







