#This code presents a least squares algorithm that selects narrowband indices to predict
#tree growth from hyperspectral and LiDAR-derived elevation data.

#PRESS (predicted residual error sum of squares) statistic
#represents leave-one-out cross-validation: analytical derivation for least squares regression
PRESS<-function (x) {sum(resid(x)^2/(1 - lm.influence(x)$hat)^2)}

#data:
#bandsAGG=hyperspectral data as a matrix with rows = to plots and columns = to number of narrowbands
#elevation=topographic data derived from LiDAR
#gro=tree diameter growth data, aggregated to the level of single-species plots
#wave=the wavelengths of the hyperspectral data, standardized by centering around mean and dividing by two sd

colnames(bandsAGG)=wave

#Function to iterate over all possible bands
LOOyou<-function (vars) {
  go<-matrix(NA,nrow=ncol(bandsAGG),ncol=ncol(bandsAGG)) #a matrix representing all possible combos of narrowbands
  colnames(go)<-wave
  rownames(go)<-wave

  goR<-matrix(NA,nrow=ncol(bandsAGG),ncol=ncol(bandsAGG))
  colnames(goR)<-wave
  rownames(goR)<-wave


  for(i in 1:ncol(bandsAGG)) {
    for(j in 1:ncol(bandsAGG)) {
      
      #this calculates narrowband indices as the normalized difference between each narrowband combo
      gulp<-((bandsAGG[,i]-bandsAGG[,j])/(bandsAGG[,i]+bandsAGG[,j]))

      #this puts the narrowband index into a regression formula
      form1<-as.formula(paste("gro"," ~ ",vars,"+","gulp",sep = ""))

      #linear regression with a particular narrowband index
      m1<-lm(form1)

      #each element of the matrix is equal to the predictive power (either R2 or PRESS statistic) of that narrowband combo
      go[i,j]<-PRESS(m1)
      goR[i,j]<-summary(m1)$r.squared
      #  print(c(i,j))
    }
  }
  return(list(Presto=go,Resto=goR))}

#the number of iterations of the least squares algorithm should be decided on using a randomization test
#in our case, the randomization test suggested 6 narrowband indices 
n_iter=6


#a list of all the matrices produced for each iteration for PRESS statistic:
MATLe<-vector("list",6)

#a list of all the matrices produced for each iteration for R2 statistic:
MATLeR<-vector("list",6)

#the list of formulas for each iteration:
nextoe<-vector("list",6)

#initialize with elevation as the first covariate:
nextoe[1]<-"elevation"

for(i in 1:6) {

  #runs the LOOyou function above for each iteration
  MATfriends<-LOOyou(nextoe[i])

  MATLe[[i]]<-MATfriends[[1]]
  MATLeR[[i]]<-MATfriends[[2]]

  use1<-which(MATLe[[i]]==min(MATLe[[i]]),arr.ind=T)[1]
  use2<-which(MATLe[[i]]==min(MATLe[[i]]),arr.ind=T)[2]
  hopes<-((bandsAGG[,use1]-bandsAGG[,use2])/(bandsAGG[,use1]+bandsAGG[,use2]))
  assign(paste("vars",use1,use2,sep="_"),hopes)
  nextoe[i+1]<-paste(nextoe[i],paste("vars",use1,use2,sep="_"),sep="+")
  print(nextoe[i+1])
}

#save the results of each iteration: this can take a while to run!
#save(nextoe,file="nextoe.Rdata")
#save(MATLe,file="MATLe.Rdata")
#save(MATLeR,file="MATLeR.Rdata")

#plotting and visualizing final results:


#use the final iteration to extract variables for use:
vars=strsplit(nextoe[[7]],"+",fixed=T)

#extract the individual narrowbands in each of the narrowband indices.
#we subtract the first element from the list because it is elevation, not a narrowband index
USE=strsplit(vars[[1]][-1],"_",fixed=T)


u1=as.numeric(unlist(lapply(USE,"[",2)))
u2=as.numeric(unlist(lapply(USE,"[",3)))

use.vars_ABS<-matrix(nrow=length(gro),ncol=length(u1))
for(i in 1:length(u1)) {
use.vars_ABS[,i]<-((bandsAGG[,u1[i]]-bandsAGG[,u2[i]])/(bandsAGG[,u1[i]]+bandsAGG[,u2[i]]))
}

dat=data.frame(gro,elevation,use.vars_ABS)

#best model from the least squares algorithm
#can look at R-squared from the summary of this model,
#but better to use leave-one-out cross-validation to assess model fit
bestMODELnow=lm(gro~.,data=dat)

