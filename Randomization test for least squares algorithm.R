#Code for a randomization test to quantify the probability that our observed
#values of the PRESS statistic from the least squares algorithm were due to random chance. 
#This code reruns the entire least squares algorithm for 1000 iterations, randomizing the order of the standardized
#growth data at each iteration, while retaining the original predictor matrix of hyperspectral and
#elevation data. During each iteration of the randomization test, 
#we allowed the least squares algorithm to select up to ten narrowband indices.
#The number of indices to select is up to the user, but be warned: this proceduce is computationally-intensive!

sample_size=1000


#data:
#bandsAGG=hyperspectral data as a matrix with rows = to plots and columns = to number of narrowbands
#elevation=topographic data derived from LiDAR
#gro=tree diameter growth data, aggregated to the level of single-species plots
#wave=the wavelengths of the hyperspectral data, standardized by centering around mean and dividing by two sd


#PRESS (predicted residual error sum of squares) statistic
#represents leave-one-out cross-validation: analytical derivation for least squares regression
PRESS<-function (x) {sum(resid(x)^2/(1 - lm.influence(x)$hat)^2)}


colnames(bandsAGG)=wave

ALLtest<-vector("list",length(1000)) # a list of randomized growth data

ALLtest[[1]]=gro #the first element of the list is the growth data with the real order

for(i in 2:sample_size){ #randomly sampled growth data. 
  ALLtest[[i]]<-sample(groAGG$GROWTH_std)
  }


#function to iterate over all possible narrowband indices and calculate R2 and PRESS statistic
LOOyou<-function (itera,vars) {
  
  gmat<-bandsAGG
  
  gro=ALLtest[[itera]]
  
  go<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
  colnames(go)<-paste("c",c(1:ncol(gmat)),sep="")
  rownames(go)<-paste("c",c(1:ncol(gmat)),sep="")
  
  goR<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
  colnames(goR)<-paste("c",c(1:ncol(gmat)),sep="")
  rownames(goR)<-paste("c",c(1:ncol(gmat)),sep="")
  
  for(i in 1:ncol(gmat)) {
    for(j in 1:ncol(gmat)) {
      gulp<-((gmat[,i]-gmat[,j])/(gmat[,i]+gmat[,j]))
      
      form1<-as.formula(paste("gro"," ~ ",vars,"+","gulp",sep = ""))
      
      m1<-lm(form1)
      
      
      go[i,j]<-PRESS(m1)
      goR[i,j]<-summary(m1)$r.squared
      #  print(c(i,j))
    }
  }
  return(list(Presto=go,Resto=goR))}


#vectors to store R-squared and PRESS statistic for each randomized iteration

ultiR<-vector("list",length(ALLtest))
ultiP<-vector("list",length(ALLtest))

for(j in 1:length(ALLtest)){

  MATLe<-vector("list",10) #again, we are letting it select up to 10 narrowband indices

  MATLeR<-vector("list",10)

  nextoe<-vector("list",10)
  nextoe[1]<-"elevation"

  for(i in 1:10) {
    MATfriends<-LOOyou(j,nextoe[i])
    MATLe[[i]]<-MATfriends[[1]]
    MATLeR[[i]]<-MATfriends[[2]]

    use1<-which(MATLe[[i]]==min(MATLe[[i]],na.rm=T),arr.ind=T)[1]
    use2<-which(MATLe[[i]]==min(MATLe[[i]],na.rm=T),arr.ind=T)[2]
    hopes<-((bandsAGG[,use1]-bandsAGG[,use2])/(bandsAGG[,use1]+bandsAGG[,use2]))
    assign(paste("vars",use1,use2,sep="_"),hopes)
    nextoe[i+1]<-paste(nextoe[i],paste("vars",use1,use2,sep="_"),sep="+")
    save(nextoe,file=paste("nextoe",j,".Rdata",sep=""))
    #print(i)
    if(i>1) {
      if(min(MATLe[[i]],na.rm=T)>=min(MATLe[[i-1]],na.rm=T)) break}
    else{next}

  }

  MATLe<-MATLe[c(1:i)]
  MATLeR<-MATLeR[c(1:i)]

  ultiP[[j]]<-unlist(lapply(MATLe,min,na.rm=T))
    getters<-unlist(lapply(MATLeR,max,na.rm=T))

  ultiR[[j]]<-getters
 print(j)

 save(ultiR,file="ultiR.Rdata") #it is important to save this at each iteration, because the code takes a while to run
 save(ultiP,file="ultiP.Rdata")
}

#save the final versions just to be sure
save(ultiR,file="ultiR.Rdata")
save(ultiP,file="ultiP.Rdata")
save(nextoe[[1]],file="")

load("ultiR.Rdata")
load("ultiP.Rdata")
load("nextoeFINAL.Rdata")

#here we are using the results of the randomization test for the PRESS statistic
#to calculate the probability that that value is lower than the 
#PRESS statistic for the real data, which is the first iteration of the randomization algorithm: ultiP[[1]]

gulpo=function(P) {
  vec1=unlist(lapply(ultiP,"[",P))
  return(length(which(vec1<ultiP[[1]][P]))/1000)
  }

#applying the function above to each of the ten slots for narrowband indices:
MYEP=rep(NA,times=10)
for(i in 1:10) {
MYEP[i]<-gulpo(i)
}

#plotting probability of observing a PRESS statistic < reality via random change:
par(mar=c(6,6,4,4))
plot(MYEP~c(1:10),main="",
xlab="Iteration of algorithm",ylab="Probability PRESS statistic due to random chance",pch=21,bg="gray",cex=1.5,
cex.lab=1.4,cex.axis=1.3,ylim=c(0,0.15))
#abline(h=0.05,col="red",lwd=5)

points(MYEP~c(1:10),pch=21,bg="gray",
cex=2.3)       
abline(h=0.05,col="gray30",lwd=12,lty=1)
points(MYEP~c(1:10),pch=21,bg="gray",
       cex=2.9)       
