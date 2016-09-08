#This is code for leave-one-out cross-validation on the least squares algorithm
#to predict tree growth from hyperspectral data.

#data:
#bandsAGG=hyperspectral data as a matrix with rows = to plots and columns = to number of narrowbands
#elevation=topographic data derived from LiDAR
#gro=tree diameter growth data, aggregated to the level of single-species plots
#wave=the wavelengths of the hyperspectral data, standardized by centering around mean and dividing by two sd


sample_size=length(gro)

ALLresp<-vector("list",length=sample_size) #response data with one observation left out
ALLpredB<-vector("list",length=sample_size) #narrowband hyperspectral data with one observation left out
ALLpredE<-vector("list",length=sample_size) #elevation data with one observation left out

testB<-vector("list",length=sample_size) #predictor variables with only one observation (the one that is left out above)
testE<-vector("list",length=sample_size)

#this for loop creates the training and test datasets for LOOCV
for(i in 1:sample_size){
  ALLresp[[i]]<-gro[-i]
  ALLpredB[[i]]<-bandsAGG[-i,]
  ALLpredE[[i]]<-elevation[-i]
  testB[[i]]<-bandsAGG[i,]
  testE[[i]]<-elevation[i]
}



#this function takes in the iteration of the LOOCV and then 
#evaluates predictive power of all possible narrowband indices with R2 and PRESS statistic
LOOyou<-function (itera,vars) {

  gmat<-ALLpredB[[itera]]

  gro=ALLresp[[itera]]

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
           }
  }
  return(list(Presto=go,Resto=goR))}


ultiR<-vector("list",length(ALLresp))
ultiP<-vector("list",length(ALLresp))

#this nested for loop evaluates the least squares algorithm for
#each of the leave one out iterations

for(j in 1:length(ALLresp)){

  MATLe<-vector("list",10)

  MATLeR<-vector("list",10)

  nextoe<-vector("list",10)

  elevation<-ALLpredE[[j]]

  nextoe[1]<-"elevation"

  for(i in 1:6) { #number of iterations determined by randomization test
    MATfriends<-LOOyou(j,nextoe[i])
    MATLe[[i]]<-MATfriends[[1]]
    MATLeR[[i]]<-MATfriends[[2]]

    use1<-which(MATLe[[i]]==min(MATLe[[i]],na.rm=T),arr.ind=T)[1]
    use2<-which(MATLe[[i]]==min(MATLe[[i]],na.rm=T),arr.ind=T)[2]
    hopes<-((ALLpredB[[j]][,use1]-ALLpredB[[j]][,use2])/(ALLpredB[[j]][,use1]+ALLpredB[[j]][,use2]))
    assign(paste("vars",use1,use2,sep="_"),hopes)
    nextoe[i+1]<-paste(nextoe[i],paste("vars",use1,use2,sep="_"),sep="+")
    save(nextoe,file=paste("nextoe",j,".Rdata",sep=""))
    #the save command labels the saved variables with "j" the iteration of LOOCV
      if(i>1) {
      if(min(MATLe[[i]],na.rm=T)>=min(MATLe[[i-1]],na.rm=T)) break}
    else{next}

  }
}

#the output of the nested for loop above is a saved vector of the final variables used in the iteration
#of LOOCV for the model (nextoe)
remove(i)
remove(j)

varLIST=matrix(NA,nrow=sample_size,ncol=6)

storedDATA=vector("list",sample_size)

#this for loop generates the matrix of predictor variables from
#each iteration of the LOOCV

for(i in 1:sample_size) {
#loads in the final narrowband indices selected by each iteration of LOOCV
    load(paste("nextoe",i,".Rdata",sep=""))
  
  NextForm=nextoe[-which(sapply(nextoe,is.null)==T)]
  splat=unlist(strsplit(NextForm[[length(NextForm)]],"+",fixed=T))[-1]
  
  varLIST[i,]=splat
  
  LEN=length(NextForm)-1
  
  banded=matrix(NA,nrow=sample_size,ncol=LEN)
  
  b1=rep(NA,times=LEN)
  b2=rep(NA,times=LEN)
  
  for(j in 1:LEN){
    b1[j]=as.numeric(unlist(strsplit(splat[j],"_",fixed=T))[2])
    b2[j]=as.numeric(unlist(strsplit(splat[j],"_",fixed=T))[3])
    banded[,j]=(bandsAGG[,b2[j]]-bandsAGG[,b1[j]])/(bandsAGG[,b2[j]]+bandsAGG[,b1[j]])
  }
  

  storedDATA[[i]]=data.frame(groAGG$gro,banded,groAGG$elevation)
  colnames(storedDATA[[i]])=paste("X",c(1:ncol(storedDATA[[i]])),sep="")
  colnames(storedDATA[[i]])[1]="gro"
}



predicted=rep(0,times=sample_size)
#predictions for each iteration of LOOCV
for(o in 1:sample_size) {
  dats.use1=storedDATA[[o]]
  dat.train=dats.use1[-o,]
  dat.test=dats.use1[o,]
  
  m1=lm(gro~.,data=dat.train)
  
  predicted[o]=predict(m1,newdata=dat.test)
  
}


plot(predicted~groAGG$gro)

#functions to compute out-of-sample R2 and RMSE from model predictions
R2<-function(m,x,logo=F) {
  mors<-na.omit(data.frame(m,x))
  m<-mors$m
  x<-mors$x
  
  if(logo==F) {return(1-RSS(m,x)/TSS(x))} else {return(1-RSS1(m,x)/TSS1(x))}
}

rmse<-function(m,x)
{
  sqrt(mean((x-m)^2))
}

#out-of-sample predictive accuracy:
R2(predicted,gro) #r-squared
rmse(predicted,gro) #root mean squared error (RMSE)

