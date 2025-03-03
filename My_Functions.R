####################################
myPercentile<-function(x, prob){
  is.array(x)
  is.numeric(prob)
  x.sort<-sort(x)
  position<-length(x.sort)*prob
  if(is.integer(position)==TRUE){
    value<-(x.sort[position]+x.sort[position+1])/2
  }else if(is.integer(position)==FALSE){
  position1<-floor(position)
  value<-x.sort[position1]}
  return(value)
}

#############################################
myECDF<-function(x, N){
  is.array(x)
  is.numeric(N)
  if(N>max(x)|N<min(x)){
    paste('Please do not do something silly')
  }else{
  x.sort<-sort(x)
  below<-length(x.sort[x.sort<=N])
  percentile<-below/length(x.sort)
  return(percentile)
  }
}
######################################
SEmean<-function(x){
  
  is.array(x)
  n<-length(x)
  
  SEmean<-sqrt(sum((x-mean(x))^2)*(1/(n*(n-1))))
  
  return(SEmean)
}
####################################
Boot.SE<-function(x, FUN='mean'){
  
  is.array(x)
  n<-length(x)
  if(FUN=='mean'){
    SEmean<-sqrt(sum((x-mean(x))^2)*(1/(n-1)))
    
    return(SEmean)
  }else if(FUN=='median'){
    SEmedian<-sqrt(sum((x-median(x))^2)*(1/(n-1)))
    
    return(SEmedian)
    
  }
}
##############################################

MyBootstrap<-function(x, time, FUN='mean'){
  is.array(x)
  n<-length(x)
  BootFUN<-c()
  
  if(FUN=='mean'){
    for(i in 1:time-1){
      rlist<-ceiling(runif(n, min=0, max=n))
      Bootlist<-x[rlist]
      BootFUN[i]<-mean(Bootlist)
    }
    BootFUN[time]<-mean(x)
    return(mean(BootFUN))
  }else if(FUN=='median'){
    for(i in 1:time-1){
      rlist<-ceiling(runif(n, min=0, max=n))
      Bootlist<-x[rlist]
      BootFUN[i]<-mean(Bootlist)
    }
    BootFUN[time]<-median(x)
    return(median(BootFUN))
  } else if(FUN=='SE(mean)'){
    for(i in 1:time-1){
      rlist<-ceiling(runif(n, min=0, max=n))
      Bootlist<-x[rlist]
      BootFUN[i]<-mean(Bootlist)
    }
    BootFUN[time]<-mean(x)
    return(Boot.SE(BootFUN, FUN='mean'))
  } else if(FUN=='SE(median)'){
    
    for(i in 1:time-1){
      rlist<-ceiling(runif(n, min=0, max=n))
      Bootlist<-x[rlist]
      BootFUN[i]<-median(Bootlist)
    }
    BootFUN[time]<-median(x)
    return(Boot.SE(BootFUN, FUN='mean'))
  }
}
###########################################

MyBootstrapArray<-function(x, time, FUN='mean'){
  is.array(x)
  n<-length(x)
  BootFUN<-c()
  
  if(FUN=='mean'){
    for(i in 1:time-1){ 
      #bootstrap needs to include the original data
      rlist<-ceiling(runif(n, min=0, max=n))
      Bootlist<-x[rlist]
      BootFUN[i]<-mean(Bootlist)
    }
    BootFUN[time]<-mean(x)
    return(BootFUN)
  }else if(FUN=='median'){
    for(i in 1:time-1){
      rlist<-ceiling(runif(n, min=0, max=n))
      Bootlist<-x[rlist]
      BootFUN[i]<-median(Bootlist)
    }
    BootFUN[time]<-median(x)
    return(BootFUN)
  }
}

######################################################

MylmCoefficient<-function(x,y){
  is.array(x)
  is.array(y)
  beta1_up<-sum((x-mean(x))*(y-mean(y)))
  beta1_down<-sum((x-mean(x))^2)
  beta1<-beta1_up/beta1_down
  beta0<-mean(y)-beta1*mean(x)
  
  return(c(beta1, beta0))
}
#################################################

MyJkArray<-function(x, FUN='mean'){
  is.array(x)
  JKlist<-c()
  if(FUN=='mean'){
    for(i in 1:length(x)){
      JKlist[i]<-mean(x[-i])
    }
    return(JKlist)
  }else if(FUN=='median'){
    for(i in 1:length(x)){
      JKlist[i]<-median(x[-i])
    }
    return(JKlist)
  }
}
###############################################
JK.SE<-function(x, FUN='mean'){
  is.array(x)
  if(FUN=='mean'){
    xlist<-MyJkArray(x, FUN='mean')
  }else if(FUN=='median'){
    xlist<-MyJkArray(x, FUN='median')
  }
  n<-length(x)
  RightPart<-(n-1)/n
  LeftPart<-sum((xlist-mean(xlist))^2)
  JKSE<-sqrt(RightPart*LeftPart)
  return(JKSE)
}
#########################################
JK.SE.forArray<-function(JKarray){
  is.array(JKarray)
  n<-length(JKarray)
  RightPart<-(n-1)/n
  LeftPart<-sum((JKarray-mean(JKarray))^2)
  JKSE<-sqrt(RightPart*LeftPart)
  return(JKSE)
}
#########################################
########Correlation p value###################
MyCorrP<-function(x, y, n=5000){
  yPermute<-myPermutation(y, n)
  CoeffRandom<-c()
  for(i in 1:5000){
    CoeffRandom[i]<-MylmCoefficient(x, as.vector(as.matrix(yPermute[i])))[1]
  }
  Beta1<-MylmCoefficient(x,y)[1]
  if(Beta1>max(CoeffRandom)|Beta1<min(CoeffRandom)){
    return(0)
  }else if(Beta1>median(CoeffRandom)){
    Pvalue<-1-myECDF(CoeffRandom, Beta1)
    return(Pvalue)
  }else if(Beta1<=median(CoeffRandom)){
    Pvalue<-myECDF(CoeffRandom, Beta1)
    return(Pvalue)
  }
  #Pvalue<-1-myECDF(CoeffRandom, Beta1)
  #return(Pvalue)
}
##################################################
############P value of difference between two groups###################
DiffMeanP<-function(x, y, n=5000){
  PermuteLs<-myPermutPairs(x, y, n)
  xPermute<-PermuteLs[[1]]
  yPermute<-PermuteLs[[2]]
  xMean<-Mysapply(xPermute, 'mean')
  yMean<-Mysapply(yPermute, 'mean')
  xyDiff<-abs(xMean-yMean)
  xyDiff.o<-mean(x)-mean(y)
  #if(xyDiff.o>max(xyDiff)|xyDiff.o<min(xyDiff)){
  #return(0)
  #}else if(xyDiff.o>median(xyDiff)){
  #return(1-myECDF(xyDiff, xyDiff.o))
  #}else if(xyDiff.o<median(xyDiff)){
  #return(myECDF(xyDiff, xyDiff.o))
  #}else if(xyDiff.o==median(xyDiff)){
  #return(1)
  #}
  #myECDF(xyDiff, xyDiff.o)
  #1-myECDF(xyDiff, abs(xyDiff.o))
  if(abs(xyDiff.o)>max(xyDiff)){
    return(0)
  }else return(1-myECDF(xyDiff, abs(xyDiff.o)))
}


###############################################
####Random sampling with replacement###########
mySampling<-function(n){
  R<-c()
  StudentLeft<-list()
  StudentLeft[[1]]<-c(1:n)
  StudentSelect<-c()
  for(i in 1:n){
    R[i]<-ceiling(runif(1, 0, length(StudentLeft[[i]])))  
    StudentLeft[[i+1]]<-StudentLeft[[i]][-R[i]]
    StudentSelect[i]<-StudentLeft[[i]][R[i]]
  }
  StudentSelect[[n]]<-StudentLeft[[n]]
  return(StudentSelect)
}
######Another function to sample with replacement########
mySampling2<-function(n){
  sample<-rank(runif(n, min=0, max=n))
  return(sample)
}
#############################################
#####Permutation for single vector#########
myPermutation<-function(x, n=5000){
  is.array(x)
  permute<-as.data.frame(matrix(nrow=length(x), ncol=n))
  for(i in 1:n){
    sequence<-mySampling(length(x))
    permute[i]<-x[sequence]
  }
  return(permute)
}
#########Permutation for paired dataset##########
myPermutPairs<-function(x, y, n=5000){
  is.array(x)
  is.array(y)
  bind<-c(x,y)
  permute<-as.data.frame(matrix(nrow=length(bind), ncol=n))
  
  for(i in 1:n){
    sequence<-mySampling(length(bind))
    permute[i]<-bind[sequence]
  }
  PermuteX<-permute[1:length(x),]
  PermuteY<-permute[(length(x)+1):length(bind), ]
  # permute.xy<-as.data.frame(matrix(nrow=length(x),
  #ncol=2*n))
  #for(i in 1:n){
  #  permute.xy[2*i-1]<-PermuteX[i]
  # permute.xy[2*i]<-PermuteY[i]
  
  #}
  
  Result<-list(X=PermuteX, Y=PermuteY)
  return(Result)
}
#########Replacing sapply##############
Mysapply<-function(x, FUN='mean'){
  is.data.frame(x)
  k<-c()
  if(FUN=='mean'){
    for(i in 1:ncol(x)){
      k[i]<-mean(x[,i])
    }
    return(k)
  }else if(FUN=='median'){
    for(i in 1:ncol(x)){
      k[i]<-median(x[,i])
    }
    return(k)
  }else if(FUN=='sum'){
    for(i in 1:ncol(x)){
      k[i]<-sum(x[,i])
    }
    return(k)
  }else if(FUN=='max'){
    for(i in 1:ncol(x)){
      k[i]<-max(x[,i])
    }
    return(k)
  }else if(FUN=='min'){
    for(i in 1:ncol(x)){
      k[i]<-min(x[,i])
    }
    return(k)
  }else if(FUN=='sd'){
    for(i in 1:ncol(x)){
      k[i]<-sd(x[,i])
    }
    return(k)
  }else if(FUN=='quantile'){
    Q<-as.data.frame(matrix(nrow=4, ncol=ncol(x)))
    for(i in 1:ncol(x)){
      Q[i]<-myPercentile(as.vector(x[,i]),
                         c( 0.25, 0.5, 0.75, 1))
    }
    row.names(Q)<-c('0.25', '0.5',
                    '0.75', '1')
    return(Q)
  }
}
Mysapply(iris[1:4], FUN='median')
#################################################
myMRPP<-function(data, dist_method, hcut_method, k, n=1000 ){
  data_dist<-vegdist(data, method=dist_method)
  tree<-hclust(data_dist, method=hcut_method)
  groups<-as.factor(cutree(tree, k=k))
  delta<-c()
  for(j in 1:n){
    seq<-sample(c(1:nrow(data)), nrow(data), replace = FALSE)
    data.permute<-data[seq,]
    delta_origin<-c()
    for(i in 1:k){
      delta_origin[i]<-mean(vegdist(data.permute[which(groups==1),], 
                                    method=dist_method))
    }
    delta[j]<-sum(table(groups)/length(groups)*delta_origin)
  }
  delta.1<-c()
  for(m in 1:k){
    delta.1[m]<-mean(vegdist(data[which(groups==m),], 
                             method=dist_method))
  }
  delta.2<-sum(table(groups)/length(groups)*delta.1)
  significance<-ecdf(delta)(delta.2)
  return(list(group.Delta=delta.1, obs.Delta=delta.2, N=factor(groups), 
              Prob=significance, 
              Delta_Permute=delta,
              expected_Delta=mean(delta)))
}
