####### cross validation in PENG and PAUL#####


CV<-function(B,phi.aux,sig,Lambda,data.list,n){
  
  ##partI: -2loglikelihood
  part.like<-loglike.all(B,phi.aux,sig,Lambda,data.list,n)
  
  ##Part for B
  part.B<-CV.B(B,phi.aux,sig,Lambda,n)
  
  ##part for (log(Lambda),log(sig^2))
  part.Lam<-CV.Lambda(B,phi.aux,sig,Lambda,data.list,n)
  
  ##
  result<-part.like+part.B+part.Lam
  return(result)
}

loglike<-function(B,phi.aux,sig,Lambda,data.list,i){
  M<-nrow(B)
  r<-ncol(B)
  t.B<-t(B)
  
  temp<-phi.aux[[i]]
  phi<-temp[[1]]
  psi<-temp[[2]]
  Q<-sig^2*diag(1/Lambda)+t.B%*%psi%*%B
  Q.inv<-solve(Q)
  m<-ncol(phi)
  
  data.temp<-data.list[[i]]
  y<-data.temp[[1]][,1]
  
  P.inv<-(diag(1,m)-t(phi)%*%B%*%Q.inv%*%t.B%*%phi)/sig^2
  log.det.P<-2*(m-r)*log(sig)+sum(log(Lambda))+log(det(Q))
  
  result<-sum(diag(P.inv%*%y%*%t(y)))+log.det.P+m*log(2*pi)
  return(result) 
}

### the mean -2*loglikelihood over all observations
loglike.all<-function(B,phi.aux,sig,Lambda,data.list,n){
  temp<-apply(matrix(1:n), MARGIN=1, loglike, B=B,phi.aux=phi.aux, sig=sig, Lambda=Lambda, data.list=data.list) 
  result<-mean(temp)
  return(result)
}


CV.B<-function(B,phi.aux,sig,Lambda,n){
  ###para: B, sig, Lambda: current parameters
  ##phi.aux: data; n--sample size
  ##return CV score: the stiefel part; averaged over all obs
  
  M<-nrow(B)
  r<-ncol(B)
  t.B<-t(B)
  index.t<-Perm.right.index(M,r)
  
  ##(i) H^{-1}: 
  temp<-apply(matrix(1:n),MARGIN=1,FB.HB,B=B,phi.aux=phi.aux,sig=sig,Lambda=Lambda,index=index.t)
  
  temp.all<-apply(temp,MARGIN=1,sum)
  temp.all<-matrix(temp.all,M*r, 1+M*r)
  FB.all<-matrix(temp.all[,1],M,r)
  gradF.all<-FB.all-B%*%t(FB.all)%*%B
  H.all<-temp.all[,-1]
  temp1<-kronecker(t.B,B)[,index.t]
  FBB.all<-H.all-temp1%*%H.all
  
  I<-(-0.5)*B%*%t(FB.all)
  J<-(-0.5)*t(FB.all)%*%B
  PI<-diag(1,M)-B%*%t.B  
  eye.M<-diag(1,M)
  PI.half<-eye.M-0.5*B%*%t.B   
  
  eye.r<-diag(1,r)
  temp1<-kronecker(eye.r,I)+kronecker(t(J),eye.M)+kronecker(J,PI)
  temp2<-0.5*kronecker(t(FB.all),B)+0.5*kronecker(t.B,FB.all)
  temp2<-temp2[,index.t]
  
  L1<-FBB.all+temp1+temp2
  L2<-kronecker(eye.r,t.B)+kronecker(t.B,eye.r)[,index.t]
  
  A<-rbind(L1,L2)
  svd.A<-svd(A)
  U.A<-svd.A$u
  V.A<-svd.A$v
  D.A<-svd.A$d
  temp.A<-V.A%*%diag(1/D.A)%*%t(U.A)   #H^{-1}
  
  ##(ii):
  partII<-numeric(n)               ##g_c(gradF_i,H^{-1}(gradF_i))
  partIII<-numeric(n)              ##H_i(H^{-1}(gradF_i),H^{-1}(gradF_i))
  
  for (i in 1:n){
    cur<-temp[,i]
    FB.c<-matrix(cur[1:(M*r)],M,r)
    HB.c<-matrix(cur[-(1:(M*r))],M*r,M*r)
    gradF.c<-FB.c-B%*%t(FB.c)%*%B           
    
    ##delta_i=H^{-1}(gradF_i)
    C.c<-c(as.vector(gradF.c),numeric(r^2))  
    delta.c<-temp.A%*%C.c
    delta.c<-matrix(delta.c,M,r)  
    
    ##FBB_i(H^{-1}(gradF_i),H^{-1}(gradF_i))
    FBB.temp<-matrix(HB.c%*%as.vector(delta.c),M,r)%*%t(delta.c)
    FBB.c<-sum(diag(FBB.temp))  
    
    ##partIII: H_i(delta_i,delta_i)
    temp1<-(t(FB.c)%*%delta.c%*%t.B+t.B%*%delta.c%*%t(FB.c))%*%delta.c
    temp2<-(t.B%*%FB.c+t(FB.c)%*%B)%*%t(delta.c)%*%PI%*%delta.c
    partIII[i]<-1.5*(FBB.c+0.5*sum(diag(temp1))-0.5*sum(diag(temp2)))   
    
    ##partII: g_c(gradF_i,delta.c)
    temp.2<-t(gradF.c)%*%PI.half%*%delta.c
    partII[i]<-sum(diag(temp.2))
  }
  
  result<-mean(partII+partIII)
  return(result)
}

####
CV.Lambda<-function(B,phi.aux,sig,Lambda,data.list,n){
  ###para: B, sig, Lambda: current parameters
  ##phi.aux: data; n--sample size
  ##return:the CV part on for c(log(Lambda),log(sig^2)); averaged over all obs
  
  r<-ncol(B)   
  temp<-apply(matrix(1:n),MARGIN=1,F.Lambda,B=B,phi.aux=phi.aux,sig=sig,Lambda=Lambda, data.list=data.list)
  temp.all<-apply(temp,1,sum)
  temp.all<-matrix(temp.all,r+1,r+2) 
  Hess.all<-temp.all[,-1] 
  
  eigen.H<-eigen(Hess.all,symmetric=TRUE)
  U.H<-eigen.H$vectors
  D.H<-eigen.H$values
  Hess.inv<-U.H%*%diag(1/D.H)%*%t(U.H) 
  
  partII<-numeric(n)
  partIII<-numeric(n)
  
  for(i in 1:n){
    cur<-temp[,i]
    cur<-matrix(cur,r+1,r+2)
    gradF.c<-cur[,1]
    Hess.c<-cur[,-1]
    
    ##delta_i=H^{-1}(gradF_i)
    delta.c<-Hess.inv%*%gradF.c
    
    ##partII: <gradF_i, delta_i>
    partII[i]<-gradF.c%*%delta.c
    
    ##partIII: t(delta.c)%*%Hess.c%*%delta.c
    partIII[i]<-1.5*(t(delta.c)%*%Hess.c%*%delta.c)  
  }
  
  result<-mean(partII+partIII) 
  return(result)
}


