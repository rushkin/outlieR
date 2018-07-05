zscore=function(x,side='b'){

  x=apply(as.matrix(x),2,function(x){x-mean(x,na.rm=TRUE)})
  indvalid=apply(x,1,function(x){!any(is.na(x))})
  if(length(indvalid)==0) return(x)

  scores=x%*%prcomp(x[indvalid])$rotation

  scores=apply(scores,2,function(x){x/sd(x,na.rm=TRUE)})

  scores=switch(side
                ,'l'=-scores
                ,'r'=scores
                ,abs(scores)
  )

return(scores)

}

grubbs1D=function(y,level=0.5,side='b'){
  N=length(which(!is.na(y)))
  res=rep(FALSE,length(y))
  if(N<3) return(res)


  suspect=which.max(y)


  t2=switch(side
    ,'l'=(qt(p=level/(N),df=N-2))^2
    ,'r'=(qt(p=level/(N),df=N-2))^2
    ,(qt(p=level/(2*N),df=N-2))^2

  )

  Gcrit=(N-1)*sqrt(t2/(N*(N-2+t2)))

  res[suspect]=(y[suspect]>Gcrit)

  return(res)
}
