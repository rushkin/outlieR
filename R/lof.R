#' Local Outlier Factor (LOF) values
#'
#' Calculates the local outlier factor (LOF) for multidimensional data
#'
#' @param x a matrix or data frame of data points. Each row is a data point, each column is a dimension.
#' @param k number of nearest neighbors in the LOF calculation (default 5).
#' @param metric distance metric to use. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param q the power of the Minkowski distance.
#' @return a vector of LOF values, in the same order as the data points in \code{x}. NA values will be returned for those data points which contained at least one NA coordinate.
#' @export
#'
#' @examples
lof=function(x, k=5, metric='euclidean', q=3){

  x=as.matrix(x)

  indices=1:nrow(x)
  indna=which(apply(x,1,function(x){any(is.na(x))}))
  N=length(indices)

  distmat=as.matrix(dist(x,method=metric,p=q))
  diag(distmat)=NA


  #Indices of k nearest neighbors. A point may have more than k values in case there are distance-ties. So the nrows of this matrix is from the longest series of indices, and missing values are filled with NAs.
  nk=apply(distmat,1,function(y){

    ind=which(rank(y,ties.method = 'min',na.last=TRUE)<=k)

    ind=c(ind,rep(NA,N-length(ind)))
    return(ind)

  })
  nk=nk[rowSums(nk,na.rm=TRUE)>0,]

  #Add to it another row, which simply contains the indices. This is convenient later. The nrows will become l.
  l=nrow(nk)+1
  nkaugmented=rbind(nk,indices)

  #k-distances of points. returns a vector (named)
  kdist=apply(nkaugmented,2,function(y){

    temp=distmat[y[l],y[-l]]
    temp=temp[!is.na(temp)]
    if(length(temp)==0) return(NA)

    return(max(temp))

  })

  #convert kdist to a matrix by replicating the same row:
  kdist=matrix(kdist,nrow = 1)
  kdist=kdist[rep(1,N),]

  #Calculate reachability distance matrix
  reachdist=pmax(kdist,distmat)

  #Local reachability density. Returns a vector (named)
  lrd=1/apply(nkaugmented,2,function(y){
    temp=reachdist[y[l],y[-l]]
    temp=temp[!is.na(temp)]
    if(length(temp)==0) return(NA)
    return(mean(temp))
  })

  #Local outlier factor. Returns a vector
  lof=apply(nk,2,function(y){
    mean(lrd[y],na.rm=TRUE)
  })/lrd
  names(lof)=NULL

  lof[indna]=NA
  return(lof)
}
