#' Flag outliers
#'
#' Flag outliers in a multidimensional data set
#'
#' @param x a matrix, data frame or vector of data points (a vector will be understood as 1D data, equivalent to a 1-column matrix). Each row is a data point and each column is a dimension. NA values are allowed and will produce NAs in the output.
#' @param level threshold for finding outliers. Meant to be from 0 to 1. Smaller values mean a higher bar for outliers and so typically detect fewer outliers.
#' For probabilistic methods (such as based on normal distribution), it is the significance level. For LOF, we set the LOF threshold to \code{1/level -1}.
#' \code{level=0} will flag no outliers. LOF method with \code{level=1} will flag all points as outliers (or as many as \code{nmax} allows).
#' Note however, that Grubbs method may still leave some points unflagged even with \code{level=1}.
#' @param nmax the maximum number of outliers to remove. If NULL, ignored.
#' @param side if set to 'left', 'right' or 'both' (can be abbreviated to one letter and case-insensitive) will flag only the outliers on the left, right or both ends of the 1D distribution. If NULL, all outliers will be flagged. If the data is not 1D, \code{side} will be ignored. Note that for the methods that only find outliers on the sides of the distribution (e.g Chauvenet) NULL and 'both' give equivalent results.
#' @param crit criterion to use for identifying outliers. Currently, can be either 'LOF' or 'Grubbs'. Any unambiguous substring can be given, case insensitive. If 'Grubbs',
#' the 1D Grubbs method will be applied along each principal axis of the data and points deemed outliers along at least one axis will be flagged.
#' @param asInt if TRUE, the flag values will be integers (1 for outlier and 0 otherwise). If FALSE, boolean
#' @param k number of nearest neighbors for the LOF calculation
#' @param metric distance metric to use. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given, case insensitive.
#' @param q the power of the Minkowski distance.
#'
#' @return a boolean or integer (depending on \code{asInt}) vector of the same length as the number of points in the data, containing 1 (TRUE) if a data point is an outlier, 0 (FALSE) if it is not and NA if a point in the data contained NA value(s), the
#' @export

flag=function(x, level=0.1, nmax=NULL, side=NULL, crit='lof', asInt=TRUE, k=5, metric='euclidean', q=3){
  criteria=c('lof','grubbs')
  crit=pmatch(tolower(crit),criteria)
  if(is.na(crit)) stop('invalid outlier criterion name')
  if(crit==-1) stop('ambiguous outlier criterion name')

  if(level<0) stop('negative level is invalid')

  if(!is.null(side)){
    side=tolower(substr(side,1,1))
    if(!(side %in% c('l','r'))) side='b'
  }

  x=as.matrix(x)


  if(is.null(nmax)){
    nmax=nrow(x)
  }
  nmax=min(nmax,nrow(x))




  if(crit==1){

    thresh=(1-level)/level

    scores=lof(x,k=k,metric=metric, q=q)

    #Implement sidedness
    if(ncol(x)==1){
      if(!is.null(side)) {

        if(side %in% c('l','r','b')){
          s=seq_along(scores)
          temp=cbind(cbind(x,scores),s)
          temp=temp[order(temp[,1]),]
          suppressWarnings({goodrange=range(which(temp[,2]>thresh))})

          if(side=='b'){
            temp[(s>goodrange[1])&(s<goodrange[2]),2]=0
          }
          if(side=='l'){
            temp[(s>goodrange[1]),2]=0
          }
          if(side=='r'){
            temp[(s<goodrange[2]),2]=0
          }

          scores=temp[order(temp[,3]),2]

        }
      }
    }

    n=length(which(!is.na(scores)))

    nmax=min(n,nmax)

    ind=order(scores,decreasing = TRUE)[1:nmax]
    ind=ind[scores[ind]>thresh]
  }#End of LOF

  if(crit==2){

    #In this method only the outliers on the sides are possible, so null side is treated simply as the two-sided test
    #If the data is multidimensional, one-sided tests are not applied
    if(is.null(side)|(ncol(x)>1)) side='b'

    ind=c()
    for( i in 1:nmax){
      y=zscore(x,side=side)

      i_to_drop=which(apply(apply(y,2,grubbs1D,level=level,side=side),1,function(x){any(x,na.rm=TRUE)}))
      if(length(i_to_drop)==0) break

      #Add the found index to the list and then populate that row of x with NA, to exclude it from next iteration.
      ind=c(ind,i_to_drop)
      x[i_to_drop,]=NA
    }

  }#End of Grubbs








  flag=rep(FALSE,nrow(x))
  flag[ind]=TRUE


  if(asInt){
    flag=as.integer(flag)
  }

  return(flag)


}




#' Impute outliers
#'
#'Impute detected outliers in a multidimensional data set
#'
#' @param x a matrix, data frame or vector of data points (a vector will be understood as 1D data, equivalent to a 1-column matrix). Each row is a data point and each column is a dimension. NA values are allowed and will produce NAs in the output.
#' @param flag a boolean or integer (0-or-1) vector flagging outliers, such as produced by the function \code{flag}. If NULL, further arguments will be used to compute it here by calling \code{flag}.
#' @param fill method for imputing (or removing) outliers. If numeric or NA, it is the value that will replace the outliers. It the data is K-dimensional, \code{fill} is expected to be a vector of length K. If longer, the first K components will be used, and if shorter, the vector will be extended by NAs.
#' Alternatively, \code{fill} can be a character string. Values 'mean' and 'median' replace outliers with the mean or (multidimensional) median of the rest of the remaining data, 'random' generates random replacement values drawn from the estimated probability distribution of the non-outlier data-points,
#' 'remove' removes the outliers by calling the function \code{purge}. Any unambiguous substring can be given, case insensitive.
#' @param level passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param nmax passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param side passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param crit passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param k passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param metric distance metric to be used in LOF (if \code{flag} is not provided) as well as for multidimensional median if \code{fill} is 'median'. A choice of 'euclidean','maximum','manhattan','canberra','minkowski', or 'binary'. Any unambiguous substring can be given, case insensitive.
#' @param q power in Minkowski metric, used if \code{fill='median'} and \code{metric='minkowski'}
#' @param ... passed to \code{Rcgmin} if the argument \code{fill} is 'median' and data is multidimensional
#'
#' @return object like \code{x} but with outliers imputed.
#' @details The output object will be a vector, a matrix or a data-frame, depending on what \code{x} was. Row names, column names or (if \code{x} was a named vector) names will be kept.
#' @export
#'
#' @examples
impute=function(x, flag=NULL, fill='mean', level=0.1, nmax=NULL, side=NULL, crit='lof', k=5, metric='euclidean', q=3, ...){

  fills=c('mean','median','random','remove')

  if(!is.character(fill)){
    fillval=fill
  }else{
    fillval=NULL

    fill=pmatch(tolower(fill),fills)
    if(is.na(fill)) stop('invalid fill option')
    if(fill==-1) stop('ambiguous fill option')
    fill=fills[fill]

    if(fill=='remove'){
      return(purge(x=x,flag=flag,level=level,nmax=nmax,side=side,crit=crit,k=k,metric=metric,q=q))
    }
  }

  if(is.matrix(x)){
    orig.format='mat'
    orig.rownames=rownames(x)
    orig.colnames=colnames(x)
  }

  if(is.data.frame(x)){
    orig.format='df'
    orig.rownames=rownames(x)
    orig.colnames=colnames(x)
  }
  if(is.null(dim(x))){
    orig.format='vec'
    orig.names=names(x)
  }

  x=as.matrix(x)


  if(is.null(flag)){
    flag=flag(x=x,level=level, nmax=nmax, side=side,crit=crit,k=k, metric=metric, q=q, asInt=FALSE)
  }

  flag=as.logical(flag)
  flag[is.na(flag)]=FALSE
  n_to_impute=length(which(flag))
  xvalid=x[!flag,,drop=FALSE]

  if(is.null(fillval)){
    fillval=switch(fill
                 ,'mean'={
                   temp=matrix(apply(xvalid,2,function(x){mean(x,na.rm=TRUE)}),ncol=ncol(x))
                   temp[rep(1,n_to_impute),]
                 }
                 ,'median'={
                   temp=matrix(mmedian(xvalid,metric=metric,q=q,...),ncol=ncol(x))
                   temp[rep(1,n_to_impute),]
                 }
                 ,'random'={
                   err=sapply(bw(xvalid),function(bw){rnorm(n_to_impute,0,bw)})
                   if(n_to_impute==1) err=matrix(err,nrow = n_to_impute)
                   ind=sample(which(apply(xvalid,1,function(x){all(!is.na(x))})),n_to_impute,replace=TRUE)
                   temp=xvalid[ind,]+err
                 }
                 )
  }else{
    fillval=fillval[1:ncol(x)]
    fillval=matrix(fillval,ncol=ncol(x))
    fillval=fillval[rep(1,n_to_impute),]
  }



  x[flag,]=fillval


  x=switch(orig.format
           ,'vec'={y=x[,1]; names(y)=orig.names; y}
           ,'mat'={y=as.matrix(x); colnames(y)=orig.colnames; rownames(y)=orig.rownames; y}
           ,'df'={y=as.data.frame(x); rownames(y)=orig.rownames; colnames(y)=orig.colnames; y}
           )

  return(x)
}

#' Remove outliers
#'
#'Remove detected outliers in a multidimensional data set
#'
#' @param x a matrix, data frame or vector of data points (a vector will be understood as 1D data, equivalent to a 1-column matrix). Each row is a data point and each column is a dimension. NA values are allowed and will produce NAs in the output.
#' @param flag a boolean or integer (0-or-1) vector flagging outliers, such as produced by the function \code{flag}. If NULL, further arguments will be used to compute it here by calling \code{flag}.
#' @param level passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param nmax passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param side passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param crit passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param k passed to the function \code{flag} if the argument 'flag' is NULL or missing
#' @param metric distance metric to be used in LOF (if \code{flag} is not provided). A choice of 'euclidean','maximum','manhattan','canberra','minkowski', or 'binary'. Any unambiguous substring can be given, case insensitive.
#' @param q power in Minkowski metric, used if \code{fill='median'} and \code{metric='minkowski'}
#'
#' @return object like \code{x} but with outliers removed.
#' @details The output object will be a vector, a matrix or a data-frame, depending on what \code{x} was. Row names, column names or (if \code{x} was a named vector) names will be kept.
#' @export
#'
#' @examples
purge=function(x, flag=NULL, level=0.1, nmax=NULL, side=NULL, crit='lof', k=5, metric='euclidean', q=3){

  if(is.null(flag)){
    flag=flag(x=x,level=level, nmax=nmax, side=side,crit=crit,k=k, metric=metric, q=q, asInt=FALSE)
  }

  flag=as.logical(flag)
  flag[is.na(flag)]=FALSE

  if(is.null(dim(x))){
    x=x[!flag]
  }else{
      x=x[!flag,,drop=FALSE]
    }

  return(x)
}




#' Density kernel bandwidth
#'
#' The Silverman's rule of thumb kernel bandwidth for Gaussian 2nd order density kernel, for multidimensional data.
#'
#' @param x data, as a matrix. Each row is a data points, each column is a dimension. NA values will be ignored.
#'
#' @return a vector of bandwidths, the length is the number of columns in x.
#'
#' @examples
bw=function(x){
  s=apply(x,2,sd,na.rm=TRUE)
  mad=apply(x,2,function(x){mean(abs(x-mean(x,na.rm=TRUE)),na.rm=TRUE)})
  iq=apply(x,2,IQR,na.rm=TRUE)
  n=apply(x,2,function(x){length(which(!is.na(x)))})
  bw=1.06*pmin(s,0.6745*mad,0.7413*iq,na.rm=TRUE)*(n^(-1/(4+length(n))))

  return(bw)

}
