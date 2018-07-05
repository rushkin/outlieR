abserror=function(par,points,metric='euc',q=2){

  a=switch(
    metric
    ,'euc'=sum(sqrt(apply(points,1,function(y){sum((y-par)^2)})))
    ,'max'=sum(apply(points,1,function(y){max(abs(y-par))}))
    ,'man'=sum(apply(points,1,function(y){sum(abs(y-par))}))
    ,'can'=sum(apply(points,1,function(y){sum(abs(y-par)/(abs(y)+abs(par)))}))
    ,'min'=sum((apply(points,1,function(y){sum(abs(y-par)^q)}))^(1/q))
    ,'bin'=sum(apply(points,1,function(y){
      yy=(y!=0)
      pp=(par!=0)
      return(length(which(xor(yy,pp)))/length(which(yy|pp)))
      }))
  )

  return(a)
}

#' Median of multidimensional data
#'
#' @param x a matrix, data frame or vector of data points (a vector will be understood as 1D data, equivalent to a 1-column matrix). Rows with NA values will be dropped.
#' @param metric a choice of 'euclidean','maximum','manhattan','canberra','minkowski', or 'binary'. Can be abbreviated to first three letters, case-insensitive.
#' @param q the power in the minkowski metric. If \code{q=2}, equivalent to 'euclidean'.
#' @param simple1d boolean, if TRUE, will resort to the standard \code{median} function in case the data is 1D. This means that only euclidean metric will be used. However, most metrics are identical to euclidean in 1D anyway, so it is usually no big loss.
#' @param ... other arguments passed to \code{Rcgmin::Rcgmin}.
#'
#'@details The median is calculated as the minimizer m of the sum of distances ||x_i-m||, where x_i are the data points.
#'In an N-dimensional space, each x_i is an N-dimensional vector, and so is m.
#'The distances are calculated according to the metric chosen (Euclidean by default). Minimization is done with \code{Rcgmin} (note that exotic metrics you might cause convergence issues).
#'In 1D case, unless \code{simple1d} is FALSE, the standard \code{median} function is used instead of this method.
#'
#'
#' @return a vector of median components.
#' @export
#'
#' @examples
mmedian=function(x,metric='euclidean',q=3,simple1d=TRUE,...){

    metric=tolower(substr(metric,1,3))

    if(is.null(dim(x))){
      x=data.frame(x=x)
    }

    #coerce to a data frame
    x=as.data.frame(x)

    #Special case: 1D
    if((ncol(x)==1)&simple1d){
      return(median(x[,1],na.rm=TRUE))
    }




    #General case:
      indna=apply(x,1,function(x){any(is.na(x))})
      x1=x[!indna,]
      if(nrow(x1)==0) return(rep(NA,ncol(x)))

    par=colSums(x1)/nrow(x1)
    suppressWarnings({ae=Rcgmin::Rcgmin(par=par,fn=abserror,points=x1,metric=metric,q=q)})
    if(ae$convergence>0){
      warning('Method did not converge!')
    }

    names(ae$par)=NULL

    return(ae$par)



  }
