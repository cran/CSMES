#' Generate predictions for all Pareto-optimal ensemble classifier candidates selected through CSMES
#'
#' This function generates predictions for all pareto-optimal ensemble classifier candidates as identified through the first training stage of CSMES (\code{CSMES.ensSel}).
#'
#' @param ensSelModel ensemble selection model (output of \code{CSMES.ensSel})
#' @param newdata data.frame or matrix containing data to be scored
#' @import rpart zoo mco
#' @importFrom ROCR prediction performance
#' @importFrom caTools colAUC
#' @importFrom graphics axis lines matplot mtext plot points text
#' @importFrom stats approx as.formula dbeta
#' @return An object of the class \code{CSMES.predictPareto} which is a list with the following two components:
#' \item{Pareto_predictions_c}{A vector with class predictions.}
#' \item{Paret_predictions_p}{A vector with probability predictions.}
#' @export
#' @author Koen W. De Bock, \email{kdebock@@audencia.com}
#' @references De Bock, K.W., Lessmann, S. And Coussement, K., Cost-sensitive business failure prediction
#' when misclassification costs are uncertain: A heterogeneous ensemble selection approach,
#' European Journal of Operational Research (2020), doi: 10.1016/j.ejor.2020.01.052.
#' @seealso \code{\link{CSMES.ensSel}}, \code{\link{CSMES.predict}}, \code{\link{CSMES.ensNomCurve}}
#' @examples
#' ##load data
#' library(rpart)
#' library(zoo)
#' library(ROCR)
#' library(mco)
#' data(BFP)
#' ##generate random order vector
#' BFP_r<-BFP[sample(nrow(BFP),nrow(BFP)),]
#' size<-nrow(BFP_r)
#' ##size<-300
#' train<-BFP_r[1:floor(size/3),]
#' val<-BFP_r[ceiling(size/3):floor(2*size/3),]
#' test<-BFP_r[ceiling(2*size/3):size,]
#' ##generate a list containing model specifications for 100 CART decisions trees varying in the cp
#' ##and minsplit parameters, and trained on bootstrap samples (bagging)
#' rpartSpecs<-list()
#' for (i in 1:100){
#'   data<-train[sample(1:ncol(train),size=ncol(train),replace=TRUE),]
#'   str<-paste("rpartSpecs$rpart",i,"=rpart(as.formula(Class~.),data,method=\"class\",
#'   control=rpart.control(minsplit=",round(runif(1, min = 1, max = 20)),",cp=",runif(1,
#'   min = 0.05, max = 0.4),"))",sep="")
#'   eval(parse(text=str))
#' }
#' ##generate predictions for these models
#' hillclimb<-mat.or.vec(nrow(val),100)
#' for (i in 1:100){
#'   str<-paste("hillclimb[,",i,"]=predict(rpartSpecs[[i]],newdata=val)[,2]",sep="")
#'   eval(parse(text=str))
#' }
#' ##score the validation set used for ensemble selection, to be used for ensemble selection
#' ESmodel<-CSMES.ensSel(hillclimb,val$Class,obj1="FNR",obj2="FPR",selType="selection",
#' generations=10,popsize=12,plot=TRUE)
#' ## Create Ensemble nomination curve
#' enc<-CSMES.ensNomCurve(ESmodel,hillclimb,val$Class,curveType="costCurve",method="classPreds",
#' plot=FALSE)
CSMES.predictPareto<-function(ensSelModel,newdata) {
  popsize<-ensSelModel$popsize
  #extra check on the weights: if sum is zero (no selection whatsoever) the first member is chosen.
  #Otherwise the Pareto prediction function will produce NaN's for that ensemble.
  weights<-ensSelModel$weights
  weights[rowSums(weights)==0,1]<-1
  cutoff<-0.500000000001

  Pareto_predictions_p <- t(t(as.matrix(newdata[,1:(ncol(newdata)-1)])%*%t(weights))/rowSums(weights))

  Pareto_predictions_c <- array(0, c(nrow(newdata),popsize))
  for (j in 1:popsize) {
    Pareto_predictions_c[,j]=cbind((Pareto_predictions_p[,j]>cutoff)*1)
  }
  ans<- list(Pareto_predictions_c=Pareto_predictions_c,Pareto_predictions_p=Pareto_predictions_p)
  class(ans) <- "CSMES.predictPareto"
  ans
}
