#' CSMES scoring: generate predictions for the optimal ensemble classifier according to CSMES in function of cost information.
#'
#' This function generates predictions for a new data set (containing candidate member library predictions) using a CSMES model. Using Pareto-optimal ensemble definitions
#' generated through \code{CSMES.ensSel} and the ensemble nomination front generated using \code{CSMES.EnsNomCurve}, final ensemble predictions are generated in function of
#' cost information known to the user at the time of model scoring. The model allows for three scenarios: (1) the candidate ensemble is nominated in function of a specific cost
#' ratio, (2) the ensemble is nominated in function of partial AUCC (or a distribution over operating points) and (3) the candidate ensemble that is
#' optimal over the entire cost space in function of area under the cost or brier curve is chosen.
#'
#' @param ensNomCurve ensemble nomination curve object (output of \code{CSMES.ensNomCurve})
#' @param ensSelModel ensemble selection model (output of \code{CSMES.ensSel})
#' @param criterion This argument specifies which criterion determines the selection of the ensemble candidate that delivers predictions. Can be one of three options: "minEMC", "minAUCC" or "minPartAUCC".
#' @param costRatio Specifies the cost ratio used to determine expected misclassification cost. Only relvant when \code{criterion} is "minEMC".
#' @param partAUCC_mu Desired mean operating condition when \code{criterion} is "minPartAUCC" (partial area under the cost/brier curve).
#' @param partAUCC_sd Desired standard deviation when \code{criterion} is "minPartAUCC" (partial area under the cost/brier curve).
#' @param newdata matrix containing ensemble library member model predictions for new data set
#' @return An list with the following components:
#' \item{pred}{A matrix with model predictions. Both class and probability predictions are delivered.}
#' \item{criterion}{The criterion specified to determine the selection of the ensemble candidate.}
#' \item{costRatio}{The cost ratio in function of which the \code{criterion} "minEMC" has selected the optimal candidate ensemble that delivered predictions}
#' @import rpart zoo mco
#' @importFrom ROCR prediction performance
#' @importFrom caTools colAUC
#' @importFrom graphics axis lines matplot mtext plot points text
#' @importFrom stats approx as.formula dbeta
#' @export
#' @author Koen W. De Bock, \email{kdebock@@audencia.com}
#' @references De Bock, K.W., Lessmann, S. And Coussement, K., Cost-sensitive business failure prediction
#' when misclassification costs are uncertain: A heterogeneous ensemble selection approach,
#' European Journal of Operational Research (2020), doi: 10.1016/j.ejor.2020.01.052.
#' @seealso \code{\link{CSMES.ensSel}}, \code{\link{CSMES.predictPareto}}, \code{\link{CSMES.ensNomCurve}}
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
CSMES.predict<-function(ensSelModel, ensNomCurve, newdata, criterion=c("minEMC","minAUCC","minPartAUCC"),
                        costRatio=5,partAUCC_mu=0.5,partAUCC_sd=0.1) {

  estBetaParams <- function(mu, sd) {
    #kudos to assumednormal at StackExchange
    alpha <- ((1 - mu) / sd**2 - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }

  newdata_paretopreds<-CSMES.predictPareto(ensSelModel,newdata)
  nr_intervals=ensNomCurve$intervals
  P_plus<-ensNomCurve$incidence
  if (criterion=="minEMC") {
    alpha<-costRatio
    PC_plus=(P_plus*alpha)/(P_plus*alpha+(1-P_plus))
    #winner <- ensNomCurve$nomcurve[which(round(ensNomCurve$nomcurve[,1],digits=log10(nr_intervals))==round(PC_plus,digits=log10(nr_intervals))),3]
    winner <- ensNomCurve$nomcurve[which.min(abs(ensNomCurve$nomcurve[,1]-PC_plus)),3]
    pred <- newdata_paretopreds$Pareto_predictions_p[,winner]
  } else if (criterion=="minAUCC") {
    #identify member with overall minimal AUCC (default ensemble)
    areas_under_curve<-mat.or.vec(1,ensSelModel$popsize)
    x<-seq(from = 0, to = 1, by = 1/(ensNomCurve$intervals-1))
    for (i in 1:ensSelModel$popsize){
      areas_under_curve[i]<- sum(diff(x)*rollmean(ensNomCurve$curves[i,],2))
    }
    winner<- which.min(areas_under_curve)

    pred <- newdata_paretopreds$Pareto_predictions_p[,winner]
  } else if (criterion=="minPartAUCC") {
    #identify member with overall minimal partial AUCC for the cost ratio and sd provided
    areas_under_curve<-mat.or.vec(1,ensSelModel$popsize)
    x<-seq(from = 0, to = 1, by = 1/(ensNomCurve$intervals-1))
    id <- order(x)
    alpha<-costRatio
    PC_plus=(P_plus*alpha)/(P_plus*alpha+(1-P_plus))
    parests<-estBetaParams(PC_plus,partAUCC_sd)
    ydist<-dbeta(x,parests[[1]],parests[[2]])
    ydist[which(is.infinite(ydist))]<-max(ydist[which(!is.infinite(ydist))])*10
    for (i in 1:ensSelModel$popsize){


      areas_under_curve[i] <- sum(diff(x[id])*rollmean(ensNomCurve$curves[i,id],2)*rollmean(ydist[id],2))
    }
    winner<- which.min(areas_under_curve)
    pred <- newdata_paretopreds$Pareto_predictions_p[,winner]
  }


  formula <- as.formula(dependent~.)
  n <- length(newdata[,1])
  vardep <- newdata[,as.character(formula[[2]])]
  tmp <- unique(vardep)
  depvalues <- tmp[order(tmp)]
  cutoff=0.5
  pred_prob <- cbind(1-pred, pred)
  colnames(pred_prob) <- as.character(depvalues)

  pred_class=(pred_prob[,2]>cutoff)*1+1
  pred_class2 <- as.character(pred_class)
  for (i in 1:nrow(as.array(unique(depvalues)))){
    pred_class2[(pred_class==i)] <- as.character(depvalues[i])
  }
  pred_class <- pred_class2
  rm(pred_class2)
  pred <- cbind(pred_class,pred_prob)

  pred_all <- cbind(newdata[,as.character(formula[[2]])], pred)

  colnames(pred_all)[1:2] = rbind("dependent","pred")
  ans<- list(pred=pred_all,criterion=criterion,costRatio=costRatio)
  class(ans) <- "CSMES"
  ans
}

