#' CSMES Training Stage 2: Extract an ensemble nomination curve (cost curve- or Brier curve-based) from a set of Pareto-optimal ensemble classifiers
#'
#' Generates an ensemble nomination curve from a set of Pareto-optimal ensemble definitions as identified through \code{CSMES.ensSel)}.
#'
#' @param ensSelModel ensemble selection model (output of \code{CSMES.ensSel})
#' @param memberPreds matrix containing ensemble member library predictions
#' @param y Vector with true class labels. Currently, a dichotomous outcome variable is supported
#' @param curveType the type of cost curve used to construct the ensemble nomination curve. Shoul be "brierCost","brierSkew" or "costCurve" (default).
#' @param method how are candidate ensemble learner predictions used to generate the ensemble nomination front? "classPreds" for class predictions (default), "probPreds" for probability predictions.
#' @param nrBootstraps optionally, the ensemble nomination curve can be generated through bootstrapping. This argument specifies the number of iterations/bootstrap samples. Default is 1.
#' @param plotting \code{TRUE} or \code{FALSE}: Should a plot be generated showing the Brier curve? Defaults to \code{FALSE}.
#' @return An object of the class \code{CSMES.ensNomCurve} which is a list with the following components:
#' \item{nomcurve}{the ensemble nomination curve}
#' \item{curves}{individual cost curves or brier curves of ensemble members}
#' \item{intervals}{resolution of the ensemble nomination curve}
#' \item{incidence}{incidence (positive rate) of the outcome variable}
#' \item{area_under_curve}{area under the ensemble nomination curve}
#' \item{method}{method used to generate the ensemble nomination front:"classPreds" for class predictions (default), "probPreds" for probability predictions}
#' \item{curveType}{the type of cost curve used to construct the ensemble nomination curve}
#' \item{nrBootstraps}{number of boostrap samples over which the ensemble nomination curve was estimated}
#' @export
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
#' @seealso \code{\link{CSMES.ensSel}}, \code{\link{CSMES.predictPareto}}, \code{\link{CSMES.predict}}
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
CSMES.ensNomCurve<-function(ensSelModel,memberPreds,y,curveType=c("costCurve","brierSkew","brierCost"),method=c("classPreds","probPreds"),plotting=FALSE,nrBootstraps=1){
  #classPreds method: every pareto ensemble delivers class predictions and thus delivers one line on the brier curve
  #probPreds method: every pareto ensemble candidate delivers prob predictions and delivers a full brier curve.
  data<-cbind(memberPreds,y)
  paretopreds<-CSMES.predictPareto(ensSelModel,data)
  popsize<-ensSelModel$popsize
  performances <- array(0,c(popsize,2))
  incidence <-sum((data[,ncol(data)]=="1")*1)/nrow(data)
  intervals=1000
  x<-seq(from = 0, to = 1, by = 1/(intervals-1))
  curves<-list()
  for (j in 1:nrBootstraps){
    curves[[j]]<-array(0,c(popsize,intervals))
    if (nrBootstraps>1) {
      sampleindices<- sample(1:nrow(data),round(nrow(data)/2),replace=TRUE)
    } else {
      sampleindices<-1:nrow(data)
    }
    Real <- factor(data[sampleindices,ncol(data)],levels=0:1)
    if (curveType=="brierCost"|curveType=="brierSkew") {
      if (method=="classPreds"){

        for (i in 1:popsize) {
          Pred <- paretopreds$Pareto_predictions_c[sampleindices,i]
          CONF2 = table(Real, Pred)

          if (length(rownames(CONF2))==1) {
            if (rownames(CONF2)==0) {
              CONF2 <- rbind(CONF2,c(0,0))
              rownames(CONF2)[2] <- 1
            } else if (rownames(CONF2)==1) {
              CONF2 <- rbind(c(0,0),CONF2)
              rownames(CONF2)[1] <- 0
            }
          }
          if (length(colnames(CONF2))==1) {
            if (colnames(CONF2)==0) {
              CONF2 <- cbind(CONF2,c(0,0))
              colnames(CONF2)[2] <- 1
            } else if (colnames(CONF2)==1) {
              CONF2 <- cbind(c(0,0),CONF2)
              colnames(CONF2)[1] <- 0
            }
          }
          FP = CONF2[1,2]
          FN = CONF2[2,1]
          FPR= CONF2[1,2]/(CONF2[1,2]+CONF2[1,1])
          FNR= CONF2[2,1]/(CONF2[2,1]+CONF2[2,2])

          performances[i,1]=FPR
          performances[i,2]=FNR
          if (curveType=="brierCost") {curves[[j]][i,]<-rbind(2*(x*(1-incidence)*FNR+(1-x)*incidence*FPR))
          } else if (curveType=="brierSkew") {curves[[j]][i,]<-FNR*(1-x)+FPR*x }

        }
      } else if (method=="probPreds"){

        for (i in 1:popsize) {
          Pred <- paretopreds$Pareto_predictions_p[sampleindices,i]
          a<-brierCurve(data[,ncol(data)],Pred,resolution=1/intervals-1)

          if (curveType=="brierCost") {curves[[j]][i,]<-a$brierCurveCost[,2]
          } else if (curveType=="brierSkew") {curves[[j]][i,]<-a$brierCurveSkew[,2] }
          rm(a)
        }
      }
    } else if (curveType=="costCurve") {
      if (method=="probPreds") {
        preds<-paretopreds$Pareto_predictions_p[sampleindices,]
      } else if (method=="classPreds") {
        preds<-paretopreds$Pareto_predictions_c[sampleindices,]
      }
      labels<-data[sampleindices,ncol(data)]
      labels <- t(do.call("rbind", rep(list(labels), ncol(preds))))
      preds_t<-ROCR::prediction(preds,labels)
      perf <- ROCR::performance(preds_t,'ecost')

      for (i in 1:length(perf@x.values)) {
        curves[[j]][i,]<-approx(perf@x.values[[i]], perf@y.values[[i]],xout=x)[[2]]
      }
    }
  }
  curves2 <- do.call(cbind, curves)
  curves2 <- array(curves2, dim=c(dim(curves[[1]]), length(curves)))
  curves<-colMeans(aperm(curves2, c(3, 1, 2)), na.rm = TRUE)
  curve<- cbind(x,as.matrix(apply( curves, 2, min)),as.matrix(apply( curves, 2, which.min)))

  #calculate score as area under thecurve
  x2 <- curve[,1]
  y2 <- curve[,2]
  id <- order(x2)
  area_under_curve<- sum(diff(x2[id])*rollmean(y2[id],2))

  if (plotting==TRUE && dim(curve[match(unique(curve[,3]),curve[,3]),,drop=FALSE])[1]>2){

    matplot(curve[,1],t(curves),type="l",ylab=curveType,xlab="Operating point",ylim=c(0,1.5*max(curve[,2])),col=3,axes=FALSE,yaxs="i", frame.plot=TRUE)
    axis(1,at=c(0,1))
    axis(2)
    lines(curve[,1],curve[,2],ylab=curveType,xlab="Operating point",ylim=c(0,1.5*max(curve[,2])),col=2)
    ccc<-cbind(curve,rbind(0,curve[1:(nrow(curve)-1),]))
    es_labels<-curve[which(ccc[,3] != ccc[,6]),,drop=FALSE]
    marg<-0.01*(1.5*max(curve[,2]))
    #abline(a=0.015,b=0)
    for (j in 1:(dim(es_labels)[1])) {
      lines(rep(es_labels[j,1],2),c(0,marg),type="l")
      lines(rep(es_labels[j,1],2),c(marg,es_labels[j,2]),type="l",pch=23, lty=3)
    }
    lines(rep(1,2),c(0,marg),type="l")
    text((rbind(es_labels[2:nrow(es_labels),1,drop=FALSE],1)+es_labels[,1,drop=FALSE])/2,0.01,labels=es_labels[,3],cex=0.7,col=2)
    #text(c(0,1),0,labels=c(0,1))
    #text(rbind(es_labels[,1,drop=FALSE],1),0,labels=round(rbind(es_labels[,1,drop=FALSE],1),digits=2),cex=0.5,col=4)
    axis(1,at=c(0,1))

    mtext(side=1,round(es_labels[2:nrow(es_labels),1,drop=FALSE],digits=2),at=es_labels[2:nrow(es_labels),1,drop=FALSE],cex=0.7,col=4,las=2,adj=1.2)
  }
  ans<- list(nomcurve=curve,curves=curves,intervals=intervals,incidence=incidence,area_under_curve=area_under_curve,method=method,curveType=curveType,nrBootstraps=nrBootstraps)
  class(ans) <- "CSMES.ensNomCurve"
  ans
}
