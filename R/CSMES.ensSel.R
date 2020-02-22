#' CSMES Training Stage 1: Cost-Sensitive Multicriteria Ensemble Selection resulting in a Pareto frontier of candidate ensemble classifiers
#'
#' This function applies the first stage in the learning process of CSMES: optimizing Cost-Sensitive Multicriteria Ensemble
#' Selection, resulting in a Pareto frontier of equivalent candidate ensemble classifiers along two objective functions. By default, cost space is optimized
#' by optimizing false positive and false negative rates simultaneously. This results in a set of optimal ensemble classifiers, varying in the tradeoff between
#' FNR and FPR. Optionally, other objective metrics can be specified. Currently, only binary classification is supported.
#'
#' @param memberPreds matrix containing ensemble member library predictions
#' @param y Vector with true class labels. Currently, a dichotomous outcome variable is supported
#' @param obj1 Specifies the first objective metric to be minimized
#' @param obj2 Specifies the second objective metric to be minimized
#' @param selType Specifies the type of ensemble selection to be applied: \code{"selection"} for basic selection, \code{"selectionWeighted"} for weighted selection, \code{"weighted"} for weighted sum
#' @param plotting \code{TRUE} or \code{FALSE}: Should a plot be generated showing objective function values throughout the optimization process?
#' @param generations the number of population generations for nsga-II. Default is 30.
#' @param popsize the population size for nsga-II. Default is 100.
#' @return An object of the class \code{CSMES.ensSel} which is a list with the following components:
#' \item{weights}{ensemble member weights for all pareto-optimal ensemble classifiers after multicriteria ensemble selection}
#' \item{obj_values}{optimization objective values}
#' \item{pareto}{overview of pareto-optimal ensemble classifiers}
#' \item{popsize}{the population size for nsga-II}
#' \item{generarations}{the number of population generations for nsga-II}
#' \item{obj1}{Specifies the first objective metric that was minimized}
#' \item{obj2}{Specifies the second objective metric that was minimized}
#' \item{selType}{the type of ensemble selection that was applied: \code{"selection"}, \code{"selectionWeighted"} or \code{"weighted"}}
#' \item{ParetoPredictions_p}{probability predictions for pareto-optimal ensemble classifiers}
#' \item{ParetoPredictions_c}{class predictions for pareto-optimal ensebmle classifiers}
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
CSMES.ensSel <- function(memberPreds,y,obj1=c("FNR","AUCC","MSE","AUC"),obj2=c("FPR","ensSize","ensSizeSq","clAmb"),
                     selType=c("selection","selectionWeighted","weighted"),plotting=TRUE,generations=30,popsize=100) {

  minfunc <- function(x,memberPredictions,y, obj1=c("FNR","AUCC","MSE","AUC"),
                      obj2=c("FPR","ensSize","ensSizeSq","clAmb"), selType=c("selection","selectionWeighted","weighted"),plotting=TRUE) {

    calculate_MSE <- function(preds,labels) {
      y = labels;
      x = preds;
      original <- as.matrix(cbind(y,x))
      inp <- as.matrix(cbind(y,x))
      n0n1 <- nrow(inp)

      MSE<-sum((inp[,1]-inp[,2])**2)/dim(inp)[1]
      ans<- list(MSE=MSE)
      class(ans) <- "MSE"
      ans
    }
    calculate_AUCC <- function(preds,labels) {
      labels <- t(do.call("rbind", rep(list(labels), ncol(preds))))
      preds_t<-ROCR::prediction(preds,labels)
      perf <- ROCR::performance(preds_t,'ecost')

      nr_intervals=100
      values<-array(0,c(length(perf@x.values),nr_intervals+1))
      seqs<-seq(from=0,to=1,by=1/nr_intervals)
      for (i in 1:length(perf@x.values)) {
        values[i,]<-approx(perf@x.values[[i]], perf@y.values[[i]],xout=seqs)[[2]]
      }
      lower_env_coordinates<-rbind(seqs,do.call(pmin, lapply(1:nrow(values), function(i)values[i,])))

      x <- lower_env_coordinates[1,]
      y <- lower_env_coordinates[2,]
      id <- order(x)

      AUCC <- sum(diff(x[id])*rollmean(y[id],2))
      ans<- list(AUCC=AUCC)
      class(ans) <- "AUCC"
      ans
    }
    calculate_clAmb<-function(memberPredictions,preds_c){
      #class ambiguity, see Dos Santos etal
      mempreds_c<-(memberPredictions>0.5)*1
      ans<-sum((mempreds_c!=preds_c)*1)/(ncol(mempreds_c)*nrow(mempreds_c))
      ans
    }
    memberPredictions<-cbind(memberPredictions,y)
    prob_cutoff<-0.5 #cutoff value used to discretize probability predictions
    selection_weight_cutoff<-0.5 #cutoff value 0<v<1 that determines whether an ensemble member is selected or not
    x_weights <-x


    if (selType=="selection") {
      weights<-(x_weights>=selection_weight_cutoff)*1
    } else if (selType=="selectionWeighted") {
      weights<-x_weights*((x_weights>=selection_weight_cutoff)*1)
    } else if (selType=="weighted") { weights <- x_weights }

    #extra check on the weights: if sum is zero (no selection whatsoever) the first member is chosen. Otherwise the Pareto prediction
    #function will produce NaN's for that ensemble.

    if (sum(weights)==0) {
      weights[1]=1
    }

    if (sum(weights)==0) {
      Pred <- array(0,c(nrow(memberPredictions),1))
    } else {
      Pred <- (as.matrix(memberPredictions[,1:(ncol(memberPredictions)-1)])%*%weights)/sum(weights)
      #Pred <- (as.matrix(memberPredictions[,1:(ncol(memberPredictions)-1)])%*%weights)/(ncol(memberPredictions)-1)
    }
   #message(paste("selected ", sum((weights>0)*1), " from ", ncol(memberPredictions)-1, " models",sep=""))

    nr1 = sum((memberPredictions[,ncol(memberPredictions)]==1)*1)

    Pred_ordered<-Pred[order(Pred),]
    Pred_c <- (Pred>prob_cutoff)*1
    Pred_c <- factor(Pred_c,levels=0:1)
    Real <- factor(memberPredictions[,ncol(memberPredictions)],levels=0:1)

    CONF2 = table(Real, Pred_c)

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

    y <- numeric(2)
    if (obj1=="AUC") {
      y[1]<- colAUC(Pred, Real,plotROC=FALSE)
    } else if (obj1=="FNR") {
      y[1]<- FNR
    } else if (obj1=="AUCC") {
      Pred_c<-(Pred>prob_cutoff)*1
      y[1]<-calculate_AUCC(cbind(Pred_c), memberPredictions[,ncol(memberPredictions)])[[1]]
    } else if (obj1=="MSE") {
      y[1]<-calculate_MSE(cbind(Pred), memberPredictions[,ncol(memberPredictions)])[[1]]
    }
    if (obj2=="ensSize") {
      y[2]<- sum(weights)
    } else if (obj2=="ensSizeSq") {
      y[2]<- (sum(weights))**2
    } else if (obj2=="FPR") {
      y[2] <- FPR
    } else if (obj2=="clAmb"){
      if (selType=="selection") {
        part1<-(memberPredictions[,as.logical(weights)]>prob_cutoff)*1
        part2<-as.vector((Pred>prob_cutoff)*1)
        clamb <- calculate_clAmb(part1,part2)
        y[2]<-1-clamb
      } else {
        warning("class ambiguity criterion only possible for ensemble selection")
        y[2]<-1
      }
    }
    if (plotting==TRUE) {
      points(y[1],y[2],pch=".")
    }

    return (y)
  }

  target_classes <- unique(y)
  target_classes_s <- target_classes[order(target_classes)]
  newy <- as.numeric(y == target_classes_s[2])

  if (plotting==TRUE) {
    if (obj2=="ensSize") { plot(1,1,ylim=c(0,ncol(memberPreds)),xlim=c(0,1),type="p",xlab=obj1,ylab=obj2)
    } else if (obj2=="ensSizeSq") { plot(1,1,ylim=c(0,ncol(memberPreds)*2),xlim=c(0,1),type="p",xlab=obj1,ylab=obj2)
    } else {plot(1,1,xlim=c(0,1),ylim=c(0,1),type="p",xlab=obj1,ylab=obj2)}
  }

  #message("Stage 1: optimizing cost space through multicriteria ensemble selection")
  nsga2mod <- nsga2(minfunc, ncol(memberPreds), 2, memberPredictions=memberPreds,y=newy,selType=selType,obj1=obj1,obj2=obj2,
                generations=generations, popsize=popsize, lower.bounds=rep(0, ncol(memberPreds)),upper.bounds=rep(1, ncol(memberPreds)))
  selection_weight_cutoffs<-as.numeric(array(0.5,popsize))
  prob_cutoffs<-array(0.5,popsize)

  if (selType=="selection") {
    weights<-((nsga2mod$par[,1:(ncol(memberPreds))]>selection_weight_cutoffs)*1)
  } else if (selType=="selectionWeighted") {
    weights<-nsga2mod$par[,1:(ncol(memberPreds))]*((nsga2mod$par[,1:(ncol(memberPreds))]>selection_weight_cutoffs)*1)
  } else if (selType=="weighted") {
    weights<-nsga2mod$par[,1:(ncol(memberPreds))]
  }

  #Generate predictions for all Pareto-equivalent ensembles
  Pareto_predictions_p <- t(t(as.matrix(memberPreds[,1:(ncol(memberPreds))])%*%t(weights))/rowSums(weights))
  Pareto_predictions_c <- array(0, c(nrow(memberPreds),popsize))
  for (j in 1:popsize) {
    Pareto_predictions_c[,j]=cbind((Pareto_predictions_p[,j]>prob_cutoffs[j])*1)
  }

  ans<- list(weights=weights, obj_values=nsga2mod$value, pareto=nsga2mod$pareto.optimal,popsize=popsize,generations=generations,obj1=obj1,obj2=obj2,selType=selType,Pareto_predictions_p=Pareto_predictions_p,Pareto_predictions_c=Pareto_predictions_c)
  class(ans) <- "CSMES.ensSel"
  ans
}



