protogarrote<-function(data.obj, center.interaction.x=0, scale.interaction.x=1, penalties=1, family="gaussian", nlambda=c(11,11), cv=10, outer=2, alpha1=0, alpha2=0.99){
  require(glmnet)
  x<-data.obj[["x"]]
  d<-data.obj[["d"]]
  y<-data.obj[["y"]]
  clinical<-data.obj[["clinical"]]
  interaction.x<-data.obj[["interaction.x"]]

  n<-nrow(x)
  k<-ncol(x)
  
  if(k != ncol(d)) stop("d does not match x in dimensions\n")
  if(n != nrow(d) | n!= length(y)) stop("Not equal sample size in variables\n")

  if(!is.null(interaction.x)) {
    int.x<-(interaction.x - center.interaction.x)/scale.interaction.x
    clinical <- cbind(clinical, interaction.x)  
    fit.int<-TRUE
  } else {
    fit.int<-FALSE
  }
  
  kclin <- ncol(clinical)
  
  xmat <- as.matrix(cbind(x, d, clinical))  
  if(fit.int) xmat <- as.matrix(cbind(x,  d, x*int.x, d*int.x, clinical))
  if(length(penalties)==1) penalties<-rep(penalties, ncol(xmat))
  # get lambdas
  nl1 <- nlambda[1]
  nl2 <- nlambda[2]
    
  prederr <- rsquare <- prederr2 <- matrix(0, nl1, nl2)
  for(outer.iter in 1:outer){
    folds <- rep(1:cv, each=ceiling(n/10))
    folds <- sample(folds)[1:n]
    for(inner in 1:cv){
      beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
      devel <- (1:n)[folds!=inner]
      test <- (1:n)[folds==inner]
      x.devel<-xmat[devel,]
      x.test <- xmat[test,]
      y.devel <- y[devel]
      y.test <- y[test]
      
      fit1 <- glmnet(y=y.devel, x=x.devel, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
      # second step: positive lasso: easier and safer in a loop  
      for(i in 1:nl1){
        #          lambda1<-lambda$lambda1[i]
        beta1 <- coef(fit1)[,i]
        if(!fit.int) {
          partial.eta <- sapply(1:k, function(j) x.devel[,j]*beta1[1+j] + x.devel[,k+j]*beta1[1+k+j])
          partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) x.devel[,2*k+j] * beta1[2*k+j+1]), nrow(x.devel), kclin,byrow=FALSE)
        }
        if(fit.int){
          partial.eta <- sapply(1:k, function(j) x.devel[,j]*beta1[1+j] 
                                + x.devel[,k+j]*beta1[k+1+j] 
                                + x.devel[,2*k+j]*beta1[2*k+1+j] 
                                + x.devel[,3*k+j]*beta1[3*k+1+j])
          partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) x.devel[,4*k+j] * beta1[4*k+j+1]), nrow(x.devel), kclin,byrow=FALSE)
        }
        xmat2 <- as.matrix(cbind(partial.eta, partial.clinical))
        fit2 <- glmnet(y=y.devel, x=xmat2, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
        for(ii in 1:nl2){
          if(!fit.int){
            beta[,nl2*(i-1)+ii]<-c(coef(fit2)[1,ii], # intercept
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                                   coef(fit2)[(k+2):(k+1+kclin),ii]*coef(fit1)[(2*k+2):(2*k+1+kclin),i]) # clinical 
          }
          if(fit.int){
            beta[,nl2*(i-1)+ii]<-c(coef(fit2)[1,ii], # intercept
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(2*k+2):(3*k+1),i], # x*int.x
                                   coef(fit2)[2:(k+1),ii]*coef(fit1)[(3*k+2):(4*k+1),i], # d*int.x
                                   coef(fit2)[(k+2):(k+1+kclin),ii]*coef(fit1)[(4*k+2):(4*k+1+kclin),i]) # clinical 
          }
          
        }
      } # now we have all nl1*nl2 betas               
      
      # validation: compute prediction error
      
      for(i in 1:nl1){
        for(ii in 1:nl2){
          yhat.test <- cbind(1,x.test) %*% beta[,nl2*(i-1)+ii]
          if(family=="binomial") yhat.test <- plogis(yhat.test)
          prederr[i,ii] <- prederr[i,ii]+mean((y.test-yhat.test)**2)/outer/cv
          prederr2[i,ii] <- prederr2[i,ii]+((mean((y.test-yhat.test)**2))**2)/outer/cv
          if(sd(yhat.test)!=0) rsquare[i,ii] <- cor(y.test, yhat.test)^2
          else rsquare[i,ii] <- 0
        }
      }
    }
  }
  index<-which(prederr==min(prederr), arr.ind=TRUE)
  if(length(index)>2) index<-tail(index,1)
  se.prederr<-sqrt(prederr2 - prederr**2)
  
  # now fit the garrote with all lambda1, lambda2, compute final betas and return object (including which lambda is best)
  
  beta<-matrix(0,ncol(xmat)+1,nl1*nl2)
  lambda <- matrix(0, nl1, nl2)
  # first step: ridge
  fit1 <- glmnet(y=y, x=xmat, family=family, alpha=alpha1, nlambda=nl1, penalty.factor =penalties)
  rownames(lambda) <- fit1$lambda
  
  # second step: positive lasso: easier and safer in a loop  
  for(i in 1:nl1){
    #          lambda1<-lambda$lambda1[i]
    beta1 <- coef(fit1)[,i]
    if(!fit.int) {
      partial.eta <- sapply(1:k, function(j) xmat[,j]*beta1[1+j] + xmat[,k+j]*beta1[1+k+j])
      partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) xmat[,2*k+j] * beta1[2*k+j+1]), nrow(xmat), kclin,byrow=FALSE)
    }
    if(fit.int){
      partial.eta <- sapply(1:k, function(j) xmat[,j]*beta1[1+j] 
                            + xmat[,k+j]*beta1[k+1+j] 
                            + xmat[,2*k+j]*beta1[2*k+1+j] 
                            + xmat[,3*k+j]*beta1[3*k+1+j])
      partial.clinical <- matrix(sapply(1:ncol(clinical), function(j) xmat[,4*k+j] * beta1[4*k+j+1]), nrow(xmat), kclin,byrow=FALSE)
    }
    xmat2 <- as.matrix(cbind(partial.eta, partial.clinical))
    fit2 <- glmnet(y=y, x=xmat2, family=family, alpha=alpha2, lower.limits=0, standardize=FALSE, nlambda=nl2)
    lambda[i,]<-fit2$lambda
    if(all(c(i, ii)==index)){
      xmat2.save <- xmat2
    }
    for(ii in 1:nl2){
      if(!fit.int){
        beta[,nl2*(i-1)+ii]<-c(coef(fit2)[1,ii], # intercept
                               coef(fit2)[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                               coef(fit2)[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                               coef(fit2)[(k+2):(k+1+kclin),ii]*coef(fit1)[(2*k+2):(2*k+1+kclin),i]) # clinical 
      }
      if(fit.int){
        beta[,nl2*(i-1)+ii]<-c(coef(fit2)[1,ii], # intercept
                               coef(fit2)[2:(k+1),ii]*coef(fit1)[2:(k+1),i], # x
                               coef(fit2)[2:(k+1),ii]*coef(fit1)[(k+2):(2*k+1),i], # d
                               coef(fit2)[2:(k+1),ii]*coef(fit1)[(2*k+2):(3*k+1),i], # x*int.x
                               coef(fit2)[2:(k+1),ii]*coef(fit1)[(3*k+2):(4*k+1),i], # d*int.x
                               coef(fit2)[(k+2):(k+1+kclin),ii]*coef(fit1)[(4*k+2):(4*k+1+kclin),i]) # clinical 
      }
    }
  } # now we have all nl1*nl2 beta vectors               
  coeffs<-beta[,nl2*(index[1]-1)+index[2]]
  if(!fit.int) {
    names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), colnames(clinical))

  }
  if(fit.int){
    names(coeffs)<-rownames(beta)<-c("(Intercept)", colnames(x), colnames(d), paste("i*", colnames(x), sep=""), paste("i*", colnames(d), sep=""), colnames(clinical))
  }
  res<-list(call=match.call(), fit.int=fit.int, family=family, lambda=lambda, coefficients=beta, glmnet.fit1=fit1, glmnet.fit2=fit2, k=k, kclin=kclin, cv.pred.err=prederr,
            cv.rsquare=rsquare, se.pred.err=se.prederr, center.interaction.x=center.interaction.x, scale.interaction.x=scale.interaction.x,
            fit=list(xmat=xmat, xmat2=xmat2, lambda=lambda, lambda.min=c(as.numeric(rownames(lambda))[index[1]], lambda[index[1],index[2]]), 
                     coefficients=coeffs, 
                     fitted.values=cbind(1,xmat) %*% beta[,nl2*(index[1]-1)+index[2]]))
  attr(res,"class")<-"protogarrotte"
  return(res)
}


predict.protogarrote<-function(obj, newdata, lambda="lambda.min", type="link"){
  if(!obj$fit.int){
    x<-newdata[["x"]]
    d<-newdata[["d"]]
    clinical <- newdata[["clinical"]]
    xmat<-cbind(x,d,clinical)
    if(lambda=="lambda.min"){
      x<-cbind(1, xmat[,names(obj$fit$coefficients[obj$fit$coefficients !=0])[-1]])
      beta <- obj$fit$coefficients[obj$fit$coefficients !=0]
      eta <- x %*% beta
    } else if (lambda=="all"){
        x <- cbind(1, xmat[,names(obj$fit$coefficients)[-1]])
        eta <- x %*% obj$coefficients
    }
  }
  if(obj$fit.int){
    # search for 'x' and 'd' in newdata:
    x<-newdata[["x"]]
    d<-newdata[["d"]]
    interaction.x <- (newdata[["interaction.x"]] - obj$center.interaction.x)/obj$scale.interaction.x
    
    ix<-x * interaction.x
    id<-d * interaction.x
    
    xmat<-cbind(x, d, ix, id, newdata[["clinical"]], newdata[["interaction.x"]])

    if(lambda=="lambda.min"){
      x<-cbind(1, xmat)
      beta <- coefficients.protogarrote(obj)
      eta <- x %*% beta
    } else if (lambda=="all"){
      x <- cbind(1, xmat)
      eta <- x %*% coefficients.protogarrote(obj)
    }
  }
  if(type=="response" & obj$family=="binomial") eta <- plogis(eta)
  return(eta)
}

coefficients.protogarrote<-function(obj, lambda="lambda.min"){
  if(lambda=="lambda.min") beta <- obj$fit$coefficients
  if(lambda=="all") beta <- obj$fit$beta
  return(beta)
}  

plot.protogarrote<-function(obj, highlight, jitter=20){
  if(missing(highlight)) highlight<-which(obj$cv.pred.err==min(obj$cv.pred.err), arr.ind=TRUE)
  yrange<-range(obj$cv.pred.err+obj$se.pred.err, obj$cv.pred.err-obj$se.pred.err)
  log.lambda1 <- log(as.numeric(rownames(obj$lambda)))
  dll1<-((1:ncol(obj$lambda))-(ncol(obj$lambda)+1)/2) *diff(log.lambda1)[1]/jitter
  xrange<-range(log(as.numeric(rownames(obj$lambda))))+range(dll1)
  plot(log.lambda1+dll1[1],obj$cv.pred.err[,1], type="l", col="grey", ylab="Prediction error", xlab="log lambda1", xlim=xrange, ylim=yrange)
  for(i in 2:ncol(obj$lambda)) lines(log.lambda1+dll1[i],obj$cv.pred.err[,i], type="l", col="grey")
  for(i in 1:ncol(obj$lambda)) {
    for(j in 1:nrow(obj$lambda)) {
      lines(rep(log.lambda1[j],2)+dll1[i],c(obj$cv.pred.err[j,i]-obj$se.pred.err[j,i], obj$cv.pred.err[j,i]+obj$se.pred.err[j,i]), type="l", col="grey")
    }
  }
  lines(rep(log.lambda1[highlight[1]],2)+dll1[highlight[2]],c(obj$cv.pred.err[highlight]-obj$se.pred.err[highlight], obj$cv.pred.err[highlight]+obj$se.pred.err[highlight]), type="l", col="black", lwd=1.5)
  lines(log.lambda1,obj$cv.pred.err[,highlight[2]], type="l", col="black", lwd=1.5)
}

    
plot.coefficients.protogarrote<-function(obj, order="none", scale=c(1,1,1,1), plot=TRUE){
  # obj: a protogarrote object
  # this function shows the proteomics- associated coefficients
  beta<-coefficients.protogarrote(obj, lambda="lambda.min")
  beta.hot<-beta[beta!=0]
  beta.x<-beta.hot[substr(names(beta.hot),1,2)=="x."]
  beta.d<-beta.hot[substr(names(beta.hot),1,2)=="d."]
  set.x<-substr(names(beta.x),3,10)
  set.d<-substr(names(beta.d),3,10)
  set.u<-union(set.x, set.d)
  if(!obj$fit.int){
    nam.x<-paste("x.", set.u, sep="")
    nam.d<-paste("d.", set.u, sep="")
    beta.x <-beta[nam.x]
    beta.d <- beta[nam.d]
    beta.x[is.na(beta.x)]<-0
    beta.d[is.na(beta.d)]<-0
    if(order == "x") ord<-order(beta.x)
    if(order == "d") ord<-order(beta.d)
    if(order == "none") ord<-1:length(beta.x)
    beta.x <- beta.x[ord]
    beta.d <- beta.d[ord]
    xrange <- range(beta.x*scale[1], beta.d*scale[2])
    plot(beta.x*scale[1], 1:length(beta.x), type="o", xlim=xrange, xlab="beta (scaled)", ylab="Coefficient")
    points(beta.d*scale[2], 1:length(beta.d), type="o", lty=2, col="red")
    for(i in 1:length(beta.x)) lines(xrange, c(i,i), lty=3, col="grey")
    legend("bottomright", pch=c("o","o"), lty=c(1,2), col=c("black","red"), legend=c("X","D"))
    abline(v=0, col="grey")
    return(list(beta=cbind(beta.x, beta.d), selected=set.u))
  }
  if(obj$fit.int){
    nam.ix<-paste("i*x.", set.u, sep="")
    nam.id<-paste("i*d.", set.u, sep="")
    beta.ix <- beta[nam.ix]
    beta.id <- beta[nam.id]
#    set.ix<-substr(names(beta.ix),5,12)
#    set.id<-substr(names(beta.id),5,12)
#    set.u<-union(union(union(set.x, set.d), set.ix), set.id)
    nam.x<-paste("x.", set.u, sep="")
    nam.d<-paste("d.", set.u, sep="")
    nam.ix<-paste("i*x.", set.u, sep="")
    nam.id<-paste("i*d.", set.u, sep="")
    beta.x <-beta[nam.x]
    beta.d <- beta[nam.d]
    beta.ix <-beta[nam.ix]
    beta.id <- beta[nam.id]
    beta.x[is.na(beta.x)]<-0
    beta.d[is.na(beta.d)]<-0
    beta.ix[is.na(beta.ix)]<-0
    beta.id[is.na(beta.id)]<-0
    if(order == "x") ord<-order(beta.x)
    if(order == "d") ord<-order(beta.d)
    if(order =="i*x") ord<-order(beta.ix)
    if(order =="i*d") ord<-order(beta.id)
    if(order == "none") ord<-1:length(beta.x)
    beta.x <- beta.x[ord]
    beta.d <- beta.d[ord]
    beta.ix <- beta.ix[ord]
    beta.id <- beta.id[ord]
    xrange <- range(beta.x*scale[1], beta.d*scale[2], beta.ix*scale[3], beta.id*scale[4])
    plot(beta.x*scale[1], 1:length(beta.x), type="o", xlim=xrange, xlab="beta (scaled)", ylab="Coefficient")
    points(beta.d*scale[2], 1:length(beta.d), type="o", lty=2, col="red")
    points(beta.ix*scale[3], 1:length(beta.ix), type="o", lty=2, col="blue")
    points(beta.id*scale[4], 1:length(beta.id), type="o", lty=2, col="green")
    for(i in 1:length(beta.x)) lines(xrange, c(i,i), lty=3, col="grey")
    legend("bottomright", pch=c("o","o"), lty=c(1,2), col=c("black","red", "blue", "green"), legend=c("X","D", "I*X", "I*D"))
    abline(v=0, col="grey")
    return(list(beta=cbind(beta.x, beta.d, beta.ix, beta.id), selected=set.u))
  }
}  

  


  
  
