#cv
	cv.jmap <- function(x,y,Group) {
	  K <- length(Group)-1
	  n <- nrow(x)
	  MU <- matrix(1, n, K)
	  betahat <- list()
	  
	 #beta calculate
	  for (j in 1:K) {
	    dj_1 <- Group[j] + 1
	    dj <- Group[j+1]
	    Gj <- x[, dj_1:dj]
	    nt <- dim(Gj)[1]
	    pt <- dim(Gj)[2]
	    
	    if (as.numeric(is.null(nt))==1)
	    {
	      Hj <- Gj %*% ginv(t(Gj) %*% Gj) %*% t(Gj)
	      betahat[[j]] <- ginv(t(Gj) %*% Gj) %*% t(Gj) %*% y
	    }else if(nt < pt){
	      Hj <- Gj %*% ginv(t(Gj) %*% Gj + diag(1, pt)) %*% t(Gj)
	      betahat[[j]] <- ginv(t(Gj) %*% Gj + diag(1, pt)) %*% t(Gj) %*% y
	    }else{
	      Hj <- Gj %*% ginv(t(Gj) %*% Gj) %*% t(Gj)
	      betahat[[j]] <- ginv(t(Gj) %*% Gj) %*% t(Gj) %*% y
	    }
	 #Leave-one-out Cross-validation  
	    Dj <- diag(1 / as.vector(diag(diag(1, n) - Hj)))
	    Hj <- Dj %*% (Hj - diag(1, n)) + diag(1, n)
	    MU[, j] <- Hj %*% y
	  }
	  
	 # w calculate
	  CV <- function(w) {
	    sum((MU %*% w - y) ^ 2)
	  }
	  w <- rep(0.1, len = K)
	  fit <-optim(w,fn = CV,method = "L-BFGS-B",lower = rep(0, len = K),upper = rep(1, len = K))
	  w <- fit$par
	  CV.score <- fit$value

	  cv.fit<-list(w=w,betahat=betahat,K=K)
	  return(cv.fit)
	}

 #predict
	jma<-function(cv.fit,testx){
	  w=cv.fit$w
	  K=cv.fit$K
	  betahat=cv.fit$betahat
	  n = length(testy)
	  MU <- matrix(0, n, K)
	  pred <- c()
	  
	  
	  for (j in 1:K) {
	    dj_1 <- Group[j] + 1
	    dj <- Group[j+1]
	    Gj <- testx[, dj_1:dj]
	    MU[, j] <- Gj %*% betahat[[j]]
	  }
	  pred <- MU %*% w
	  return(pred)
	}
