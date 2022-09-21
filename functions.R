get_diff <- compiler::cmpfun(function(v) {
  v[2:length(v)]-v[1:(length(v)-1L)]
})


penalty_matrix_v1<-function(inner_knots,order=4,method="simpson"){
  
  #extend by order-1
  r<-order-1
  
  knots <- sort(c(rep(range(inner_knots), r), inner_knots),
                method="quick")
  d <- get_diff(inner_knots)
  g_ab <- splineDesign(knots, inner_knots, derivs=2,ord=order)
  
  #Using Simpson's rule 
  if(method == "simpson"){
    
    
    knots_mid <- inner_knots[-length(inner_knots)] + d / 2
    g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2, ord=order)
    g_a <- g_ab[-nrow(g_ab), ]
    g_b <- g_ab[-1, ]
    
    (crossprod(d * g_a,  g_a) + 
        4 * crossprod(d * g_ab_mid, g_ab_mid) + 
        crossprod(d * g_b, g_b)) / 6 
    
    #Using right Riemann sums
  }else if(method == "riemann"){
    crossprod(g_ab[-1,]) * d
    
  }else{
    stop(paste0("Method \"",method,"\" is not available."))
  }
}



doSVD <- function(Phi,Omega){
  Phi_SVD <- svd(Phi)
  #Omega_tilde <- t(crossprod(Phi_SVD$v, Omega %*% Phi_SVD$v)) / Phi_SVD$d
  #Omega_tilde <- t(Omega_tilde) / Phi_SVD$d
  Omega_tilde <- t(t(crossprod(Phi_SVD$v, Omega %*% Phi_SVD$v)) / Phi_SVD$d)/Phi_SVD$d
  Omega_tilde_SVD <- svd(Omega_tilde)
  U_tilde <- Phi_SVD$u %*% Omega_tilde_SVD$u
  
  structure(list(U_tilde = U_tilde,
                 Gamma = Omega_tilde_SVD$d),
            class ="svdecomp")
}

optimize_lambda_loocv_svd_v1 <- function(y,D,m, md=0.1){
  
  best_loss<-Inf
  best_lambda<-NULL
  best_f_hat<-NULL
  
  LL<-seq(0,m,md)
  for(l in seq_along(LL)){
    
    lambda<-LL[[l]]
    
    S_lambda<-matrix(1/ (1 + (lambda * D$Gamma)))
    
    #TODO below
    Sii <- D$U_tilde  %*% S_lambda
    f_hat <- D$U_tilde %*% (D$U_tilde_y / (1 + (lambda * D$Gamma) ))
    
    loss_i<-mean(((y - f_hat) / (1 - Sii))^2, na.rm=TRUE)
    
    if(loss_i < best_loss){
      best_loss<-loss_i
      best_lambda<-lambda
      best_f_hat<-f_hat
    }
  }
  structure(list(f_hat=best_f_hat, lambda=best_lambda, loss=best_loss))
}

optimize_lambda_loocv_v1<-function(y,Phi,Omega,m,md=0.1){
  
  best_loss<-Inf
  best_lambda<-NULL
  best_f_hat<-NULL
  
  LL<-seq(0,m,md)
  for(l in seq_along(LL)){
    lambda<-LL[[l]]
    
    S <- Phi %*% solve(
      crossprod(Phi) + lambda * Omega, 
      t(Phi))
    
    f_hat<-S %*% y
    Sii <- diag(S)
    loss_i <- mean(((y - f_hat ) / (1 - Sii))^2, na.rm = TRUE)
    
    if(loss_i < best_loss){
      best_loss<-loss_i
      best_lambda<-lambda
      best_f_hat<-f_hat
    }
  }
  structure(list(f_hat=best_f_hat, lambda=best_lambda, loss=best_loss))
}


smoother_v1<-function(x,y, lambda,order=4, subsample=TRUE, useSVD = FALSE, m=10){
  
  r<-order-1
  
  if(subsample){
    # p   number of B-spline basis functions
    p <- .nknots.smspl(length(x))
    # number of inner knots
    p_in <- p - 2
    inner_knots<-seq(min(x),max(x),length.out=p_in)
  }else{
    inner_knots<-x
  }
  
  
  knots <- c(rep(range(inner_knots), r), inner_knots)
  Phi<-splineDesign(knots, x)
  Omega<-penalty_matrix_v1(inner_knots)
  
  if (useSVD){
    D<-doSVD(Phi, Omega)
    D$U_tilde_y<-t(D$U_tilde)%*%y
    
    if (missing(lambda)){
      f_hat<-optimize_lambda_loocv_svd_v1(y=y,D=D,m=m)$f_hat
    }else{
      f_hat <- D$U_tilde_y / (1 + (lambda * D$Gamma))
      f_hat <- D$U_tilde %*% f_hat
    }
    
  }else{
    
    if (missing(lambda)){
      f_hat<-optimize_lambda_loocv_v1(y=y, Phi=Phi, Omega=Omega, m=m)$f_hat
    }else{
      f_hat<-Phi %*% solve(crossprod(Phi) + lambda * Omega,t(Phi) %*% y)
    }
  }    
  structure(list(xs=x,ys=f_hat),class="smoother")
}



smoother_v2<-function(data, x,y, lambda,order=4, subsample=TRUE, useSVD = FALSE, m=10){
  
  if(class(data) == "simulation_data"){
    x<-data$df$x
    y<-data$df$y
  }
  r<-order-1
  
  if(subsample){
    # p   number of B-spline basis functions
    p <- .nknots.smspl(length(x))
    # number of inner knots
    p_in <- p - 2
    inner_knots<-seq(min(x),max(x),length.out=p_in)
  }else{
    inner_knots<-x
  }
  
  
  knots <- c(rep(range(inner_knots), r), inner_knots)
  Phi<-splineDesign(knots, x)
  Omega<-penalty_matrix_v1(inner_knots)
  
  if (useSVD){
    D<-doSVD(Phi, Omega)
    D$U_tilde_y<-t(D$U_tilde)%*%y
    
    if (missing(lambda)){
      f_hat<-optimize_lambda_loocv_svd_v1(y=y,D=D,m=m)$f_hat
    }else{
      f_hat <- D$U_tilde_y / (1 + (lambda * D$Gamma))
      f_hat <- D$U_tilde %*% f_hat
    }
    
  }else{
    
    if (missing(lambda)){
      f_hat<-optimize_lambda_loocv_v1(y=y, Phi=Phi, Omega=Omega, m=m)$f_hat
    }else{
      f_hat<-Phi %*% solve(crossprod(Phi) + lambda * Omega,t(Phi) %*% y)
    }
  }    
  structure(list(xs=x,ys=f_hat),class="smoother")
}

