Newtons_QP <- function(Q=NULL, p=NULL, A=NULL, b=NULL, C=NULL, d=NULL, tolerance){
  #####Predefined function
  #(1)Equality constrained
  Newtons_QP_EC <- function(Q, p, A, b, tolerance){
    ### Hessian
    Hessian <- Q
    ### KKT matrix 
    dimension <-  nrow(Hessian) + nrow(A)
    KKT <- matrix(0, nrow = dimension, ncol = dimension)
    KKT[1:nrow(Hessian), 1:ncol(Hessian)] <- Hessian
    KKT[(nrow(Hessian)+1):dimension, 1:ncol(Hessian)] <- A
    KKT[1:nrow(Hessian), (ncol(Hessian)+1):dimension] <- t(A)
    
    #initial x,nu (infeasible point is acceptable)
    x <- matrix(1, nrow = nrow(Q))
    nu <- matrix(1, nrow = nrow(A))
    ### Iteration times
    Iteration = 1
    t = 1
    Trigger = TRUE
    while(Trigger){
      ### Gradient
      gradient <- Q %*% x + p
      ###Residual 
      Residual <- matrix(0,nrow = dimension)
      Block_1 <- gradient + t(A)%*%nu
      Block_2 <- A%*%x - b
      Residual[1:nrow(Block_1),] <- Block_1
      Residual[(nrow(Block_1)+1):dimension,] <- Block_2
      direction <- solve(KKT, -Residual)
      direction_pri <- direction[1:nrow(x),]
      direction_du <- direction[(nrow(x)+1):nrow(direction),]
      
      ###Updated
      x = x + t*direction_pri
      nu = nu + t*direction_du
      ###Stopping criterion
      condition_1 <- (Q %*% x + p) + t(A)%*%nu
      condition_2 <- A%*%x - b
      check_1 <- all(condition_1 <= tolerance)
      check_2 <- all(condition_2 <= tolerance)
      stopping <- all(c(check_1, check_2))
      if(stopping){Trigger = F}
      else{Iteration = Iteration+1}
    }
    optimal_primal = as.vector(x)
    optimal_dual = as.vector(nu)
    optimal_value = as.vector((1/2)  * t(x) %*% Q %*% x + t(p) %*% x)
    result = list(Primal = optimal_primal, Dual = optimal_dual, Optimal_value = optimal_value, Iteration =Iteration)
    return(result)
  } 
  #(2)Inequality constrained
  Newtons_QP_IC <- function(Q, p, C, d, tolerance){
    ###Predifined function
    #diag_function 
    diag_func <- function(values){
      if(length(values) == 1){
        return(values)
      }
      else{return(diag(values))}
    }
    #ini_x_function
    random_generate_ini_x <- function(Matrix = C, Vector = d){
      triger = T
      while(triger){
        rnum = runif(ncol(C),-5,5)
        compute = C %*% rnum - d
        if(all(compute<0)){triger = F}
      }
      return(rnum)
    }
    ###Initial setting
    #initial x,nu,lambda 
    iequl_num <- nrow(C)
    x <- random_generate_ini_x(Matrix = C, Vector = d)
    lambda <- -1 / (C %*% x -d)
    # Hessian
    Hessian <- Q
    ### KKT matrix 
    dimension <-  nrow(Hessian) + nrow(lambda) 
    KKT <- matrix(0, nrow = dimension, ncol = dimension)
    #fixed_part
    KKT[1:nrow(Hessian), 1:nrow(Hessian)] <- Hessian #1-1
    KKT[1:nrow(Hessian), (nrow(Hessian)+1):(nrow(Hessian)+nrow(lambda))] <- t(C) #1-2
    ###Residual 
    Residual <- matrix(0,nrow = dimension)
    #surrogate duality gap
    eta <- iequl_num
    mu <- 50
    # Iteration times
    Iteration = 1
    #Line search
    Alpha <- 0.01
    Beta <- 0.5
    #Optimization iteration
    Trigger = TRUE
    while(Trigger){
      t = mu * iequl_num / eta
      # iterated KKT parts
      #2-1
      KKT[(nrow(Hessian)+1):(nrow(Hessian) + iequl_num), 1:nrow(Hessian)] <- -diag_func(as.vector(lambda))%*% C
      #2-2
      KKT[(nrow(Hessian)+1):(nrow(Hessian) + iequl_num), (nrow(Hessian)+1):(nrow(Hessian)+nrow(lambda))] <- -diag_func(as.vector(C%*%x-d) )
      ### Gradient
      gradient <- Hessian %*% x + p
      ### Residual
      r_dual <- gradient + t(C)%*%lambda 
      r_cent <- - diag_func(as.vector(lambda)) %*% (C%*%x-d) - ((1/t) * matrix(1,nrow = iequl_num))
      Residual[1:nrow(r_dual),] <- r_dual
      Residual[(nrow(r_dual)+1):(nrow(r_dual)+nrow(r_cent)),] <- r_cent
      direction <- solve(KKT, -Residual)
      direction_x <- direction[1:nrow(r_dual),]
      direction_lambda <- direction[(nrow(r_dual)+1):(nrow(r_dual)+nrow(r_cent)),]
      
      ### line search
      norm_current <- sqrt(sum(Residual*Residual))
      # Check lambda number
      if(length(lambda)==1){
        if(sum(which(direction_lambda<0))==0){
          s_max <- 1
          s <- 0.99*s_max
        }else{
          s_max <- min(1, -lambda/direction_lambda)
          s <- 0.99*s_max
        }
      }else{
        direction_lambda_minus_index <- which(direction_lambda<0)
        min_1 <- min(-lambda[direction_lambda_minus_index]/direction_lambda[direction_lambda_minus_index])
        s_max <- min(1, min_1)
        s <- 0.99*s_max
      }
      #Step_1
      s_triger_1 <- T
      while(s_triger_1){
        x_s_1 = x + s*direction_x
        stop = all( (C%*%x_s_1-d) < 0)
        if(stop){s_triger_1 = F}
        else{s = Beta * s}
      }
      #Step_2
      s_triger_2 <- T
      while(s_triger_2){
        x_s = x + s*direction_x
        lambda_s = lambda + s*direction_lambda
        r_dual_s <- (Hessian %*% x_s + p) + t(C)%*%lambda_s 
        r_cent_s <- - diag_func(as.vector(lambda_s)) %*% (C%*%x_s-d) - ((1/t) * matrix(1,nrow = iequl_num))
        r_next_step <- c(r_dual_s,r_cent_s)
        norm_next_step <- sqrt(sum(r_next_step*r_next_step))
        if(norm_next_step <= (1-Alpha*s)*norm_current){s_triger_2 = F}
        else{s = Beta * s}
      }
      ###Updated current point
      x = x + s*direction_x
      lambda = lambda + s*direction_lambda
      eta = as.vector(-t(C%*%x-d) %*% lambda)
      
      ###Stopping criterion
      stop_rpri <- sqrt(t(Q %*% x + p + t(C)%*%lambda) %*% (Q %*% x + p + t(C)%*%lambda) )
      check_1 <- (stop_rpri <= tolerance)
      check_2 <- eta <= tolerance
      stopping <- all(c(check_1, check_2))
      if(stopping){
        Trigger = F
      }
      else{
        Iteration = Iteration+1
      }
    }
    #Final result
    optimal_primal = x
    optimal_centrality= lambda
    optimal_value = (1/2)  * t(x) %*% Q %*% x + t(p) %*% x
    result = list(Primal = optimal_primal, Lambda = as.vector(optimal_centrality), Optimal_value = optimal_value,
                  Iteration =Iteration, duality_gap =eta)
    return(result)
  } 
  #(3)Both
  Newtons_QP_Both <- function(Q, p, A, b, C, d, tolerance){
    ###Predifined function
    #diag_function 
    diag_func <- function(values){
      if(length(values) == 1){
        return(values)
      }
      else{return(diag(values))}
    }
    #ini_x_function
    random_generate_ini_x <- function(Matrix = C, Vector = d){
      triger = T
      while(triger){
        rnum = runif(ncol(C),-5,5)
        compute = C %*% rnum - d
        if(all(compute<0)){triger = F}
      }
      return(rnum)
    }
    ###Initial setting
    #initial x,nu,lambda 
    iequl_num <- nrow(C)
    x <- random_generate_ini_x(Matrix = C, Vector = d)
    nu <- matrix(1, nrow = nrow(A))
    lambda <- -1 / (C %*% x -d)
    # Hessian
    Hessian <- Q
    ### KKT matrix 
    dimension <-  nrow(Hessian) + nrow(lambda) + nrow(A)
    KKT <- matrix(0, nrow = dimension, ncol = dimension)
    #fixed_part
    KKT[1:nrow(Hessian), 1:nrow(Hessian)] <- Hessian #1-1
    KKT[1:nrow(Hessian), (nrow(Hessian)+1):(nrow(Hessian)+nrow(lambda))] <- t(C) #1-2
    KKT[1:nrow(Hessian), (nrow(Hessian)+nrow(lambda)+1):dimension] <- t(A) #1-3
    KKT[(nrow(Hessian)+nrow(lambda)+1):dimension, 1:ncol(Hessian)] <- A #3-1
    ###Residual 
    Residual <- matrix(0,nrow = dimension)
    #surrogate duality gap
    eta <- iequl_num
    mu <- 50
    # Iteration times
    Iteration = 1
    #Line search
    Alpha <- 0.01
    Beta <- 0.5
    #Optimization iteration
    Trigger = TRUE
    while(Trigger){
      t = mu * iequl_num / eta
      # iterated KKT parts
      #2-1
      KKT[(nrow(Hessian)+1):(nrow(Hessian) + iequl_num), 1:nrow(Hessian)] <- -diag_func(as.vector(lambda))%*% C
      #2-2
      KKT[(nrow(Hessian)+1):(nrow(Hessian) + iequl_num), (nrow(Hessian)+1):(nrow(Hessian)+nrow(lambda))] <- -diag_func(as.vector(C%*%x-d) )
      ### Gradient
      gradient <- Hessian %*% x + p
      ### Residual
      r_dual <- gradient + t(C)%*%lambda + t(A) %*% nu
      r_cent <- - diag_func(as.vector(lambda)) %*% (C%*%x-d) - ((1/t) * matrix(1,nrow = iequl_num))
      r_pri <- A%*%x - b
      Residual[1:nrow(r_dual),] <- r_dual
      Residual[(nrow(r_dual)+1):(nrow(r_dual)+nrow(r_cent)),] <- r_cent
      Residual[(nrow(r_dual)+nrow(r_cent)+1):nrow(Residual),] <- r_pri
      direction <- solve(KKT, -Residual)
      direction_x <- direction[1:nrow(r_dual),]
      direction_lambda <- direction[(nrow(r_dual)+1):(nrow(r_dual)+nrow(r_cent)),]
      direction_nu <- direction[(nrow(r_dual)+nrow(r_cent)+1):nrow(Residual),]
      
      ### line search
      norm_current <- sqrt(sum(Residual*Residual))
      # Check lambda number
      if(length(lambda)==1){
        if(sum(which(direction_lambda<0))==0){
          s_max <- 1
          s <- 0.99*s_max
        }else{
          s_max <- min(1, -lambda/direction_lambda)
          s <- 0.99*s_max
        }
      }else{
        direction_lambda_minus_index <- which(direction_lambda<0)
        min_1 <- min(-lambda[direction_lambda_minus_index]/direction_lambda[direction_lambda_minus_index])
        s_max <- min(1, min_1)
        s <- 0.99*s_max
      }
      #Step_1
      s_triger_1 <- T
      while(s_triger_1){
        x_s_1 = x + s*direction_x
        stop = all( (C%*%x_s_1-d) < 0)
        if(stop){s_triger_1 = F}
        else{s = Beta * s}
      }
      #Step_2
      s_triger_2 <- T
      while(s_triger_2){
        x_s = x + s*direction_x
        lambda_s = lambda + s*direction_lambda
        nu_s = nu + s*direction_nu
        r_dual_s <- (Hessian %*% x_s + p) + t(C)%*%lambda_s + t(A) %*% nu_s
        r_cent_s <- - diag_func(as.vector(lambda_s)) %*% (C%*%x_s-d) - ((1/t) * matrix(1,nrow = iequl_num))
        r_pri_s <- A%*%x_s - b
        r_next_step <- c(r_dual_s,r_cent_s,r_pri_s)
        norm_next_step <- sqrt(sum(r_next_step*r_next_step))
        if(norm_next_step <= (1-Alpha*s)*norm_current){s_triger_2 = F}
        else{s = Beta * s}
      }
      ###Updated current point
      x = x + s*direction_x
      lambda = lambda + s*direction_lambda
      nu = nu + s*direction_nu
      eta = as.vector(-t(C%*%x-d) %*% lambda)
      
      ###Stopping criterion
      stop_rpri <- sqrt(t(Q %*% x + p + t(C)%*%lambda + t(A) %*% nu) %*% (Q %*% x + p + t(C)%*%lambda + t(A) %*% nu) )
      stop_rdual <- sqrt(t(A%*%x - b) %*% (A%*%x - b) )
      check_1 <- (stop_rpri <= tolerance)
      check_2 <- (stop_rdual <= tolerance)
      check_3 <- eta <= tolerance
      stopping <- all(c(check_1, check_2, check_3))
      if(stopping){
        Trigger = F
      }
      else{
        Iteration = Iteration+1
      }
    }
    #Final result
    optimal_primal = x
    optimal_dual = nu
    optimal_value = (1/2)  * t(x) %*% Q %*% x + t(p) %*% x
    result = list(Primal = optimal_primal, Dual = optimal_dual, Optimal_value = optimal_value,
                  Iteration =Iteration, duality_gap =eta)
    return(result)
  } 
  if(is.null(C) & is.null(d) ){
      Newtons_QP_EC(Q, p, A, b, tolerance)
  }else if(is.null(A) & is.null(b)){
      Newtons_QP_IC(Q, p, C, d, tolerance)
  }else{
    Newtons_QP_Both(Q, p, A, b, C, d, tolerance)
  }
} 

#1.Both
Q <- matrix(c(1,0,0, 0,2,0, 1,0,3),byrow = T, ncol = 3)
p <- matrix(c(1,-2,4) ,byrow = T, nrow = 3)
A <- matrix(c(2,3,4),ncol = 3)
b <- matrix(5,byrow = F)
C <- matrix(c(3,4,-2, 3,-2,-1, -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1), byrow = T,ncol = 3)
d <- matrix(c(10,-2,0,5,-1,5,0,5),byrow = F)
tolerance <- 10^(-5)
#Ine
Q <- matrix(c(6,2,2,2),byrow = T, ncol = 2)
p <- matrix(c(1,6), nrow = 2)
C <- matrix(c(-2,-3,-1,0,0,-1), byrow = T,ncol = 2)
d <- matrix(c(-4,0,0),byrow = F)
tolerance <- 10^(-5)

#e
Q <- matrix(c(6,2,1,2,5,2,1,2,0), 3, 3)
p <- matrix(c(8,-3,-3),nrow = 3)
A <- matrix(c(1,0,0,1,1,1),2,3)
b <- matrix(c(3,0), nrow = 2)
tolerance <- 10^(-5)


Newtons_QP(Q=Q,p=p,A=A,b=b,tolerance=tolerance)

