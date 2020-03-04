### Main function
Newtons_method <- function(variables, object_expression,initial_point,tolerance){
  ### Gradient
  #expression
  funcD_ex <- lapply(variables, function(v) D(object_expression, v))
  ### Hessian
  #expression
  funcDD_ex <- matrix(list(), length(variables), length(variables))
  for (i in 1:length(variables)) 
    funcDD_ex[,i] = lapply(variables, function(v) D(funcD_ex[[i]], v))
  
  ### Trigger
  ### Iteration times
  Iteration = 1
  Trigger = TRUE
  x <- initial_point
  while(Trigger){
    ### Plug-in value (Gradient, Hessian)
    values <- as.list(x)
    names(values) <- variables
    
    #Gradient
    funcD_values <- as.numeric(lapply(funcD_ex, eval, env=values))
    Gradient <- as.matrix(funcD_values)
    #Hessian
    funcDD_values <- as.numeric(lapply(funcDD_ex, eval, env=values))
    Hessian <- matrix(funcDD_values, length(variables))
    
    ###Before 
    ###Newton's method
    LT <- chol(Hessian)
    L <- t(LT)
    Newton_direction <- -solve(LT) %*% solve(L) %*% Gradient
    lamda_x_sqaure <- t(solve(L) %*% Gradient) %*% (solve(L) %*% Gradient)
    Stopping_criterion <- 1/2 * lamda_x_sqaure
    
    ###Step
    t = 1
    
    if(Stopping_criterion <= tolerance){
      Trigger = FALSE
    }
    else{
      Iteration = Iteration + 1
      x = x + t*as.vector(Newton_direction)
    }
  }
  optimal_point = x
  optimal_value = eval(object_expression, envir= values)
  result = list(optimal_point = x, optimal_value = optimal_value, Iteration = Iteration)
  return(result)
} 

#Example
variables <- c("x", "y")
object_expression <- expression(25*x^2+4*y^2-20*x+4*y+5)
initial_point <- c(0.1,0.1)
tolerance <- 10^(-5)
Newtons_method(variables, object_expression,initial_point,tolerance)


