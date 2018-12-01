#' sigmaone function
#' 
#' This gives the sigmaone part in the variance estiamte
#' @param fit the fitting object after fitting our model
#' @return the sigma_one in the paper
sigmaone <- function(fit)
{
  N <- fit$byprod$N
  score_process <- fit$byprod$score_process
  #sum of squares of score_process
  sigma_one <- t(score_process) %*% score_process
  sigma_one <- sigma_one / N
  fit$byprod$sigma_one <- sigma_one
  sigma_one
}

#' psihat function
#' 
#' This gives the psihat part in the sigmatwo estiamte
#' @param fit the fitting object after fitting our model
#' @return psi_hat in the paper
#' @importFrom stats fitted
psihat <- function(fit)
{
  comp <- fit$byprod$comp
  binary <- fit$byprod$binary
  N <- fit$byprod$N
  X <- fit$byprod$X
  firstfit <- fit$byprod$firstfit
  Z <- fit$byprod$Z
  survtime <- fit$byprod$survtime
  Z_int <- fit$byprod$Z_int
  coef_est <- fit$coef
  #there are differences whether using binary, the derivative part
  if(binary) sum_part <- (X * (exp(-fitted(firstfit)) / (1 + exp(-fitted(firstfit))) ** 2))
  else sum_part <- X
  #a big difference between whether competing or not
  if(!comp) integral_part <- Z * survtime - Z_int
  else{
    cause <- fit$byprod$cause
    censorsurv <- fit$byprod$censorsurv
    G_int <- fit$byprod$G_int
    GZ_int <- fit$byprod$GZ_int
    integral_part <- Z * survtime - Z_int + (cause >= 2) / censorsurv * (Z * G_int - GZ_int)
  }
  #
  psi_hat <- matrix(rowSums(sapply(1:N, function(i) outer(integral_part[i, ], sum_part[i, ]))), ncol = ncol(sum_part)) * coef_est[2]
  psi_hat <- psi_hat / N
  psi_hat
}

#' sigmatwo function
#' 
#' This gives the sigmatwo part in the variance estiamte
#' @param fit the fitting object after fitting our model
#' @return the sigma_two in the paper
#' @importFrom stats vcov
sigmatwo <- function(fit)
{
  N <- fit$byprod$N
  psi_hat <- psihat(fit)
  firstfit <- fit$byprod$firstfit
  sigma_two <- N * psi_hat %*% vcov(firstfit) %*% t(psi_hat)
  sigma_two
}

#' pihatt function
#' 
#' This gives the pihatt part in the sigmathree estiamte
#' @param fit the fitting object after fitting our model
#' @return pi_hatt part in sigma_three, only appears in competing risks setting, in the paper
pihatt <- function(fit)
{
  survtime <- fit$byprod$survtime
  pi_hatt <- seq(length(survtime), 1)
  pi_hatt
}

#' qhatt function
#' 
#' This gives the qhatt part in the sigmathree estiamte
#' @param fit the fitting object after fitting our model
#' @return q_hatt part in sigma_three, only appears in competing risks setting, in the paper
qhatt <- function(fit)
{
  comp <- fit$byprod$comp
  if(!comp) q_hatt <- 0
  else{
    #the following functions are a decompostion of q_hatt
    #the main motivation of the following functions can be checked in the help file
    #please contact the author for help file
    leadingfirstpart <- function(cause, Z)
    {
      leadingfirst_part <- 0
      leadingfirst_part
    }
    
    leadingsecondpart <- function(cause, Z_bar)
    {
      leadingsecond_part <- 0
      leadingsecond_part
    }
    
    firstpart <- function(s_zero, Z, cause, censorsurv)
    {
      sum_part <- apply(Z * (cause >= 2) / censorsurv, 2, cumsum)
      sum_part <- rbind(rep(0, ncol(sum_part)), sum_part[1:(nrow(sum_part) - 1), ])
      integral_part <- rev(cumsum(rev(censorsurv * (cause == 1) / s_zero)))
      first_part <- sum_part * integral_part
      first_part
    }
    
    secondpart <- function(s_zero, Z_bar, cause, censorsurv)
    {
      sum_part <- cumsum((cause >= 2) / censorsurv)
      sum_part <- c(0, sum_part[1:(length(sum_part) - 1)])
      integral_part <- apply(apply(apply(Z_bar * censorsurv * (cause == 1) / s_zero, 2, rev), 2, cumsum), 2, rev)
      second_part <- sum_part * integral_part
      second_part
    }
    
    thirdpart <- function(survtime, coef_est, Z, cause, censorsurv)
    {
      sum_part <- apply(Z * as.numeric(Z %*% coef_est) * (cause >= 2) / censorsurv, 2, cumsum)
      sum_part <- rbind(rep(0, ncol(sum_part)), sum_part[1:(nrow(sum_part) - 1), ])
      integral_part <- Gint(survtime, censorsurv)
      third_part <- sum_part * integral_part
      third_part
    }
    
    fourthpart <- function(survtime, coef_est, Z, Z_bar, cause, censorsurv)
    {
      sum_part <- apply(Z * (cause >= 2) / censorsurv, 2, cumsum)
      sum_part <- rbind(rep(0, ncol(sum_part)), sum_part[1:(nrow(sum_part) - 1), ])
      timediff <- c(survtime[1], diff(survtime))
      integral_part <- cumsum(rev(censorsurv * as.numeric(Z_bar %*% coef_est) * timediff))
      integral_part <- rev(c(0, integral_part[1:(length(integral_part) - 1)]))
      fourth_part <- sum_part * integral_part
      fourth_part
    }
    
    fifthpart <- function(survtime, coef_est, Z_bar, cause, censorsurv)
    {
      sum_part <- cumsum((cause >= 2) / censorsurv)
      sum_part <- c(0, sum_part[1:(length(sum_part) - 1)])
      timediff <- c(survtime[1], diff(survtime))
      integral_part <- apply(apply(censorsurv * Z_bar * as.numeric(Z_bar %*% coef_est) * timediff, 2, rev), 2, cumsum)
      integral_part <- rbind(rep(0, ncol(integral_part)), integral_part[1:(nrow(integral_part) - 1), ])
      integral_part <- apply(integral_part, 2, rev)
      fifth_part <- sum_part * integral_part
      fifth_part
    }
    cause <- fit$byprod$cause
    Z <- fit$byprod$Z
    Z_bar <- fit$byprod$Z_bar
    s_zero <- fit$byprod$s_zero
    censorsurv <- fit$byprod$censorsurv
    survtime <- fit$byprod$survtime
    coef_est <- fit$coef
    leadingfirst_part <- leadingfirstpart(cause, Z)
    leadingsecond_part <- leadingsecondpart(cause, Z_bar)
    first_part <- firstpart(s_zero, Z, cause, censorsurv)
    second_part <- secondpart(s_zero, Z_bar, cause, censorsurv)
    third_part <- thirdpart(survtime, coef_est, Z, cause, censorsurv)
    fourth_part <- fourthpart(survtime, coef_est, Z, Z_bar, cause, censorsurv)
    fifth_part <- fifthpart(survtime, coef_est, Z_bar, cause, censorsurv)
    q_hatt <- leadingfirst_part - leadingsecond_part - (first_part - second_part + third_part - 2 * fourth_part + fifth_part)
  }
  q_hatt
}

#' sigmathree function
#' 
#' This gives the sigmathree part in the variance estiamte
#' @param fit the fitting object after fitting our model
#' @return the sigma_three in the paper
sigmathree <- function(fit)
{
  N <- fit$byprod$N
  pZ <- fit$byprod$pZ
  cause <- fit$byprod$cause
  sigma_three <- qhatt(fit) / pihatt(fit)
  sigma_three <- matrix(rowSums(apply(sigma_three, 1, function(x) outer(x, x)) * (cause == 0)), ncol = pZ)
  sigma_three <- sigma_three / N
  sigma_three
}

#' betavarest function
#' 
#' This gives the variance estiamte for heta
#' @param fit the fitting object after fitting our model
#' @return the variance matrix of finite dimensional parameter part
betavarest <- function(fit)
{
  N <- fit$byprod$N
  useIV <- fit$byprod$useIV
  comp <- fit$byprod$comp
  omega_inv <- fit$byprod$omega_inv
  if(!useIV && !comp){
    sigma_one <- sigmaone(fit)
    sigma_two <- 0
    sigma_three <- 0
  }
  else if(useIV && !comp){
    sigma_one <- sigmaone(fit)
    sigma_two <- sigmatwo(fit)
    sigma_three <- 0
  }
  else if(!useIV && comp){
    sigma_one <- sigmaone(fit)
    sigma_two <- 0
    sigma_three <- sigmathree(fit)
  } 
  else{
    sigma_one <- sigmaone(fit)
    sigma_two <- sigmatwo(fit)
    sigma_three <- sigmathree(fit)
  }
  pararesult <- list(sigma_one = sigma_one,
                     sigma_two = sigma_two,
                     sigma_three = sigma_three)
  fit$byprod <- append(fit$byprod, pararesult)
  fit$vcov <- omega_inv %*% (sigma_one + sigma_two + sigma_three) %*% omega_inv / N
  fit
}

#' leadpart function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return the leading part in our variance estimate
leadpart <- function(fit)
{
  N <- fit$byprod$N
  survtime <- fit$byprod$survtime
  cause <- fit$byprod$cause
  s_zero <- fit$byprod$s_zero
  lead_part = N * cumsum((cause == 1) / s_zero ^ 2)
}

#' Dhatt function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return D_hat part in the paper
Dhatt <- function(fit)
{
  Z <- fit$byprod$Z
  Z_bar <- fit$byprod$Z_bar
  cause <- fit$byprod$cause
  s_zero <- fit$byprod$s_zero
  D_hatt <- apply((Z - Z_bar) * (cause == 1) / s_zero, 2, cumsum)
  D_hatt
}

#' YGint function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return YG_int part in the paper
YGint <- function(fit)
{
  #compute the integration of G(t) over s_zero
  comp <- fit$byprod$comp
  if(comp){
    survtime <- fit$byprod$survtime
    censorsurv <- fit$byprod$censorsurv
    s_zero <- fit$byprod$s_zero
    timediff <- c(survtime[1], diff(survtime))
    temp <- cumsum(rev(timediff * censorsurv / s_zero))
    temp <- c(0, temp[1:(length(temp) - 1)])
    YG_int <- rev(temp)
    
  }
  else YG_int <- 0
  YG_int
}

#' szeroint function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return szero_int part in the paper
szeroint <- function(fit)
{
  survtime <- fit$byprod$survtime
  s_zero <- fit$byprod$s_zero
  timediff <- c(survtime[1], diff(survtime))
  szero_int <- cumsum(timediff / s_zero)
  szero_int
}

#' Yint function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return Y_int part in the paper
Yint <- function(fit)
{
  #compute the integration of reciprocal of s_zero
  szero_int <- szeroint(fit)
  YG_int <- YGint(fit)
  Y_int <- szeroint + YG_int
  Y_int
}

#' Ehatt function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @param i the ith round of our data
#' @return an integrated function with speed O(n) by recording each time result
Ehatt <- function(fit, i)
{
  coef_est <- fit$coef
  firstfit <- fit$byprod$firstfit
  X <- fit$byprod$X
  N <- fit$byprod$N
  binary <- fit$byprod$binary
  censorsurv <- fit$byprod$censorsurv
  cause <- fit$byprod$cause
  
  if(binary) A <- (X * (exp(-fitted(firstfit)) / (1 + exp(-fitted(firstfit))) ** 2))
  else A <- X
  #the following functions are mainly motivated by the help file
  #you may contact the author for the help file
  szero_int <- szeroint(fit)
  YG_int <- YGint(fit)
  firstpart <- function(fit, i){
    if(is.null(fit$byprod$Ehattfirst)){
      sum_part <- coef_est[2] * A
      integral_part <- szero_int
      integral_part_diff <- c(integral_part[1], diff(integral_part))
      first_part <- colSums(sum_part * integral_part)
      fit$byprod$Ehattfirst <- list(sum_part = sum_part,
                                    integral_part_diff = integral_part_diff,
                                    first_part = first_part)
    }
    else{
      sum_part <- fit$byprod$Ehattfirst$sum_part
      temp <- c(rep(0, i), rep(1, N - i))
      sum_part <- sum_part * temp
      integral_part_diff <- fit$byprod$Ehattfirst$integral_part_diff
      first_part <- fit$byprod$Ehattfirst$first_part
      first_part <- first_part - colSums(sum_part * integral_part_diff[i + 1])
      fit$byprod$Ehattfirst$first_part <- first_part
    }
    fit
  }
  secondpart <- function(fit, i){
    comp <- fit$byprod$comp
    if(comp){
      if(is.null(fit$byprod$Ehattsecond)){
        sum_part <- coef_est[2] * A * (cause >= 2) / censorsurv
        integral_part <- YG_int
        integral_part_diff <- c(rev(diff(rev(integral_part))), integral_part[N])
        second_part <- colSums(sum_part * integral_part)
        fit$byprod$Ehattsecond <- list(sum_part = sum_part,
                                       integral_part_diff = integral_part_diff,
                                       second_part = second_part)
      }
      else{
        fit$byprod$Ehattsecond$sum_part[i:N] <- 0
        sum_part <- fit$byprod$Ehattsecond$sum_part
        integral_part_diff <- fit$byprod$Ehattsecond$integral_part_diff
        second_part <- fit$byprod$Ehattsecond$second_part
        second_part <- second_part - colSums(sum_part * integral_part_diff[i])
        fit$byprod$Ehattsecond$second_part <- second_part
      }
    }
    else fit$byprod$Ehattsecond$second_part <- 0
    fit
  }
  fit <- firstpart(fit, i)
  fit <- secondpart(fit, i)
  first_part <- fit$byprod$Ehattfirst$first_part
  second_part <- fit$byprod$Ehattsecond$second_part
  fit$byprod$E_hatt <- first_part + second_part
  fit
}

#' EThetaEpart function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return E_Theta_E part in the paper
#' @importFrom stats vcov
EThetaEpart <- function(fit)
{
  N <- fit$byprod$N
  useIV <- fit$byprod$useIV
  comp <- fit$byprod$comp
  if(!useIV) EThetaE_part <- 0
  else{
    if(is.null(fit$byprod$EThetaE_part)){
      firstfit <- fit$byprod$firstfit
      firstvcov <- vcov(firstfit)
      EThetaE_part <- c()
      for(i in N:1){
        if(!is.null(fit$byprod$Ehattfirst)){
          first <- !fit$byprod$Ehattfirst$integral_part_diff[i + 1]
          if(comp) second <- !fit$byprod$Ehattsecond$integral_part_diff[i]
          else second <- TRUE
          if(first & second){
            EThetaE_part <- c(EThetaE_part, temp)
            next
          }
        }
        fit <- Ehatt(fit, i)
        E_hatt <- fit$byprod$E_hatt
        temp <- t(E_hatt) %*% firstvcov %*% E_hatt
        EThetaE_part <- c(EThetaE_part, temp)
      }
      EThetaE_part <- rev(EThetaE_part)
      fit$byprod$EThetaE_part <- EThetaE_part
    }
    else EThetaE_part <- fit$byprod$EThetaE_part
  }
  fit <<- fit
  EThetaE_part
}

#' Ghatt function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @param newobsz the new obtained Z value
#' @return G_hat in the paper, the only part that changes with newobsz
Ghatt <- function(fit, newobsz)
{
  survtime <- fit$byprod$survtime
  Z_int <- fit$byprod$Z_int
  G_hatt <- survtime %*% t(newobsz) - Z_int
}

#' qprimehatt function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @param i the i-th round
#' @return update the fit, add qprime_hatt
qprimehatt <- function(fit, i)
{
  #this is the q_t(u) function in the paper for estimating the variance of hazard function
  #i is the index of event time
  #the following functions are mainly motivated by the help file
  #you may contact the author for the help file
  #the first part is automatically zero
  secondpart <- function(fit, i)
  {
    N <- fit$byprod$N
    cause <- fit$byprod$cause
    s_zero <- fit$byprod$s_zero
    censorsurv <- fit$byprod$censorsurv
    #return a vector of v, the time t is fixed and decreasing each time it is called
    if(is.null(fit$byprod$qprimehattsecond)){
      sum_part <- cumsum((cause >= 2) / censorsurv)
      sum_part <- c(0, sum_part[1:(length(sum_part) - 1)])
      integral_part <- rev(cumsum(rev(censorsurv * (cause == 1) / s_zero ^ 2)))
      integral_part_diff <- censorsurv * (cause == 1) / s_zero ^ 2
      second_part <- sum_part * integral_part
      fit$byprod$qprimehattsecond <- list(sum_part = sum_part,
                                          integral_part_diff = integral_part_diff,
                                          second_part = second_part)
    }
    else{
      if(i < N){
        if(i < N - 1){
          fit$byprod$qprimehattsecond$sum_part[(i + 2):N] <- 0
        }
        sum_part <- fit$byprod$qprimehattsecond$sum_part
        integral_part_diff <- fit$byprod$qprimehattsecond$integral_part_diff
        second_part <- fit$byprod$qprimehattsecond$second_part
        second_part <- second_part - sum_part * integral_part_diff[i + 1]
        fit$byprod$qprimehattsecond$sum_part[i + 1] <- 0
        fit$byprod$qprimehattsecond$second_part <- second_part
      } 
    }
    fit <<- fit
  }
  thirdpart <- function(fit, i)
  {
    N <- fit$byprod$N
    coef_est <- fit$coef
    cause <- fit$byprod$cause
    Z <- fit$byprod$Z
    YG_int <- YGint(fit)
    censorsurv <- fit$byprod$censorsurv
    if(is.null(fit$byprod$qprimehattthird)){
      sum_part <- cumsum(Z %*% coef_est * (cause == 2) / censorsurv)
      sum_part <- c(0, sum_part[1:(length(sum_part) - 1)])
      integral_part <- YG_int
      integral_part_diff <- c(rev(diff(rev(integral_part))), YG_int[N])
      third_part <- sum_part * integral_part
      fit$byprod$qprimehattthird <- list(sum_part = sum_part,
                                         integral_part_diff = integral_part_diff,
                                         third_part = third_part)
    }
    else{
      if(i < N){
        fit$byprod$qprimehattthird$sum_part[(i + 1):N] <- 0
        sum_part <- fit$byprod$qprimehattthird$sum_part
        integral_part_diff <- fit$byprod$qprimehattthird$integral_part_diff
        third_part <- fit$byprod$qprimehattthird$third_part
        third_part <- third_part - sum_part * integral_part_diff[i + 1]
        fit$byprod$qprimehattthird$third_part <- third_part
      } 
    }
    fit <<- fit
  }
  YGZbetaint <- function(fit)
  {
    survtime <- fit$byprod$survtime
    censorsurv <- fit$byprod$censorsurv
    s_zero <- fit$byprod$s_zero
    timediff <- c(survtime[1], diff(survtime))
    coef_est <- fit$coef
    Z_bar <- fit$byprod$Z_bar
    temp <- cumsum(rev(timediff * Z_bar %*% coef_est * censorsurv / s_zero))
    temp <- c(0, temp[1:(length(temp) - 1)])
    YGZbeta_int <- rev(temp)
    YGZbeta_int
  }
  
  fourthpart <- function(fit, i)
  {
    N <- fit$byprod$N
    cause <- fit$byprod$cause
    censorsurv <- fit$byprod$censorsurv
    if(is.null(fit$byprod$qprimehattfourth)){
      sum_part <- cumsum((cause >= 2) / censorsurv)
      sum_part <- c(0, sum_part[1:(length(sum_part) - 1)])
      integral_part <- YGZbetaint(fit)
      integral_part_diff <- c(rev(diff(rev(integral_part))), integral_part[N])
      fourth_part <- sum_part * integral_part
      fit$byprod$qprimehattfourth <- list(sum_part = sum_part,
                                          integral_part_diff = integral_part_diff,
                                          fourth_part = fourth_part)
    }
    else{
      if(i < N){
        fit$byprod$qprimehattfourth$sum_part[(i + 1):N] <- 0
        sum_part <- fit$byprod$qprimehattfourth$sum_part
        integral_part_diff <- fit$byprod$qprimehattfourth$integral_part_diff
        fourth_part <- fit$byprod$qprimehattfourth$fourth_part
        fourth_part <- fourth_part - sum_part * integral_part_diff[i + 1]
        fit$byprod$qprimehattfourth$fourth_part <- fourth_part
      } 
    }
    fit <<- fit
  }
  secondpart(fit, i)
  thirdpart(fit, i)
  fourthpart(fit, i)
  second_part <- fit$byprod$qprimehattsecond$second_part
  third_part <- fit$byprod$qprimehattthird$third_part
  fourth_part <- fit$byprod$qprimehattfourth$fourth_part
  qprime_hatt <- 0 - second_part - third_part + fourth_part
  fit$byprod$qprime_hatt <- qprime_hatt
  fit <<- fit
}

#' qprimepihatt function
#' 
#' This prepares for the variance estiamte of baseline hazards function
#' @param fit the fitting object after fitting our model
#' @return the integration of squares of q_t(u) over pi(u) in the paper, introduced by competing risks
qprimepihatt <- function(fit)
{
  #compute the q_t(u) part in the paper
  if(is.null(fit$byprod$qprimepi_hatt)){
    comp = fit$byprod$comp
    if(comp){
      N <- fit$byprod$N
      cause <- fit$byprod$cause
      pi_hatt <- pihatt(fit)
      qprimepi_hatt <- c()
      for(i in N:1){
        if(!is.null(fit$byprod$qprimehattsecond)){
          second <- !fit$byprod$qprimehattsecond$integral_part_diff[i + 1]
          third <- !fit$byprod$qprimehattthird$integral_part_diff[i + 1]
          fourth <- !fit$byprod$qprimehattfourth$integral_part_diff[i + 1]
          if(second&third&fourth) {
            qprimepi_hatt <- c(qprimepi_hatt, temp)
            next
          }
        }
        fit <- qprimehatt(fit, i)
        temp <- fit$byprod$qprime_hatt / pi_hatt
        temp <- temp ^ 2
        temp <- sum(temp * (cause == 0))
        qprimepi_hatt <- c(qprimepi_hatt, temp)
      }
      qprimepi_hatt <- rev(qprimepi_hatt)
      qprimepi_hatt <- N * qprimepi_hatt
    }
    else qprimepi_hatt <- 0
    fit$byprod$qprimepi_hatt <- qprimepi_hatt
  }
  else qprimepi_hatt <- fit$byprod$qprimepi_hatt
  fit <<- fit
  qprimepi_hatt
} 

#' hazardpredvarest function
#' 
#' This gives the variance estimate for our prediction of the hazard function
#' @param newobsz the new obtained Z value
#' @param fit the fitting object after fitting our model
#' @return the variance of hazard function at each time point
hazardpredvarest <- function(newobsz, fit = NULL)
{
  N <- fit$byprod$N
  omega_inv <- fit$byprod$omega_inv
  qprimepi_hatt <- qprimepihatt(fit)
  betavar_est <- fit$vcov
  lead_part <- leadpart(fit)
  D_hatt <- Dhatt(fit)
  EThetaE_part <- EThetaEpart(fit)
  G_hatt <- Ghatt(fit, newobsz)
  hazard_aux <- list(lead_part = lead_part,
                     qprimepi_hatt = qprimepi_hatt,
                     D_hatt = D_hatt,
                     EThetaE_part = EThetaE_part,
                     G_hatt = G_hatt)
  fit$byprod <- append(fit$byprod, hazard_aux)
  #the variance estimate at each point for the hazard function
  hazardpredvar_est <- (lead_part + qprimepi_hatt + N * rowSums(G_hatt %*% betavar_est * G_hatt)
                        + EThetaE_part + 2 * rowSums(G_hatt %*% omega_inv * D_hatt)) / N
  fit$hazardpredvar_est <- hazardpredvar_est
  fit
}
