#' outbreaksim is a tool to simulate an outbreak based on simple assupmtions
#' on the transmission, sampling and mutation
#'
#' @title outbreaksim
#'
#' \code{outbreak simulation functions}
#'
#' @author yixuan ye
#'
#' @rdname outbreaksim
#'
#' @aliases outbreaksim
#'
#' @export
#'
#' @param sim simulation pattern, when \code{sim=1}, the rate of exponential distribution (dist. of transmission interval) beta is set as constant, when \code{sim=2}, assume log(beta) follows a normal distribution
#' @param n_sample sample size
#' @param ita ita, fixed value, \code{beta = exp(ita + theta * z)}
#' @param theta theta, fixed value, \code{beta = exp(ita + theta * z)}
#' @param z_m mean for z, \code{beta = exp(ita + theta * z)}
#' @param z_sig deviation for z, \code{beta = exp(ita + theta * z)}
#' @param as shape for sampling time interval gamma distribution
#' @param ms mean for sampling time interval gamma distribution
#' @param fixed_beta when \code{sim=1}, beta is set as a fixed value
#' @param plotsir plot the ground truth of SIR model when set as True
#' @param seq_length length of simulated sequence
#' @param transi_m rate of transition mutation
#' @param transv_m rate of transvertion mutation

# Outbreak_simulation
outbreaksim <- function(n_sample = 20, ita = 1e-1, theta = 1e-1, z_m = 0, z_sig = 1,
                        as = 3, ms = 1, fixed_beta = 1, plotsir = T,sim = 2,
                        seq_length = 1e4, transi_m = NULL, transv_m = NULL){

  # import library
  library(stringr)
  library(ape)
  library(phylogram)
  library(phyloch)

  # set parameters
  if(sim==1){
    ## if beta is a constant
    beta <- rep(fixed_beta,n_sample)
  }else{
    ## if beta follows a log normal distribution
    z <- rnorm(n_sample, z_m, z_sig)
    beta <- exp(ita + theta * z)
  }

  # simulate infection times t
  # based on exponential distrition for the interval of infectious times
  # t[i] - t[i-1] ~ Exponential(beta*S(t[i-1])*I(t[i-1]))
  t <- array() # infection times
  t[1] <- 0
  S <- (n_sample-1):1 # number of susceptible individuals
  I <- 1:(n_sample-1) # number of infectors
  for (i in 2:n_sample) { t[i] <- t[i-1] + rexp(1,beta[i]*S[i-1]*I[i-1])}
  #t <- t*365

  # simulate sampling time
  # based on gamma distribution with shape as, mean ms
  t_sam <- t + rgamma(n_sample, as, as/ms)

  # simulate who infect whom
  # randomly choose previous patients as infectors
  h <- array() # infector hosts
  h[1] <- 0
  for (i in 2:n_sample) {h[i] <- sample(1:(i-1), 1, replace=T)}

  # plot SIR (ground truth)
  if(plotsir==T){
    #png(file="figure/SIR(ground_truth).png",res=1200,height=1200*6,width = 1200*10)
    matplot(t,cbind(c(n_sample-1,S),c(1,I)),lty=c(1,1),type="S",col = c('black','red'),
            xlab = "Time", ylab = "Susceptible and Infected", main = "SIR Model (ground truth)")
    legend(0,n_sample/2+2.5,c("Susceptible", "Infected"),pch=1,col=c('black','red'))
    #dev.off()
  }

  # simulate the sequence-related factor x
  if(is.null(transi_m)){ transi_m <- 1/seq_length }
  if(is.null(transv_m)){ transv_m <- transi_m/2 }
  X <- seqsim(seq_length, transi_m, transv_m, n_sample, t, h)

  # save data
  data <- list(1:n_sample,t_sam,X,t,h)
  names(data) <- c('index','sampling_time', 'sequence','infection_time','infector')
  row.names(data$sequence) <- 1:n_sample

  return(data)
}

