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
#' @param sim when sim=1, the beta is set as constant, when sim=2, beta follows a distribution
#' @param n_sample sample size
#' @param ita_m mean for ita
#' @param ita_sig deviation for ita
#' @param theta_m mean for theta
#' @param theta_sig deviation for theta
#' @param as shape for sampling time interval distribution
#' @param ms mean for sampling time interval distribution
#' @param fixed_beta when sim=1, the value for beta
#' @param plotsir plot the ground truth of SIR model or not
#' @param seq_length length of simulated sequence
#' @param transi_m rate of transition mutation
#' @param transv_m rate of transvertion mutation

# Outbreak_simulation
outbreaksim <- function(sim = 2,n_sample = 20, ita_m = 0, ita_sig=1e-2, theta_m = 0, theta_sig = 1e-2,
                        as = 3, ms = 1, fixed_beta = 1, plotsir = T,
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
    ## if beta follows a distribution
    ita <- rnorm(n_sample,ita_m,ita_sig)
    theta <- rnorm(n_sample,theta_m,theta_sig)
    beta <- exp(ita + theta * (1:n_sample))
  }

  # simulate infection times t
  # based on exponential distrition for the interval of infectious times
  # t[i] - t[i-1] ~ Exponential(beta*S(t[i-1])*I(t[i-1]))
  t <- array() # infection times
  t[1] <- 0
  S <- (n_sample-1):1 # number of susceptible individuals
  I <- 1:(n_sample-1) # number of infectors
  for (i in 2:n_sample) { t[i] <- t[i-1] + rexp(1,beta[i]*S[i-1]*I[i-1])}
  t <- t*365

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

# sequence_simulation
# here x is the number of mutation, based on binormal distribution
# x[[i]][t] - x[[i]][s] ~ BN(0,sigma^2(t-s));
# x[[j]][t[j]] = x[[i]][t[j]], i injects j at time t[j]
# this part of code is modified from https://github.com/cran/adegenet/blob/master/R/haploGen.R
seqsim <- function(seq_length, transi_m, transv_m, n_sample, t, h){

  # set parameters
  TRANSISET <- list('a'=as.DNAbin('g'), 'g'=as.DNAbin('a'), 'c'=as.DNAbin('t'), 't'=as.DNAbin('c'))
  TRANSVSET <- list('a'=as.DNAbin(c('c','t')), 'g'=as.DNAbin(c('c','t')), 'c'=as.DNAbin(c('a','g')), 't'=as.DNAbin(c('a','g')))

  # generate initial DNA sequence
  gendna <- function(seq_length){
    s <- sample(c('a','c','g','t'),seq_length,replace = T)
    return(s)
  }

  # create transitions for defined SNPs
  transi <- function(snp){
    res <- unlist(TRANSISET[as.character(snp)])
    class(res) <- "DNAbin"
    res <- as.character(res)
    return(res)
  }

  # create transversions for defined SNPs
  transv <- function(snp){
    res <- sapply(TRANSVSET[as.character(snp)],sample,1)
    class(res) <- "DNAbin"
    res <- as.character(res)
    return(res)
  }

  # duplicate a sequence (including possible mutations)
  # Time is the number of time units between ancestor and decendent
  seq_dupli <- function(seq, Time){
    # transitions #
    n_transi <- rbinom(n=1, size=round(seq_length*Time), prob=transi_m) # total number of transitions
    if(n_transi>0) {
      idx <- sample(1:seq_length, size=n_transi, replace=FALSE)
      seq[idx] <- transi(seq[idx])
    }
    # transversions #
    n_transv <- rbinom(n=1, size=round(seq_length*Time), prob=transv_m) # total number of transitions
    if(n_transv>0) {
      idx <- sample(1:seq_length, size=n_transv, replace=FALSE)
      seq[idx] <- transv(seq[idx])
    }
    return(seq)
  }

  # x[[j]][t[j]] = x[[i]][t[j]], i injects j at time t[j]
  x <- list() # sequence-related factor
  x[[1]] <- list()
  x[[1]][[1]] <- gendna(seq_length)
  if(length(which(h==1))>=1){
    x[[1]][[2]] <- seq_dupli(x[[1]][[1]],(t[which(h==1)[1]]-0))
    if(length(which(h==1))>1){
      for (i in 3:(length(which(h==1))+1)){
        x[[1]][[i]] <- seq_dupli(x[[1]][[i-1]], (t[which(h==1)[i-1]]-t[which(h==1)[i-2]]))
      }
    }
  }

  for (j in 2:n_sample){
    x[[j]] <- list()
    x[[j]][[1]] <- x[[h[j]]][[which(which(h==h[j])==j)]]
    if(length(which(h==j))>=1){
      x[[j]][[2]] <- seq_dupli(x[[j]][[1]], (t[which(h==j)[1]]-t[j]))
      if(length(which(h==j))>1){
        for (i in 3:(length(which(h==j))+1)){
          x[[j]][[i]] <- seq_dupli(x[[j]][[i-1]], (t[which(h==j)[i-1]]-t[which(h==j)[i-2]]))
        }
      }
    }
  }

  # sequence related factor at the observed time
  X <- sapply(x,function(v) return(v[[length(v)]]))
  X <- as.DNAbin(t(X))

  return(X)
}
