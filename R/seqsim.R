#' function to simulate sequence during an outbreak
#'
#' @title seqsim
#'
#' @param seq_length length of simulated sequence
#' @param transi_m rate of transition mutation
#' @param transv_m rate of transvertion mutation
#' @param n_sample sample size
#' @param t true infection time
#' @param h true infector
#'
#' @export

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
