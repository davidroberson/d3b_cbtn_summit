# This script is to calculate cluster stats for each k
suppressPackageStartupMessages({
  library(IntNMF)
  library(tidyverse)
})

run_clusterstats <- function(dat, output_dir, k_value) {
  # if weight is not assigned, use default
  if (is.null(wt)) {
    wt = if (is.list(dat))
      rep(1, length(dat))
    else
      1
  }
  
  # do this for each cluster
  # run Nonnegative Matrix Factorization of Multiple data using Nonnegative Alternating Least Square
  nmf_fname <- file.path(output_dir, "intnmf_fit_all.rds")
  max_k <- nmf.opt.k(dat = dat, 
            n.runs = 30, 
            n.fold = 5, 
            k.range = 2:15, 
            result = TRUE,
            make.plot = TRUE, 
            progress = TRUE, 
            st.count = 10, 
            maxiter = 100, 
            wt=if(is.list(dat)) 
              rep(1,length(dat)) else 1)

  nmf_output <- nmf.mnnals(dat = dat, 
                           k = max_k, 
                           maxiter = 200, 
                           st.count = 20, 
                           n.ini = 30, 
                           ini.nndsvd = TRUE, 
                           seed = TRUE,
                           wt=if(is.list(dat)) rep(1,length(dat)) else 1)
  
  # add modality names instead of H1, H2... for better clarity in downstream plotting
  names(nmf_output[[max_k]]$H) <- names(dat)
  
  # save best fit to file
  saveRDS(object = nmf_output[[max_k]],
          file = file.path(output_dir, "intnmf_best_fit.rds"))
  
  # return the nmf output corresponding to the most optimal k for downstream analyses
  return(nmf_output[[max_k]])
}
