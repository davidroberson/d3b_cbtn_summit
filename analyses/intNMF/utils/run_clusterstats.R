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
  
  
  # Normalize by each omics type's frobenius norm
  count_data_norm <- count_data / norm(as.matrix(count_data), type="F")
  methyl_data_norm <- methyl_data / norm(as.matrix(methyl_data), type="F")
  splice_data_norm <- splice_data / norm(as.matrix(splice_data), type="F")
  snv_data_norm <- snv_data / norm(as.matrix(snv_data), type="F")
  cnv_data_norm <- cnv_data/ norm(as.matrix(cnv_data), type = "F")

  nmf_output <- nmf.mnnals(dat = list(count_data_norm, methyl_data_norm, splice_data_norm), 
                           k = k_value, 
                           maxiter = 200, 
                           st.count = 20, 
                           n.ini = 30, 
                           ini.nndsvd = TRUE, 
                           seed = TRUE,
                           wt=if(is.list(dat)) rep(1,length(dat)) else 1)
  
  # add modality names instead of H1, H2... for better clarity in downstream plotting
  names(nmf_output$H) <- names(dat)
  
  # save best fit to file
  saveRDS(object = nmf_output,
          file = file.path(output_dir, "intnmf_best_fit.rds"))
  
  # return the nmf output corresponding to the most optimal k for downstream analyses
  return(nmf_output)
}
