#' use SCimpute to impute dropout values in scRNA-seq data
#'
#' @param out_dir A character specifying the full path of the output directory, 
#' which is used to store all intermdediate and final outputs.
#' @param drop_thre A number between 0 and 1, 
#' specifying the threshold to determine dropout values.
#' @param ncores A integer specifying the number of cores used for parallel computation.

library(Matrix)
library(foreach)
library(parallel)
library(glmnet)
library(stats)
library(utils)


scimpute <- function (raw_count = raw_count, drop_thre = 0.5, celltype = FALSE, labels = NULL, ncores = detectCores()) 
{
    out_dir <- paste0(tempdir(),"/")
    
    if(celltype == TRUE & is.null(labels)){
      print("'labels' must be specified when 'celltype = TRUE'!"); stop()
    }
    print("reading in raw count matrix ...")
    count_lnorm = read_count(raw_count = raw_count, out_dir = out_dir)
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)
    print("estimating mixture models ...")

    get_mix_parameters(count = count_lnorm, point = log10(1.01), 
        path = paste0(out_dir, "parslist.rds"), ncores = ncores)
    parslist = readRDS(paste0(out_dir, "parslist.rds"))
    
    print("imputing dropout values ...")
    if (celltype == FALSE){
      count_imp = imputation_model1(count = count_lnorm, point = log10(1.01), 
                                    parslist, drop_thre = drop_thre, method = 2, ncores = ncores)
    }else{
      count_imp = imputation_model1_bytype(count = count_lnorm, labels, point = log10(1.01), parslist, 
                                           drop_thre = drop_thre, method = 2, ncores = ncores)
    }
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    return(write_count(count_imp, out_dir = out_dir))
}



#' quick re-run of SCimpute with a different \code{drop_thre}
#'
#' @param out_dir A character specifying the full path of the output directory, 
#' which is used to store all intermdediate and final outputs.
#' @param drop_thre A number between 0 and 1, 
#' specifying the threshold to determine dropout values.
#' @param celltype A logical value indicating whether cell type information is available.
#' \code{labels} must be specified if \code{celltype = TRUE}.
#' @param labels A character vector specifying the cell type of 
#' each column in the raw count matrix. Only needed when \code{celltype = TRUE}.
#' Each cell type should have at least two cells for imputation.
#' @param ncores A integer specifying the number of cores used for parallel computation.
#' @return Save the imputed count matrix to SCimpute.csv or SCimpute.txt 
#' (depending on \code{outfile}) to \code{out_dir}.
#' @export
#' @import parallel
#' @import glmnet
#' @import stats
#' @import utils
scimpute_quick <-
  function (raw_count, drop_thre = 0.5, celltype = FALSE, labels = NULL, ncores = 5) 
  {
    out_dir <- paste0(tempdir(),"/")
      
    if(celltype == TRUE & is.null(labels)){
      print("'labels' must be specified when 'celltype = TRUE'!"); stop()
    }
    print(drop_thre)
    print("reading in raw count matrix ...")
    count_lnorm = read_count(raw_count = raw_count, out_dir = out_dir)
    genenames = rownames(count_lnorm)
    cellnames = colnames(count_lnorm)
    parslist = readRDS(paste0(out_dir, "parslist.rds"))
    print("imputing dropout values ...")
    if (celltype == FALSE){
      count_imp = imputation_model1(count = count_lnorm, point = log10(1.01), 
                                    parslist, drop_thre = drop_thre, method = 2, ncores = ncores)
    }else{
      count_imp = imputation_model1_bytype(count = count_lnorm, labels = labels, point = log10(1.01), parslist, 
                                           drop_thre = drop_thre, method = 2, ncores = ncores)
    }
    count_imp = 10^count_imp - 1.01
    rownames(count_imp) = genenames
    colnames(count_imp) = cellnames
    print("writing imputed count matrix ...")
    return(write_count(count_imp, out_dir = out_dir))
  }
