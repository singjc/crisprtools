# ---- Import libraries ----
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(reticulate)
library(tictoc)

if(!library(jacks, logical.return = T)){
  rjacks_zip <- paste(getwd(),'utils/JACKS-master/rjacks/jacks_0.1.0.tar.gz',sep='')
  install.packages(rjacks_zip, repos = NULL, type = "source")
  library(jacks)
} else{library(jacks)}

#' Current_CRISPR_algos
#'
#' This function allows the user to call current CRISPR algorithms to analyze their data, it includes drugZ, BAGEL and JACKS
#' drugZ (python3: <https://github.com/hart-lab/drugz>) and BAGEL (python2: <https://github.com/hart-lab/bagel>) are two python scripts,
#' developed in the Hart lab <http://www.hart-lab.org>.
#' drugZ is used to rank the resistant and sensitive genes comparing a treated sample to non-treated sample, returning a normalzied Z-score.
#' BAGEL using a bayesian approach to identify essential and non-essential genes.
#' JACKS (<https://github.com/felicityallen/JACKS/tree/master/2018_paper_materials>) is implemented in R and python, used for the inference of essentiallity
#' of genes from one or many CRISPR screens taking into account largely on how well each guide performs. JACKS was developed in the Parts Group at Sanger Institute.
#'
#'
#' @param algo_call Which algorithm to use, can take one or all. Accepted values c("drugZ", "BAGEL", "JACKS")
#' @param py_path (optional) python path to use
#' @param jacks_target_genes (Default = NULL) A vector of gene names to run JACKS on if not all genes are required. Default NULL (all genes will be used).
#' @param jacks_n_iter (Default = 5) An integer number of iterations of JACKS inference performed. Default 5.
#' @param jacks_reference_library (Default = NA) Name of the gRNA library to be used in inference ("avana","gecko2","yusa_v1o", or the path to a local grna results file). If this is specified, gRNA efficacies are not estimated, which greatly speeds up the computation. Default NULL (recalculate gRNA efficacies).
#' @param ... Other optional params to pass to either drugZ, BAGEL or JACKS
#'
#' @return Returns specfic results for either drugZ, BAGEL or JACKS
#'
#' @author drugZ was developed in the Traver Hart lab, by Gang Wang \url{<https://github.com/hart-lab/drugz>} \url{https://www.biorxiv.org/content/early/2017/12/12/232736}
#' @author BAGEL was developed in the Traver Hart lab, by Traver Hart \url{<https://github.com/hart-lab/bagel>} \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1015-8}
#' @author JACKS was developed in the Leopold Parts lab, by Fellicity Allen \url{<https://github.com/felicityallen/JACKS/tree/master/2018_paper_materials>} \url{https://www.biorxiv.org/content/early/2018/04/12/285114}
#'
#' @author Justin Sing, I take no credit for the slgorithms drugZ, BAGEL and JACKS. I merely put this script together to allow for the analysis of drugZ, BAGEL and JACKS in one framework hub. \url{https://github.com/singjc}
#'
#' @export
#'
#'
Current_CRISPR_algos <- function(algo_call, ...) {


  logParams_ <- function(params){

    tmp2 <- unlist(params)
    tmp2 <- cbind.data.frame(names(tmp2),tmp2, stringsAsFactors = F)
    cat('------- Arugments Used ------------------------------------------------\n')
    apply(tmp2,1,function(x) cat('Argument:',x[1],' | Value:',x[2],'\n\n',sep=' '))
    cat('-----------------------------------------------------------------------\n')
  }

  # Optional Arguments Check
  params = list(...)
  if (any(!grepl('py_path',names(params)))) {
    py_path=NULL
  }
  if (any(!grepl('jacks_target_genes',names(params)))) {
    jacks_target_genes=NULL
  }
  if (any(!grepl('jacks_n_iter',names(params)))) {
    jacks_n_iter = 5
  }
  if (any(!grepl('jacks_reference_library',names(params)))) {
    jacks_reference_library = NA
  }

  logParams_(params)

  available_python <- py_discover_config()

  # Get original environment path
  Org_PATH <- Sys.getenv('PATH')
  # Get python versions available
  for (i in available_python$python_versions) {
    if (any(grepl('python3.\\d$',list.files(path = dirname(i), pattern = 'python')))) {
      py3 <-  paste(dirname(i), '/',grep('python3.\\d$',list.files(path = dirname(i), pattern = 'python'), value = T),sep='')[1]
    }
    if (any(grepl('python2.\\d$',list.files(path = dirname(i), pattern = 'python')))) {
      py2 <-  paste(dirname(i), '/',grep('python2.\\d$',list.files(path = dirname(i), pattern = 'python'), value = T),sep='')[1]
    }
  }

  if ((grep('drugZ',algo_call,ignore.case = T))&&(length(grep('drugZ',algo_call,ignore.case = T))>0)) {
    if ((length(py3)>0) || (!is.null(py_path))) {
      if ((length(py3)>0)) {
        py_path <- py3
      }
      # Set Python Interpreter
      Sys.setenv(PATH = paste(py_path, Sys.getenv("PATH"), sep=":"))
      cat('Using python version:\n')
      system(paste(py_path,'--version'))

      drugZ_dir <- list.files(path = paste(getwd(),'/utils/drugz',sep=''), pattern = 'drugz.py', recursive = T, full.names = T)
      python_pipe <- paste('python', drugZ_dir, '-i', readcount_input_file,'-o',drugZ_output_file,'-f',drugZ_foldchange_file, '-c', drugZ_control, '-x', drugZ_treatmnet, sep=' ')
      cat('Running drugZ python script...\n')
      system(python_pipe, wait=T)
    } else {
      cat('There seems to be no python version 3 installed on this system. You can install python using Anacondas!\n Or pass the path (i.e. py_path = "/anaconda3/bin") to the python 3 executable if there is a python 3 avaible on this system.\n')
    }
    Sys.unsetenv("PATH")
    Sys.setenv(PATH = Org_PATH)
    cat('\n')
  }

  if ((grep('BAGEL',algo_call,ignore.case = T))&&(length(grep('BAGEL',algo_call,ignore.case = T))>0)) {

    if ((length(py2)>0) || (!is.null(py_path))) {
      if ((length(py2)>0)) {
        py_path <- py2
      }
      # Set Python Interpreter
      Sys.setenv(PATH = paste(py_path, Sys.getenv("PATH"), sep=":"))
      cat('Using python version:\n')
      system(paste(py_path,'--version'))

      BAGEL_dir <-  list.files(path = paste(getwd(),'/utils/bagel-master',sep=''), pattern = 'BAGEL.py', recursive = T, full.names = T)
      cat('Running BAGEL python script...\n')
      cat('--- Calculating foldchange...\n')
      # Run Fold Change
      python_pipe <- paste('python',BAGEL_dir,'fc -i',readcount_input_file,'-o',foldchange_out_file,'-c',T0_control_col,sep=' ')
      system(python_pipe, wait=T)
      cat('--- Calculating bayes factor...\n')
      # Bayes Factor
      reference_nonessentials <- list.files(path = paste(getwd(),'/utils/bagel-master',sep=''), pattern = '^nonessential.txt', recursive = T, full.names = T)
      reference_essentials <- list.files(path = paste(getwd(),'/utils/bagel-master',sep=''), pattern = '^essentials.txt', recursive = T, full.names = T)
      python_pipe <- paste('python',BAGEL_dir,'bf -i',paste(foldchange_out_file,'.foldchange.txt',sep=''),'-o',bf_out_file,'-e',reference_essentials,'-n',reference_nonessentials,'-c',bf_columns,sep=' ')
      system(python_pipe, wait=T)

    } else {
      cat('There seems to be no python version 2 installed on this system. You can install python using Anacondas!\n Or pass the path (i.e. py_path = "/anaconda3/envs/py27/bin") to the python 2 executable if there is a python 2 avaible on this system.\n')
    }
    Sys.unsetenv("PATH")
    Sys.setenv(PATH = Org_PATH)
    cat('\n')
  }

  if ((grep('JACKS',algo_call,ignore.case = T))&&(length(grep('JACKS',algo_call,ignore.case = T))>0)) {

    cat('JACKS: Calculating foldchange...\n')
    lfc = read_counts_from_spec_files(readcount_input_file, repmapfile, replicate_col=jacks_replicate_col,
                                      sample_col=jacks_sample_col, gene_spec_file=readcount_input_file, grna_col=jacks_grna_col,
                                      gene_col=jacks_gene_col, count_prior=jacks_count_prior, normalization=jacks_normalization, window=jacks_window,
                                      reference_sample=jacks_reference_sample)

    cat('JACKS: Inferring JACKS decomposition...\n')
    tic('JACKS: Inference decompostion Done!!: ')
    result = infer_jacks(lfc, target_genes = jacks_target_genes, n_iter = jacks_n_iter, reference_library = jacks_reference_library)
    toc()

    for (i in 1:result@colData@nrows) {
      condition <- result@colData@rownames[i]
      jacks_w <- t(result@metadata[["jacks_w"]])[,i]
      jacks_sdw <- t(result@metadata[["jacks_sdw"]])[,i]
      jacks_neg_pval <- t(result@metadata[["jacks_neg_pval"]])[,i]
      jacks_pos_pval <- t(result@metadata[["jacks_pos_pval"]])[,i]
      jacks_neg_fdr <- t(result@metadata[["jacks_neg_fdr"]])[,i]
      jacks_pos_fdr <- t(result@metadata[["jacks_pos_fdr"]])[,i]
      JACKS_df <- cbind.data.frame(jacks_w,jacks_sdw,jacks_neg_pval,jacks_pos_pval,jacks_neg_fdr,jacks_pos_fdr)
      cat('JACKS: Writing ', condition, ' results to disk...\n',sep='')
      write.table(x = JACKS_df, file = paste(jacks_saveDir,'/',condition,'_JACKS.txt',sep=''), quote = F, sep = '\t', row.names = F, col.names = T)
    }
    cat('JACKS: Writing JACKS results SummarizedExperiment object to disk...\n')
    filename <- paste(jacks_saveDir,'/JACKS_Results.rds', sep='')
    saveRDS(result, file = filename)
    cat('\n')
  }

  if (all(!grepl('drugZ|BAGEL|JACKS',algo_call,ignore.case = T))) {
    cat('You did not specify a correct algorithm call, algo_call only accepts either of the following: c("drugZ", "BAGEL", "JACKS")')
  }


}



