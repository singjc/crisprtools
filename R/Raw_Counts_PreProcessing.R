### Pre-Process raw counts ###
# Description:
#
# This script is a callable function script used to visualize raw counts data
# to check for any outliers when comparing different guides of the same gene.
# There is some natural variation one would expect, but if there is a large
# variation for one guide compared to the next, then this will cause a bias
# towards this gene being over expressed or significant. The best case scenario
# is to remove these outlying guides.
#

# Load Libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(data.table)
if(!library(tictoc, logical.return = T)){
  devtools::install_github("collectivemedia/tictoc")
} else{library(tictoc)}

# Initiate Parallel Processing
library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1


# Main Function
#' Raw Counts PreProcessing
#'
#' This script is a callable function script used to visualize raw counts data
#' to check for any outliers when comparing different guides of the same gene.
#' There is some natural variation one would expect, but if there is a large
#' variation for one guide compared to the next, then this will cause a bias
#' towards this gene being over expressed or significant. The best case scenario
#' is to remove these outlying guides.
#'
#' @param df data.frame of count data containing the sgRNA, Genes columns and 2 replicate columns for each condition
#' @param meta_cols the two columns that contain the information for the sgRNA ids and the gene ids
#' @param col_with_replicates the columns in df that contain replicates. Replicate samples need to be right next to each other for this script to work
#' @param nreplicates the number of replicates per sample
#' @param n number of quantiles/vectors to make for rnorm distribution
#' @param std the standard deviation to create a rnorm distribution
#' @param t_stat_rep_alpha the p-value cutt-off for t-test, to decide whether the replicate distributions are significantly different or not
#' @param data_process_partition the amount you want to split your data into, to make processing quicker, especially when using mclapply
#' @param consistency_cutoff the ratio of acceptable inconsistent sgRNAs to allow
#' @param printPlots option to print plots
#' @param filename File name to save filtered Results as
#' @param saveDir Directory to save results. Ex. /User/Documents
#' @param savePlots option to save plots. T or F
#' @param nrdiPlotstoSave (Not Working Yet) The number of every nth  replicate distribution intersection plot to save
#'
#' @return Writes filtered counts matrix to dataframe to disk, and sgRNA counts removed dataframe to disk, and returns both stored in a list variable
#'
#' @author Justin Sing, \url{https://github.com/singjc}
#'
#' @export
#'
#'
Raw_Counts_PreProcessing <- function(df, # data.frame of count data containing the sgRNA, Genes columns and 2 replicate columns for each condition
                                     meta_cols = c(1,2), # the two columns that contain the information for the sgRNA ids and the gene ids
                                     col_with_replicates = NULL, # the columns in df that contain replicates. Replicate samples need to be right next to each other for this script to work
                                     nreplicates = 2, # the number of replicates per sample
                                     n = 100, # number of quantiles/vectors to make for rnorm distribution
                                     std = 100, # the standard deviation to create a rnorm distribution
                                     t_stat_rep_alpha = 5e-7, # the p-value cutt-off for t-test, to decide whether the replicate distributions are significantly different or not
                                     data_process_partition = 10, # the amount you want to split your data into, to make processing quicker, especially when using mclapply
                                     consistency_cutoff = 0.3, # the ratio of acceptable inconsistent sgRNAs to allow
                                     printPlots = F, # option to print plots
                                     filename = '', # File name to save filtered Results as
                                     saveDir = NULL, # Directory to save results. Ex. /User/Documents
                                     savePlots = F, # option to save plots. T or F
                                     nrdiPlotstoSave = 1500 # (Not Working Yet) The number of every nth  replicate distribution intersection plot to save
){

  # Calculate the Bhattacharyya coefficient doing a permutation for more than 2 replicates
  bhatt.coeff.permute <- function(distributions){
    combTotest <- combn(names(distributions), m = 2)
    bc <- numeric(length = ncol(combTotest))
    for (i in 1:ncol(combTotest)) {
      bc[i] <- bhatt.coeff(distributions[[combTotest[1,i]]],distributions[[combTotest[2,i]]])
      names(bc)[i] <- paste(combTotest[1,i],'-',combTotest[2,i], sep='')
    }
    return(bc)
  } # end of bhatt.coeff.permute

  # Determine which replicates per guide have a decent amount of consistency for counts in a similar range between the replicates
  replicate_distribution_intersection <- function(x,nreplicates,n,std,saveDir,savePlots){

    sgRNA_id <- as.matrix(x[1])
    # sgRNA <- rownames(x)
    cat('Working on',sgRNA_id,'guide...\n',sep = ' ')
    # Remove sgRNA column from matrix
    x <- x[-c(1)]
    condition_names <- names(x)
    # Convert matrix to numeric data
    x <- t(as.numeric(as.matrix(x)))

    # Determine whether a plots folder and Replicate_Distribution_Intersection folder exits, otherwise, make dir
    if(savePlots==T){
      cat('Called...Saved Plots...(This could take a while)\n')
      if (!is.null(saveDir)) {
        Main_Plots_Dir <- paste(saveDir,'/Plots/',sep='')
        if (!file.exists(Main_Plots_Dir)) {
          cat('Creating Plots folder to store plots in...\n')
          dir.create(file.path(Main_Plots_Dir))
        }
        Distribution_Plots_Dir <- paste(saveDir,'/Plots/Replicate_Distribution_Intersection/',sep='')
        if (!file.exists(Distribution_Plots_Dir)) {
          cat('Creating Distribution_Plots_Dir folder to store plots in...\n')
          dir.create(file.path(Distribution_Plots_Dir))
        }
      } else {
        stop()
        cat('You did not provide a directory to save the result plots!!\n')
      }
    }

    # convert matrix to dataframe for plotting
    x_df <- (as.data.frame(x))
    # list to store plots for grid plotting
    p <- list()
    # out var to store results of true/false boolean results
    out <- list(); count = 1
    for (i in seq(from=1,to=length(x),by=nreplicates)) {

      # create random normal distribution centred around count mean with a std of user defined arg
      # dat1 <- rnorm(n,mean = x[i],sd = std)
      # dat2 <- rnorm(n,mean = x[i+1],sd = std)
      # cat('-- Current Replicates:',paste(condition_names[(i:(i+nreplicates))], collapse = ', '),'...\n',sep = ' ')
      df3 <- apply(as.matrix(x[(i:(i+nreplicates-1))]), 1, rnorm, n=n,sd = std)
      df3 <- as.data.frame(df3)
      rep_ids <- sprintf("Rep%d", 1:nreplicates)
      colnames(df3) <- c(rep_ids)
      out[[count]] <- as.list(df3)
      # Plotting action if user decides to and Print Plots
      if (savePlots==T) {
        # df3 <- cbind.data.frame(dat1,dat2)
        # colnames(df3) <- c('Rep1','Rep2')
        data.plot <- melt(df3, measure.vars = c(rep_ids))
        title <- paste(condition_names[(i:(i+nreplicates-1))], collapse = ', ')
        p[[count]] <- ggplot(data.plot, aes(x=value, fill=variable)) +
          geom_density(alpha=0.25) +
          ggtitle(title) + theme(plot.title = element_text(size=12, hjust = 0.5))
      }

      count = count + 1
    }

    # Calculate bhatt coefficient for distribution overlap percentatage, and then calculate t-test statistic to determine how significnatly different
    OV <- lapply(out,bhatt.coeff.permute)


    if (nreplicates==2) {
      pVal <- lapply(OV, function(x) 2*pnorm(-abs((((x*100)-100)/1)),sd = 10, lower.tail = T))
      pVal <- t(unlist(pVal))
      out <- apply(pVal,2,function(x) if (!is.nan(x)&&x>t_stat_rep_alpha) {return(T)}else(return(F)))
    } else {
      pVal <- lapply(OV,t.test,mu=0,alternative = "two.sided")
      pVal <- lapply(pVal,function(x) x[["p.value"]])
      pVal <- t(unlist(pVal))
      # If the two distributions intersect, then assign T, otherwise assign false to out matrix
      out <- apply(pVal,2,function(x) if (!is.nan(x)&&x<t_stat_rep_alpha) {return(T)}else(return(F)))
    }



    # # Determine whether the replicate random distributions interset
    # # cat('-- Determining if replicates have an intersection\n',sep = ' ')
    # list_df <- as.list(df3)
    # list_df <- lapply(list_df, as.integer)
    # # Rep_Intersect <- Reduce(intersect,list_df)
    # Rep_Intersect <- overlap(x=list_df,nbins=1000,plot=F)
    #
    # mean(Rep_Intersect[["OV"]])
    # sd(Rep_Intersect[["OV"]])
    #
    # pVal <- t.test(Rep_Intersect[["OV"]],mu=0,alternative = "two.sided")[['p.value']]

    # # If the two distributions intersect, then assign T, otherwise assign false to out matrix
    # out <- apply(pVal,2,function(x) if (!is.nan(x)&&x<t_stat_rep_alpha) {return(T)}else(return(F)))

    # Print Grid Plot and Save plot to directory
    if (savePlots==T) {
      org_wd <- getwd()
      setwd(Distribution_Plots_Dir)
      plotname <- sgRNA_id
      png(plotname, width = 1000, height = 1000,
          units = "px", pointsize = 12, bg = "white", res = NA)
      do.call(grid.arrange,c(p,list(top = textGrob(sgRNA_id,gp=gpar(fontsize=20,font=3)))))
      dev.off()
      setwd(org_wd)
    }
    # cat('Returning boolean result!...\n')
    return(out)
  } # End replicate_distribution_intersection Function

  # Find the main condition names between replicates
  main_condition_names <- function(col.names,nreplicates){

    out <- matrix(nrow = 1, ncol = (length(col.names)/nreplicates)); count = 1
    for (i in seq(from=1,to=length(col.names),by=nreplicates)) {
      # Get the condition name of both replicates
      s1 = col.names[i]
      s2 = col.names[i+1]
      # split the two replicate string names
      s3 <- strsplit(c(s1,s2),split="")
      # Find what characters are exactly the same to extract only the unique condition name
      new.colname <- paste((s3[[1]][s3[[1]]==s3[[2]]]),collapse = '')
      out[,count] <- new.colname
      count = count + 1
    }
    return(out)
  } # End main_condition_names Function

  # If save plots was called, check to see if Plots/ folder exists, if not create one
  if (savePlots==T) {
    if (!is.null(saveDir)) {
      Main_Plots_Dir <- paste(saveDir,'/Plots/',sep='')
      if (!file.exists(Main_Plots_Dir)) {
        cat('Creating Plots folder to store plots in...\n')
        dir.create(file.path(Main_Plots_Dir))
      }
    }else {
      stop()
      cat('You did not provide a directory to save the result plots!!\n')
    }
  }

  # Find and separate replcates from non-replicates
  if (!is.null(col_with_replicates)) {
    cat('Separating non-replicate columns from replicate columns...\n')
    org_df <- df
    meta_df <- subset(df, select = meta_cols)
    non_rep_cols <- subset(df, select = -col_with_replicates)
    genes <- non_rep_cols[,c(colnames(non_rep_cols)=='Gene' | colnames(non_rep_cols)=='GENE')]
    sgRNA <- non_rep_cols$sgRNA
    df <- subset(df, select = col_with_replicates)
    df <- cbind.data.frame(meta_df$sgRNA, df)
    colnames(df)[1] <- 'sgRNA'
    rownames(df) <- sgRNA
  } else {
    df_with_rep_cols <- df
  }

  # Boxplot of all screens before filtering
  if (printPlots==T) {
    # visualize boxplots of all the daa to look for outliers per genes per samples
    # Prepare data for plotting
    tmp <- df
    tmp$sgRNA <- NULL
    condition_names <- colnames(tmp)
    tmp <- t(tmp)
    colnames(tmp) <- sgRNA
    PlotData <- melt(tmp)
    graph_name <- "Boxplot of all screens conditions before any filtering"
    Norm <- ggplot(PlotData, aes(x=Var1, y=value)) +
      geom_boxplot() +
      theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
      ggtitle(graph_name) +
      xlab("Screen Condition") +
      theme(plot.title = element_text(hjust = 0.5,size = 16)) +
      ylab("Raw-Counts")
    if (savePlots==T) {
      org_wd <- getwd()
      setwd(Main_Plots_Dir)
      plotname <- graph_name
      png(plotname, width = 1500, height = 1000,
          units = "px", pointsize = 12, bg = "white", res = NA)
      print(Norm)
      dev.off()
      setwd(org_wd)
    } else {print(Norm)}
  }


  tmp <- df
  # rownames(tmp) <- sgRNA
  tmp.list <- as.list(as.data.frame(t(tmp)))

  seq.tmp.list <- c(1:length(tmp.list))
  split.idxs <- split(seq.tmp.list, cut(seq.tmp.list, data_process_partition, labels = FALSE))

  # Check to see whether each replicate per guide is consistent, return true or false matrix
  # true indicates that the two replicate counts are in a good range from each other to be considered
  # consistent for that one unique sample condition
  cat('Computing boolean replicate intersection...\n')

  grand_bool_list <- list()
  tic('Replicate Overlap Calculation Time ...')
  for (i in 1:length(split.idxs )) {
    tic(paste('Working on sub-split data ',i, ' ...', sep=''))
    tmp.list.split <- tmp.list[split.idxs[[i]]]
    Bool_Rep_Count_Check <- t(mclapply(tmp.list.split,replicate_distribution_intersection,nreplicates=nreplicates,n=n,std=std,saveDir=saveDir,savePlots=F, mc.cores = no_cores))
    grand_bool_list <-append(grand_bool_list,Bool_Rep_Count_Check)
    toc()
  }
  toc()
  tic.clear()

  # # Check to see whether each replicate per guide is consistent, return true or false matrix
  # # true indicates that the two replicate counts are in a good range from each other to be considered
  # # consistent for that one unique sample condition
  # cat('Computing boolean replicate intersection...\n')
  # system.time(
  #   Bool_Rep_Count_Check <- t(mclapply(tmp.list.split,replicate_distribution_intersection,nreplicates=nreplicates,n=n,std=std,saveDir=saveDir,savePlots=F, mc.cores = no_cores))
  # )

  bool_df <- as.data.frame(matrix(unlist(grand_bool_list),ncol = length(grand_bool_list[[1]]),byrow=T))

  # Obtain the unique condition names including the replicates
  condition_names <- colnames(tmp)[-c(1)]
  new.colnames <- main_condition_names(condition_names,nreplicates = nreplicates)
  colnames(bool_df) <- new.colnames

  # Calculate the ratio of how many unique sample conditions per guide have been deemed as consistent
  # over the total unique sample conditions
  cat('Computing consistency ratio...\n')
  consistency_ratio <- apply(bool_df,1,function(x) sum(x==T)/length(x))
  consistency_ratio <- as.data.frame(consistency_ratio)

  # Consistency Distribution plot
  if (printPlots==T) {
    graph_name <- 'Distribution of sgRNA consistency between replicates'
    sgRNA_distribution <- ggplot(consistency_ratio,aes(x=consistency_ratio)) + geom_density(adjust = 2) +
      ggtitle(graph_name) +
      theme(plot.title = element_text(hjust = 0.5,size = 16)) +
      geom_vline(xintercept = (consistency_cutoff), colour="#BB0000",linetype = "dashed") +
      geom_rect(aes(xmin=consistency_cutoff, xmax=max(consistency_ratio), ymin=0, ymax=(max(table(consistency_ratio))/10000)+0.3), fill = "green", alpha = 0.01)
    if (savePlots==T) {
      org_wd <- getwd()
      setwd(Main_Plots_Dir)
      plotname <- graph_name
      png(plotname, width = 1500, height = 1000,
          units = "px", pointsize = 12, bg = "white", res = NA)
      print(sgRNA_distribution)
      dev.off()
      setwd(org_wd)
    }else{print(sgRNA_distribution)}
  }

  Results <- cbind.data.frame(sgRNA,genes,bool_df,consistency_ratio)

  # Amount of guides to remove
  Percent_sgRNA_to_remove <- (sum(consistency_ratio$consistency_ratio<consistency_cutoff)/length(consistency_ratio$consistency_ratio))*100
  cat('Going to remove ',floor(Percent_sgRNA_to_remove), '% of guides from the total number of guides in the screen.\n',sep = '')
  # Determine which guides need to be removed based on a set threshold of acceptable consistency to let through
  sgRNA_to_remove <- Results[consistency_ratio$consistency_ratio<consistency_cutoff,c("sgRNA","genes")]

  # Filter out the guides that need to be removed
  cat('Filtering out sgRNA guides to be removed...\n')
  nbefore_removal <- nrow(org_df)
  Removed_sgRNAs <- df[(df$sgRNA %in% sgRNA_to_remove$sgRNA),]
  filtered.out <- df[!(df$sgRNA %in% sgRNA_to_remove$sgRNA),]
  nafter_removal <- nrow(filtered.out)
  cat('There are ', nafter_removal,' sgRNA guides left from ', nbefore_removal,' guides, after removing ', nrow(sgRNA_to_remove), ' guides.\n')

  if (!is.null(col_with_replicates)) {
    # Filter out the sgRNAs to be removed from the columns that don't have replicates
    non_rep_cols_filtered <- non_rep_cols[!(non_rep_cols$sgRNA %in% sgRNA_to_remove$sgRNA),]
    # combine the non-replicate columns with the replicate columns
    final_filtered <- cbind.data.frame(non_rep_cols_filtered,filtered.out[,-c(1)])

    # Filter out the sgRNAs to be removed from the columns that don't have replicates
    non_rep_cols_sgRNAs_removed <- non_rep_cols[(non_rep_cols$sgRNA %in% sgRNA_to_remove$sgRNA),]
    # combine the non-replicate columns with the replicate columns
    sgRNAs_removed_final <- cbind.data.frame(non_rep_cols_sgRNAs_removed,Removed_sgRNAs[,-c(1)])

  } else {
    final_filtered <- filtered.out
  }

  # Boxplot of all screens after filtering
  if (printPlots==T) {
    # visualize boxplots of all the daa to look for outliers per genes per samples
    # Prepare data for plotting
    tmp <- filtered.out
    tmp$sgRNA <- NULL
    condition_names <- colnames(tmp)
    tmp <- t(tmp)
    # colnames(tmp) <- sgRNA
    PlotData <- melt(tmp)
    graph_name <- "Boxplot of all screens conditions after filtering"
    Norm <- ggplot(PlotData, aes(x=Var1, y=value)) +
      geom_boxplot() +
      theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
      ggtitle(graph_name) +
      xlab("Screen Condition") +
      theme(plot.title = element_text(hjust = 0.5,size = 16)) +
      ylab("Raw-Counts")
    if (savePlots==T) {
      org_wd <- getwd()
      setwd(Main_Plots_Dir)
      plotname <- graph_name
      png(plotname, width = 1500, height = 1000,
          units = "px", pointsize = 12, bg = "white", res = NA)
      print(Norm)
      dev.off()
      setwd(org_wd)
    } else {print(Norm)}
  }

  # Write to disk the filtered count results
  if (!is.null(saveDir)){
    cat('Writing results to disk...\n')
    filtered_file <- paste(filename,'_filtered_count.txt',sep='')
    input_file <- paste(saveDir,'/',filtered_file,sep='')
    write.table(final_filtered, input_file, append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)

    sgRNA_removed_file <- paste(filename,'_removed_sgRNA.txt')
    input_file <- paste(saveDir,'/',sgRNA_removed_file,sep='')
    write.table(sgRNAs_removed_final, input_file, append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
  }

  # stopCluster(cl)
  return(list(final_filtered,sgRNA_to_remove))

}
