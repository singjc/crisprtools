### algo comparison scatter plot single ###
# Description:
#
# This script can used to compare two different datasets of the same nature.
# Meaning they can be from the same treatment/condition and have something that
# they can be compared on, i.e. Genes.
#
# This will create a scatter plot comparing two metrics from the two datasets.
# For example, comparing JACKS w essentiality score vs BAGELS baye factor score.
#

# Check for required packages existence and load or install
if(!library(data.table, logical.return = T)){
  install.packages("data.table")
} else{library(data.table)}
if(!library(ggplot2, logical.return = T)){
  install.packages("ggplot2")
} else{library(ggplot2)}
if(!library(ggrepel, logical.return = T)){
  install.packages("ggrepel")
} else{library(ggrepel)}
if(!library(ggpmisc, logical.return = T)){
  install.packages("ggpmisc")
} else{library(ggpmisc)}
if(!library(plyr, logical.return = T)){
  install.packages("plyr")
} else{library(plyr)}

# ---- Main ----
#' Scatter Plot for comparing CRISPR screen algorithms or metrics
#'
#'This script can used to compare two different datasets of the same nature.
#'Meaning they can be from the same treatment/condition and have something that
#'they can be compared on, i.e. Genes.
#'
#'This will create a scatter plot comparing two metrics from the two datasets.
#'For example, comparing JACKS w essentiality score vs BAGELS baye factor score.
#'
#' @param item1_file A full filename path to  dataset one
#' @param item2_file A full filename path to dataset two
#' @param join_on What do you want to join your two datasets on to make a comparison. i.e. 'GENE'
#' @param extract1 Which column of values to extract from dataset1. i.e. 'jacks_w'
#' @param extract2 Which column of values to extract from dataset2. i.e. 'BF'
#' @param main_title Main title of plot. i.e. 'BAGEL Bayes factor score vs JACKS essentiality score'
#' @param sub_title Dataset name for subtitle and savename i.e. "CIS_MO_T18_tkoV2_2017_04"
#' @param fraction_text fraction of outlying point identifiers taken from join_on column to plot on graph. lower value prints less labels
#' @param regression_line add a regression line to plot
#' @param col_based_on color points based on which dataset, 1 = dataset1, 2 = dataset2
#' @param column_point Which column(s) to use to identifier different points. Ex. c('jacks_neg_fdr','jacks_pos_fdr'), should be a numeric colour
#' @param alpha Threshold for what to colour different points as based on. If using FDR column, then FDR < 0.1 --> will have different colour from rest
#' @param printPlt Print plot to R graphics
#' @param savePlot Option to save plot to disk
#' @param saveDir directory you wish to save plots to
#' @param fileType type of plots you want to save
#'
#' @return A plot comparing the same dataset analysed for different metrics or algorithms, and returns scatter plot object
#'
#' @author Justin Sing, \url{https://github.com/singjc}
#'
#' @export
#'
#'
algo_comparison_scatter_plot_single <- function(item1_file,
                                                item2_file,
                                                join_on, # What do you want to join your two datasets on to make a comparison. i.e. 'GENE'
                                                extract1, # Which column of values to extract from dataset1. i.e. 'jacks_w'
                                                extract2, # Which column of values to extract from dataset2. i.e. 'BF'
                                                main_title, # Main title of plot. i.e. 'BAGEL Bayes factor score vs JACKS essentiality score'
                                                sub_title, # Dataset name for subtitle and savename i.e. "CIS_MO_T18_tkoV2_2017_04"
                                                fraction_text = 0.01, # fraction of outlying point identifiers taken from join_on column to plot on graph. lower value prints less labels
                                                regression_line =  F, # add a regression line to plot
                                                col_based_on = 1, # color points based on which dataset, 1 = dataset1, 2 = dataset2
                                                column_point, # Which column(s) to use to identifier different points. Ex. c('jacks_neg_fdr','jacks_pos_fdr'), should be a numeric colour
                                                alpha = 0.1, # Threshold for what to colour different points as based on. If using FDR column, then FDR < 0.1 --> will have different colour from rest
                                                printPlot = T, # Print plot to R graphics
                                                savePlot = F, # Option to save plot to disk
                                                saveDir = NULL, # directory you wish to save plots to
                                                fileType = c('pdf','png','tiff')) # type of plots you want to save
{

  # Set up Save Directory
  org_wd <- getwd()

  if (savePlot==T) {
    if (!is.null(saveDir)) {
      setwd(saveDir)
    } else {
      cat('Saving to Default Directory:', ord_wd, sep = " ")
    }
  }

  item = item1_file

  current_item <- gsub(ident1,"",basename(item))
  cat('--- Working on: ', current_item,'\n',sep='')
  item2 <- item2_file

  # Read in the two datasets
  dt1 <- fread(item)
  dt2 <- fread(item2)

  # Create new column to identify what points should be coloured what based on metric col_based_on
  if (!(is.null(col_based_on))) {
    if( col_based_on==1){
      dt3 <- fread(item)
    } else{
      dt3 <- fread(item2)
    }

    dt3 <- subset(dt3, select = c(join_on,column_point))
    bool_alpha <- dt3[,-c(1)] < alpha
    color_mat <- matrix('not a hit',nrow = nrow(bool_alpha),ncol = 1)
    for (col.i in 1:ncol(bool_alpha)) {
      # col.i=1
      if (col.i == 1) {
        color_mat[bool_alpha[,col.i]==T,1] <- paste(column_point[col.i],'hit',sep=' ')
      } else if (col.i == 2) {
        color_mat[bool_alpha[,col.i]==T,1] <- paste(column_point[col.i],'hit',sep=' ')
      } else{
        color_mat[bool_alpha[,col.i]==T,1] <- 'not a hit'
      }

    }

    col_df <- cbind.data.frame(dt3[,c(1)],color_mat)
    colnames(col_df)[2] <- 'Color'

  } else {
    col_df <- cbind.data.frame(subset(dt1,select = join_on),matrix('black',nrow = nrow(dt1),ncol = 1))
    colnames(col_df)[2] <- 'Color'
  }

  # Determine what colour points should be
  if (length(sort(unique(col_df$Color)))==3){
    man_col <- c("blue1","green3","black")
  } else if (length(sort(unique(col_df$Color)))==2){
    man_col <- c("blue1","black")
  } else {
    man_col <- c("black")
  }

  dt1_sub <- subset(dt1, select = c(join_on,extract1))
  dt2_sub <- subset(dt2, select = c(join_on,extract2))
  joined_dt <- join_all(list(dt1_sub,dt2_sub,col_df), by = join_on, type = 'left', match = 'all') # Join the two datasets on selected column
  joined_dt <- na.omit(joined_dt) # Remove any NA rows

  cat('--- Creating Plot...\n')

  # Create Scatter Plot
  scatter_plot <- ggplot(joined_dt,aes_string(x=extract2,y=extract1,label = join_on, colour = 'Color')) +
    geom_point() + # Call specifically geom_point to create scatter plot
    scale_colour_manual(values=man_col) + # manually colour points
    stat_dens2d_filter(geom = "text_repel", keep.fraction = fraction_text) + # Add point labels, repel text and only plot a fraction of labels
    ggtitle(main_title) +
    labs(subtitle = current_item) +
    theme(
      plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
      plot.subtitle = element_text(color="black", hjust = 0.5),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold")
    )

  # Add regression line to plot if TRUE
  if (regression_line == T) {
    scatter_plot <- scatter_plot +
      geom_smooth(method="auto", se=TRUE, linetype = 'dashed', color="red", fill="orange", fullrange=FALSE, level=0.95)
  }

  if (savePlot==T) {
    # Save Plot to disk
    cat('Saving plot: ',current_item,'\n',sep='')

    if (any(grepl('png',fileType,ignore.case = T))) {
      plotname <- paste(current_item,'.png',sep='')
      png(plotname, width = 1500, height = 1000,units = "px", pointsize = 12, bg = "white", res = NA)
      print(scatter_plot)
      dev.off()
    }

    if (any(grepl('tiff',fileType,ignore.case = T))) {
      plotname <- paste(current_item,'.tiff',sep='')
      tiff(plotname, width = 1500, height = 1000,units = "px", pointsize = 12, bg = "white", res = NA)
      print(scatter_plot)
      dev.off()
    }

    if (any(grepl('pdf',fileType,ignore.case = T))) {
      plotname <- paste(current_item,'.pdf',sep='')
      pdf(plotname,width=11,height=8.5)
      print(scatter_plot)
      dev.off()
    }
    cat('Saved location: ',saveDir,'\n',sep='')
  }

  # Set back the original working Directory
  setwd(org_wd)
  cat('\n')

  return(scatter_plot)
}




