# ---------------------------------
#
# DRUGZ:  Identify drug-gene interactions in paired sample genomic perturbation screens
#
# This script was adpated from the original drugZ python version developed in the Hart Lab.
# Original drugZ python algorithm can be found: https://github.com/hart-lab/drugz
#
# For more information about drugZ's alogrithm, the bioRxiv paper can be found:
# https://www.biorxiv.org/content/biorxiv/early/2017/12/12/232736.full.pdf
#
# or you can visit Dr. Traver Hart's lab website: http://www.hart-lab.org
#
# At the time of drugZ's python script conversion to R:
#
# drugZ Python Version and Build:
# VERSION = "1.1.0.2"
# BUILD   = 110
#
# Script Conversion Author@ Justin Sing
# Date Modified: 2018-11-18
#
# ---------------------------------

# Load required packages
library(data.table)
library(dplyr)
library(plyr)
library(Matrix)
library(scales)

# ------------------------------------
# Arguments:
#
# readfile - an absolute path with filename to counts matrix data (string arrary)
# drugz_outfile - an absolute path with file name to save drugZ output (string arrary)
# control_samples - vector of column names used for controls (string vector, i.e. c('CTRL_repA', 'CTRL_repB'))
# drug_samples - vector of column names used for treatment conditions (string vector, i.e. c('Drug_repA', 'Drug_repB'))
# fc_outfile - logical T or F (default = F) to save foldchange results to disk
# remove_genes - vector of gene names (default = NULL) used to remove control genes (string vector, i.e. c('LacZ', 'luciferase', etc..))
# pseudocount - numeric value (default = 5) for the amount of pseudo counts to add to raw read counts
# minObs = numeric value (default = 1) for the minimum amount of allowable observations per gene
# half_window_size - numeric value (default = 500) for calculating empiracle bayes estimate of standard deviation iteratively scanning through a window for calculation
# index_column - numeric value (default = 1) to indicate which column contains sgRNA identifiers
# verbose - logical T or F (default = T) to print progress information to console
#
# ------------------------------------

#' Main Function
#'
#' rdrugZ
#'
#' DRUGZ:  Identify drug-gene interactions in paired sample genomic perturbation screens
#' This script was adpated from the original drugZ python version developed in the Hart Lab.
#'
#' @param readfile an absolute path with filename to counts matrix data (string arrary)
#' @param drugz_outfile an absolute path with file name to save drugZ output (string arrary)
#' @param control_samples vector of column names used for controls (string vector, i.e. c('CTRL_repA', 'CTRL_repB'))
#' @param drug_samples vector of column names used for treatment conditions (string vector, i.e. c('Drug_repA', 'Drug_repB'))
#' @param fc_outfile logical T or F (default = F) to save foldchange results to disk
#' @param remove_genes vector of gene names (default = NULL) used to remove control genes (string vector, i.e. c('LacZ', 'luciferase', etc..))
#' @param pseudocount numeric value (default = 5) for the amount of pseudo counts to add to raw read counts
#' @param minObs numeric value (default = 1) for the minimum amount of allowable observations per gene
#' @param half_window_size numeric value (default = 500) for calculating empiracle bayes estimate of standard deviation iteratively scanning through a window for calculation
#' @param index_column numeric value (default = 1) to indicate which column contains sgRNA identifiers
#' @param verbose logical T or F (default = T) to print progress information to console
#'
#' @return Writes drugZ results to disk
#'
#' @author drugZ algorithm was developed in the Traver Hart Lab, by Gang Wang \url{http://www.hart-lab.org}
#' @references Publication on drugZ can be found at: \url{https://www.biorxiv.org/content/biorxiv/early/2017/12/12/232736.full.pdf}
#' @references Original Python source code can be found at: \url{https://github.com/hart-lab/drugz}
#'
#' @author Script conversion to R: Justin Sing \url{https://github.com/singjc}
#'
#' @export
#'
#'

rdrugZ <- function( readfile, drugz_outfile, control_samples, drug_samples, fc_outfile=F,
                   remove_genes=NULL, pseudocount=5, minObs=1, half_window_size=500, index_column=1, verbose=T ) {

  # ------------------------------------
  # constants
  norm_value  = 1e7
  min_reads_thresh = 1
  # ------------------------------------

  log_ <- function( msg ){
    if ( verbose == T ) {
      cat( msg,'\n',sep='' )
    }
  }

  # Argument validation
  if ( !missing( readfile ) ) {
    if ( !file.exists( readfile ) ) { stop( paste( 'Provided readfile:', readfile, 'does not exist!', sep =' ' ) ); geterrmessage() }
  } else { stop( 'Readile argument was not supplied!!' ); geterrmessage() }
  if ( !missing( drugz_outfile ) ) {
    if ( !grepl( '*.txt', drugz_outfile ) ) { drugz_outfile = paste( drugz_outfile, '.txt', sep='' ) }
  } else { drugz_outfile = paste( dirname( readfile), '/', paste( gsub( ':|[[:space:]]', '.', Sys.time() ), 'drugZ.txt', sep = '_' ), sep = '' ) }
  if ( missing( control_samples ) ) { stop( 'Control_samples argument was not supplied!!' ); geterrmessage() }
  if ( missing( drug_samples ) ) { stop( 'Drug_samples argument was not supplied!!' ); geterrmessage() }

  num_replicates = length(control_samples)

  log_( paste( 'Control Samples:', control_samples,sep=" " ) )
  log_( paste( 'Treated Samples:', control_samples,sep=" " ) )

  # read sgRNA reads counts file
  log_( 'Reading in readfile' )
  reads <- fread( readfile )
  rownames( reads ) <- as.matrix( subset( reads, select = c( index_column ) ) )
  reads[ , c( index_column ) ] <- NULL

  # remove control genes
  # e.g. TKOv1 genes ['chr10Promiscuous','chr10Rand','chr10','EGFP','LacZ','luciferase']
  # TKOv3: 'EGFP','LacZ','luciferase'
  if ( length( remove_genes ) > 0 ) {
    log_( paste( 'Removing genes: ', remove_genes, paste = "" ) )
    reads = reads[ !( as.matrix( reads[,1] ) %in% remove_genes ), ]
  }
  numGuides = nrow( reads )
  numSamples = ncol( reads )

  # normalize to norm_value reads
  log_( 'Normalizing read counts' )
  normed <- sweep(
    x = ( norm_value * ( as.matrix( subset( reads, select = c( control_samples,drug_samples ) ) ) ) ),
    MARGIN = 2,  FUN = '/',
    STATS = ( colSums( as.matrix( subset( reads, select = c( control_samples,drug_samples ) ) ) ) ) )
  rownames( normed ) <- rownames( reads )
  ##
  # Caculate fold change with normalized reads + pseudocount
  # maintain raw read counts for future filtering
  ##
  log_( 'Processing data' )
  fc = data.frame( reads[,1] )
  fc$Order = matrix( c( 1:numGuides ), nrow = numGuides, ncol = 1 )
  rownames( fc ) = rownames( reads )
  colnames( fc )[1] <- 'GENE'
  for ( k in 1:length( control_samples ) ) {
    log_( paste( 'Calculating raw fold change for replicate', k, sep=" " ) )
    fc[ , control_samples[k] ] = subset( reads, select = control_samples[k] )
    fc[ , drug_samples[k] ] = subset( reads, select = drug_samples[k] )
    fc[ ,paste( 'fc_',( k-1 ), sep='') ] = log2( ( subset( normed, select = drug_samples[k] ) + pseudocount ) / ( subset( normed, select = control_samples[k] ) + pseudocount) )
    #
    # sort guides by readcount, descending:
    #
    fc <- fc[ order( subset( fc, select = control_samples[k] ), decreasing = T), ]
    #
    # add columns for eb_mean, eb_std
    #
    # eb_mean_samplid = 'eb_mean_{0}'.format(k)
    eb_std_samplid  = paste( 'eb_std_', (k-1),sep='' )
    # fc[eb_mean_samplid] = np.zeros(numGuides)
    fc[ eb_std_samplid ] = matrix( 0, nrow = numGuides )
    #
    # get mean, std of fold changes based on 800 nearest fc
    #
    log_( paste( 'Caculating smoothed Epirical Bayes estimates of stdev for replicate', (k), sep = " " ) )
    #
    # initialize element at index 250
    #
    # do not mean-center. fc of 0 should be z=score of 0.
    # ebmean = fc.iloc[0:500]['fc_{0}'.format(k)].mean()
    # fc[eb_mean_samplid][0:250] = ebmean
    #
    # parameter "half_window_size" added
    half_window_size = half_window_size
    ebstd  = sd( fc[ (1):( half_window_size*2 ), paste( 'fc_',( k-1 ),sep='' ) ] )
    fc[1:half_window_size,eb_std_samplid]  = ebstd
    #
    # from 250..(end-250), calculate mean/std, update if >= previous (monotone smoothing)
    #
    for ( i in seq( half_window_size, numGuides-half_window_size+25, 25 ) ) {
      print(i)
      # every 25th guide, calculate stdev. binning/smoothing approach.
      if ( ( i+half_window_size ) > numGuides ) {
        ebstd  = sd( as.matrix( subset( fc, subset = ( 1:numGuides ) %in% ( ( i-half_window_size+1 ):( numGuides ) ), select = paste( 'fc_',( k-1 ), sep = '' ) ) ) )
      } else {
        ebstd  = sd( as.matrix( subset( fc, subset = ( 1:numGuides ) %in% ( ( i-half_window_size ):( i+half_window_size ) ), select = paste( 'fc_',( k-1 ), sep = '' ) ) ) )
      }
      if ( ebstd >= fc[i-1, eb_std_samplid] ) {
        fc[ i:(i+25), eb_std_samplid ] = ebstd # set new std in whole step size (25)
      } else {
        fc[ i:(i+25), eb_std_samplid ] = as.numeric( fc[i-1, eb_std_samplid] )
      }
    }
    #
    # set ebmean, ebstd for bottom half-window set of guides
    #
    #log_('Smoothing estimated std for replicate {0}'.format(k+1))
    #fc[eb_mean_samplid][numGuides-250:] = fc.iloc[numGuides-251][eb_mean_samplid]
    fc[ ( numGuides-half_window_size ):numGuides,eb_std_samplid ] = as.numeric( fc[ ( numGuides-( half_window_size+1 ) ),eb_std_samplid ] )
    #
    # calc z score of guide
    #
    log_( paste( 'Caculating Zscores for replicate', (k), sep = " " ) )
    #fc['Zlog_fc_{0}'.format(k)] = (fc['fc_{0}'.format(k)] - fc[eb_mean_samplid]) / fc[eb_std_samplid]
    fc[ paste( 'Zlog_fc_',k-1,sep="" ) ] = fc[ paste( 'fc_',(k-1),sep='' ) ] / fc[ eb_std_samplid ]

    fc <- fc[ order( fc$Order ), ]
  }
  fc$Order <- NULL
  ##
  # write fc file as intermediate output
  ##
  if ( !(isTRUE( fc_outfile )) ) { log_( 'Writing fold-change file to disk' ); write.table(x = fc, file = paste( gsub( '*.txt','',drugz_outfile ),'-foldchange.txt', sep = '' ), quote = F, sep='\t', dec = '.', row.names = F, col.names = T) }
  ##
  # sum guide-level zscores to gene-level drugz scores. keep track of how many elements (fold change observations) were summed.
  ##
  log_( 'Combining drugZ scores' )
  # get unique list of genes in the data set
  usedColumns = sprintf( "Zlog_fc_%d", 0:(num_replicates-1) )
  drugz = ddply( melt( fc, id.vars = 'GENE', measure.vars = usedColumns ), .( GENE ), summarize, sumZ = sum( value ), numObs = nnzero( value ) )
  #
  #
  log_( 'Writing output file' )
  #
  # calculate normZ, pvals (from normal dist), and fdrs (by benjamini & hochberg).
  #
  drugz_minobs = drugz[ drugz$numObs>=minObs, ]
  numGenes = nrow( drugz_minobs )
  numCols = ncol( drugz_minobs )
  ##
  # normalize sumZ by number of observations
  ##
  drugz_minobs[,'normZ'] = drugz_minobs[,'sumZ'] / sqrt( drugz_minobs[,'numObs'])
  #
  # renormalize to fit uniform distribution of null p-vals
  #
  drugz_minobs['normZ'] = ( drugz_minobs$normZ - mean( drugz_minobs$normZ ) ) / sd( drugz_minobs$normZ )
  #
  # rank by normZ (ascending) to identify synthetic interactions
  #
  drugz_minobs = drugz_minobs[ order( subset( drugz_minobs, select = normZ ), decreasing = F ), ]
  drugz_minobs['pval_synth'] = pnorm( (drugz_minobs$normZ * -1), lower.tail = F )
  drugz_minobs['rank_synth'] = seq( from = 1, to = numGenes, by = 1)
  drugz_minobs['fdr_synth'] = ( drugz_minobs$pval_synth*numGenes ) / ( drugz_minobs$rank_synth )
  #
  # rerank by normZ (descending) to identify suppressor interactions
  #
  drugz_minobs = drugz_minobs[ order( subset( drugz_minobs, select = normZ ), decreasing = T ), ]
  drugz_minobs['pval_supp'] = pnorm( (drugz_minobs$normZ), lower.tail = F )
  drugz_minobs['rank_supp'] = seq( from = 1, to = numGenes, by = 1)
  drugz_minobs['fdr_supp']  = ( drugz_minobs$pval_supp*numGenes ) / ( drugz_minobs$rank_supp )
  drugz_minobs = drugz_minobs[ order( subset( drugz_minobs, select = normZ ), decreasing = F ), ]
  # Formatting numeric data
  setDT(drugz_minobs)
  drugz_minobs[ , sumZ := round( sumZ, digits = 2) ]
  drugz_minobs[ , numObs := round( numObs, digits = 0) ]
  drugz_minobs[ , normZ := round( normZ, digits = 2) ]
  drugz_minobs[ , pval_synth := scientific( pval_synth, digits = 3)]
  drugz_minobs[ , rank_synth := round( rank_synth, digits = 0)]
  drugz_minobs[ , fdr_synth := scientific( fdr_synth, digits = 3)]
  drugz_minobs[ , pval_supp := scientific( pval_supp, digits = 3)]
  drugz_minobs[ , rank_supp := round( rank_supp, digits = 0)]
  drugz_minobs[ , fdr_supp := scientific( fdr_supp, digits = 3)]
  #
  # write output file
  #
  write.table( x = drugz_minobs, file = drugz_outfile , append = F, quote = F, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = TRUE )
  log_( 'Finished writing output file to disk' )
}
