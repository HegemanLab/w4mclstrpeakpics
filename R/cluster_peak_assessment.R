#library(grDevices)
#library(reshape2)
#library(sqldf)


#' @title
#' W4M Cluster Peak Pics
#'
#' @description
#' Visualize feature-peaks among a cluster of samples identified through W4M XCMS preprocessing.
#'
#' @param sample_selector_column_name string input: column of W4M/XCMS sampleMetadata holding selector string values (default: "sampleType")
#' @param sample_selector_value       string input: value within selector column to identify samples for analysis (default: "pool")
#' @param sample_metadata_path        string input: path to W4M/XCMS sampleMetadata tab-separated values file 
#' @param variable_metadata_path      string input: path to W4M/XCMS variableMetadata tab-separated values file 
#' @param data_matrix_path            string input: path to W4M/XCMS dataMatrix tab-separated values file
#' @param output_pdf                  string output: path to write assessment figure PDF
#' @param output_tsv                  string output: path to write assessment summary tab-separated values file
#' @param output_rdata                string output: (optional) path to write RData containing all processing and plotting intermediates
#' @param failure_action              function(msg): action to take when a failure occurs (defalut: "print")
#'
#' @return environment: env containing 'df_result' that was written to output_tsv plus all intermediate and plotted variables (see source code)
#' 
#' @author Art Eschenlauer, \email{esch0041@@umn.edu}
#' @concept w4m workflow4metabolomics
#' @seealso \url{https://github.com/HegemanLab/w4mclstrpeakpics}
#' @seealso \url{http://workflow4metabolomics.org/}
#'
#' @importFrom grDevices dev.list dev.off pdf rgb
#' @importFrom graphics legend matplot par plot title
#' @importFrom utils str write.table
#'
#' @export
cluster_peak_assessment <- function(
  sample_selector_column_name = "sampleType"
, sample_selector_value       = "pool"
, sample_metadata_path  
, variable_metadata_path
, data_matrix_path      
, output_pdf            
, output_tsv            
, output_rdata                = ""            
, failure_action              = print
) {
  peak_assessment_env                        <- new.env()
  peak_assessment_env$sample_selector_value  <- sample_selector_value 
  peak_assessment_env$sample_selector        <- sample_selector_column_name
  peak_assessment_env$sample_metadata_path   <- sample_metadata_path  
  peak_assessment_env$variable_metadata_path <- variable_metadata_path
  peak_assessment_env$data_matrix_path       <- data_matrix_path      
  peak_assessment_env$output_pdf             <- output_pdf            
  peak_assessment_env$output_tsv             <- output_tsv            
  peak_assessment_env$output_rdata           <- output_rdata          

  compute_plots(peak_assessment_env, failure_action)

  # write results
  result <- peak_assessment_env$df_result
  prevalence <- result$prevalence
  write_env <- write_result(
    result = result
  , file_path = peak_assessment_env$output_tsv
  , kind_string = "output table"
  , failure_action
  )
  if (!is.environment(write_env)) {
    stop(sprintf("write_result failed, producing result %s",as.character(write_env)))
  }
  if ( write_env$success != TRUE ) {
    stop(write_env$msg)
  }

  result_env <- new.env()
  result_env$result <- FALSE
  tcr <- tryCatch.W.E({
    result <<- TRUE
    # save RData if requested
    if (nchar(peak_assessment_env$output_rdata) > 0) {
      save( peak_assessment_env, file = peak_assessment_env$output_rdata )
    }
    # capture plot and write to PDF; then close any devices opened in the process
    plot2pdf(
      file_name     = peak_assessment_env$output_pdf
    , plot_function = function() { draw_plots(peak_assessment_env) }
    , width         = 12
    , height        = 12
    )
    result_env$result <- TRUE
  })
  if (!result_env$result) {
    str(tcr)
  }

  return (result_env$result)
}

 # from 'demo(error.catching)'
tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w){
    # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(
    value = withCallingHandlers(
      tryCatch(expr, error = function(e) e)
    , warning = w.handler
    )
  , warning = W
  )
}


write_result <- function(result, file_path, kind_string, failure_action = print) {
  my.env <- new.env()
  my.env$success <- FALSE
  my.env$msg <- sprintf("no message writing %s", kind_string)
  tryCatch(
    expr = {
      write.table(
        x = result
      , sep = "\t"
      , file = file_path
      , quote = FALSE
      , row.names = FALSE
      )
      my.env$success <- TRUE
    }
  , error = function(e) {
     my.env$msg <- sprintf("%s write failed", kind_string)
    }
  )
  if (!my.env$success) {
    failure_action(my.env$msg)
    return ( FALSE )
  }
  return (my.env)
}

# read_data_frame - read a w4m data frame, with error handling
#   e.g., data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
read_data_frame <- function(file_path, kind_string, rdf_failure_action = failure_action) {
  my.env <- new.env()
  my.env$success <- FALSE
  my.env$msg <- sprintf("no message reading %s", kind_string)
  tryCatch(
    expr = {
      my.env$data    <- utils::read.delim( fill = FALSE, file = file_path )
      my.env$success <- TRUE
    }
  , error = function(e) {
     my.env$ msg <- sprintf("%s read failed", kind_string)
    }
  )
  if (!my.env$success) {
    rdf_failure_action(my.env$msg)
    return ( FALSE )
  }
  return (my.env)
}

compute_plots <- function(plot_env, failure_action = print) {

  # read data structures

  # read in the sample metadata
  smpl_metadata_input_env <- read_data_frame(plot_env$sample_metadata_path, "sample metadata input")
  if (!smpl_metadata_input_env$success) {
    failure_action(smpl_metadata_input_env$msg)
    return ( FALSE )
  }
  sampleMetadata <- smpl_metadata_input_env$data

  # read in the variable metadata
  vrbl_metadata_input_env <- read_data_frame(plot_env$variable_metadata_path, "variable metadata input")
  if (!vrbl_metadata_input_env$success) {
    failure_action(vrbl_metadata_input_env$msg)
    return ( FALSE )
  }
  variableMetadata <- vrbl_metadata_input_env$data
  
  # read in the data matrix
  #   dataMatrix <- read_data_frame(
  #     file_path = plot_env$data_matrix_path
  #   )
  data_matrix_input_env <- read_data_frame(plot_env$data_matrix_path, "data matrix input")
  if (!data_matrix_input_env$success) {
    failure_action(data_matrix_input_env$msg)
    return ( FALSE )
  }
  dataMatrix <- data_matrix_input_env$data
  
  # identify names of pooled samples

  selected_samples <- sampleMetadata[sampleMetadata[plot_env$sample_selector] == plot_env$sample_selector_value,1]


  # extract matrix of intensities for only pooled samples

  m <- as.matrix(dataMatrix[,colnames(dataMatrix) %in% selected_samples])
  rownames(m) <- unname(unlist(dataMatrix[,1]))


  # compute and plot "how many samples have count" versus "count of samples with feature"

  feature_intensity_counts <- sapply( X = rownames(m), FUN = function(r) sum(m[r,]>0) )
  max_counts <- max(feature_intensity_counts)
  seq_1_to_maxcounts <- 1:max_counts
  how_many_features_have_count <- sapply(X = seq_1_to_maxcounts, FUN = function(x) sum(feature_intensity_counts == x))
  how_likely_are_features <- ( how_many_features_have_count * seq_1_to_maxcounts ) / max_counts
  proportionate_likelihood <- how_likely_are_features / sum(how_likely_are_features)

  pch    <- c(1,17)
  my_col <- c("black","olivedrab")
  plot_env$l_matplot_1 <- list(
    x = seq_1_to_maxcounts
  , y = data.frame (
          how_many   = how_many_features_have_count
        , likelihood = how_likely_are_features
        )
  , pch  = pch
  , ylim = c( 0, max(1.1 * max(how_many_features_have_count)) )
  , xlab = "Count of Samples having Feature"
  , ylab = "Number of, or Relative Likelihood of, Features having Count"
  , main = "Feature Number and Likelihood"
  , sub = "(Prevalence)"
  , col = my_col
  )
  plot_env$l_legend_1 <- list(
    x = 4
  , y = 1.1 * max(how_many_features_have_count)
  , legend = c("Number of features", "Likelihood of features")
  , pch  = pch
  , col = my_col
  , cex = 0.75
  )

  plot_env$df_result <- data.frame(
    prevalence = seq_1_to_maxcounts,
    number     = how_many_features_have_count,
    likelihood = how_likely_are_features,
    proportion = proportionate_likelihood
  )
  # sorting by decreasing prevalence is straightforward with SQL but obtuse with pure R
  #   with( plot_env, df_result <- sqldf::sqldf("select * from df_result order by prevalence desc") )
  # compare to the shorter but almost illegible R:
  with( plot_env, df_result <- df_result[ seq( dim(df_result)[1], 1 ), ] )

  # compute and plot intensity of each feature for each sample versus prevalence of the feature among the samples

  lut_count_to_feature_names <- lapply(
    X = seq_1_to_maxcounts
  , FUN = function(i) names(feature_intensity_counts[feature_intensity_counts == i])
  )

  lut_count_to_feature_intensities <- sapply(
    X = seq_1_to_maxcounts
  , FUN = function(x) { tmp <-m[lut_count_to_feature_names[[x]],]; tmp[tmp>0] }
  )

  tmp <- sapply(X = seq_1_to_maxcounts, FUN = function(i) rep(i, length(lut_count_to_feature_intensities[[i]])) )
  tmp <- unlist(tmp)


  plot_env$l_plot_2 <- list(
    x = jitter(unlist(tmp))
  , y = log10(unlist(lut_count_to_feature_intensities))
  , xlab = "Count of Samples having Feature"
  , ylab = "Intensity of Peak"
  , main = "Peak Intensity"
  , sub = "(Prevalence)"
  , col = rgb(0, 0, 0, 0.2)
  , pch = 16
  )

  ############

  mm <- reshape2::melt( m, varnames = c("feature", "sample"), value.name="intensity" )
  mm <- mm[ mm$intensity > 0, ]

  mzmargin <- with( variableMetadata, max(mz) - min(mz) ) / 20.0
  mzlim    <- with( variableMetadata, c( min(mz) - mzmargin, max(mz) + mzmargin ) )
  rtmargin <- with( variableMetadata, max(rt) - min(rt) ) / 20.0
  rtlim    <- with( variableMetadata, c( min(rt) - rtmargin, max(rt) + rtmargin ) )

  # create data.frame(mz, rt, intensity, count)
  heatplotdata <- sqldf::sqldf("
    select mz,
           rt,
           intensity,
           counts.sample_count as `count`
    from   variableMetadata, 
           mm,
           (select feature, count(intensity) as sample_count from mm group by feature) as counts
    where variableMetadata.variableMetadata = mm.feature
      and variableMetadata.variableMetadata = counts.feature
    order by mz, rt
  ")

  # diameter reflects intensity, hue reflects prevalence
  with(
    heatplotdata
  , {
      cexmaxsq <- 4.0 / max(mm$intensity)
      alpha    <- 1.0  / max_counts
      plot_env$l_plot_3 <<- list(
        x = rt
      , y = mz
      , xlab = "retention time"
      , ylab = "m/z"
      , xlim = c(min(rt), min(rt) + 1.1 * (max(rt) - min(rt)))
      , cex = sqrt(intensity * cexmaxsq)
      , col = grDevices::hcl(
                h = 180 * (1 + count / max_counts)
              , alpha = alpha
              , l = 85
              , c = 400
              )
      , pch = 16
      , main = "Symbol area/intensity reflect ion intensity"
      , sub = "Symbol hue subtly reflects prevalence"
      )
      lcol <- grDevices::hcl(
                h = 180 * (1 + 1:10 / max_counts)
              , alpha = 1.0 
              , l = 85
              , c = 400
              )
      plot_env$l_legend_3 <<- list(
        x = max(rt)
      , y = max(mz)
      , legend = 10:1
      , pch = 16
      , cex = 0.75
      , col = lcol[10:1]
      )
    }
  )

  # symbol size/shape reflects prevalence, color vividness reflects intensity
  with(
      sqldf::sqldf("
        select mz, rt, count , avg(intensity)*count as total_intensity
        from heatplotdata 
        group by mz, rt, count
      ")
      , {
          cexscalesq <- 2.0 / max(max_counts)
          # note the use of pmin, operates properly on a vector
          pch <- 18 - pmin(3, as.integer(max_counts - 1:max_counts))
          scale <- unlist(list(1,3.5/5,2.75/3,3.5/4)[19 - pch])
          cex <- 1:max_counts * cexscalesq * scale
          intensity <- total_intensity / max(total_intensity)
          plot_env$l_plot_4 <<- list(
            x = rt
          , y = mz
          , xlab = "retention time"
          , ylab = "m/z"
          , xlim = c(min(rt), min(rt) + 1.1 * (max(rt) - min(rt)))
          , cex = cex[count]
          , col = grDevices::hcl(
                    h = 180 * (1 + intensity / max(intensity))
                  , alpha = 0.5
                  , l = 85
                  , c = 35 + 800 * intensity / max(intensity)
                  )
          , pch = pch[count]
          , main = "Symbol size/shape reflects prevalence"
          , sub = "Symbol color vividness subtly reflects intensity"
          )
          plot_env$l_legend_4 <<- list(
            x = max(rt)
          , y = max(mz)
          , legend = 10:1
          , pch = pch[10:1]
          , pt.cex = cex[10:1]
          , cex=0.75
          )
      }
  )
}

dev.off.all <- function() {
  while (!is.null(dev.list())) { dev.off() }
}

# capture plot and write to PDF; then close any devices opened in the process
plot2pdf <- function(
  file_name
, plot_function
, width = 12
, height = 12
) {
  cur.dev <- dev.list()
  filename <- file_name
  pdf(file = filename, width = width, height = height)
  plot_function()
  dev.off()
  if (is.null(cur.dev)) {
      dev.off.all()
  } else {
      while ( length(dev.list()) > length(cur.dev) ) { dev.off() }
  }
}

draw_plots <- function(plot_env) {
  # show results as a 2 x 2 matrix of plots

  par(mfcol=c(2,2))

  matplot(
    x    = plot_env$l_matplot_1$x
  , y    = plot_env$l_matplot_1$y
  , pch  = plot_env$l_matplot_1$pch
  , col  = plot_env$l_matplot_1$col
  , ylim = plot_env$l_matplot_1$ylim
  , xlab = ""
  , ylab = ""
  , main = ""
  , sub  = ""
  )
  title(
    xlab = plot_env$l_matplot_1$xlab
  , ylab = plot_env$l_matplot_1$ylab
  , main = plot_env$l_matplot_1$main
  , sub  = plot_env$l_matplot_1$sub 
  , mgp  = c(2,1,0)
  , cex.lab = 1.2
  , cex.sub = 1.2
  , cex.main = 1.2
  )
  legend(
    x      = plot_env$l_legend_1$x
  , y      = plot_env$l_legend_1$y
  , legend = plot_env$l_legend_1$legend
  , pch    = plot_env$l_legend_1$pch
  , col    = plot_env$l_legend_1$col
  , cex    = plot_env$l_legend_1$cex
  )

  plot(
    x    = plot_env$l_plot_2$x   
  , y    = plot_env$l_plot_2$y   
  , xlab = ""
  , ylab = ""
  , main = ""
  , sub  = ""
  , col  = plot_env$l_plot_2$col 
  , pch  = plot_env$l_plot_2$pch 
  )
  title(
    xlab = plot_env$l_plot_2$xlab
  , ylab = plot_env$l_plot_2$ylab
  , main = plot_env$l_plot_2$main
  , sub  = plot_env$l_plot_2$sub 
  , mgp  = c(2,1,0)
  , cex.lab = 1.2
  , cex.sub = 1.2
  , cex.main = 1.2
  )

  plot(
    x    = plot_env$l_plot_3$x
  , y    = plot_env$l_plot_3$y
  , xlim = plot_env$l_plot_3$xlim
  , cex  = plot_env$l_plot_3$cex
  , col  = plot_env$l_plot_3$col
  , pch  = plot_env$l_plot_3$pch
  , xlab = ""
  , ylab = ""
  , main = ""
  , sub  = ""
  )
  title(
    xlab = plot_env$l_plot_3$xlab
  , ylab = plot_env$l_plot_3$ylab
  , main = plot_env$l_plot_3$main
  , sub  = plot_env$l_plot_3$sub
  , mgp  = c(2,1,0)
  , cex.lab = 1.2
  , cex.sub = 1.2
  , cex.main = 1.2
  )
  legend(
    x      = plot_env$l_legend_3$x
  , y      = plot_env$l_legend_3$y
  , legend = plot_env$l_legend_3$legend
  , pch    = plot_env$l_legend_3$pch
  , cex    = plot_env$l_legend_3$cex
  , col    = plot_env$l_legend_3$col
  )

  plot(
    x    = plot_env$l_plot_4$x
  , y    = plot_env$l_plot_4$y   
  , xlim = plot_env$l_plot_4$xlim
  , cex  = plot_env$l_plot_4$cex 
  , col  = plot_env$l_plot_4$col 
  , pch  = plot_env$l_plot_4$pch 
  , xlab = ""
  , ylab = ""
  , main = ""
  , sub  = ""
  )
  title(
    xlab = plot_env$l_plot_4$xlab
  , ylab = plot_env$l_plot_4$ylab
  , main = plot_env$l_plot_4$main
  , sub  = plot_env$l_plot_4$sub 
  , mgp  = c(2,1,0)
  , cex.lab = 1.2
  , cex.sub = 1.2
  , cex.main = 1.2
  )
  legend(
    x      = plot_env$l_legend_4$x     
  , y      = plot_env$l_legend_4$y     
  , legend = plot_env$l_legend_4$legend
  , pch    = plot_env$l_legend_4$pch   
  , pt.cex = plot_env$l_legend_4$pt.cex
  , cex    = plot_env$l_legend_4$cex   
  )

}

smoke_test_pool_peak_assessment <- function() {
  smoke_test_result <<- pool_peak_assessment(
    sample_selector_value       = "pool"
  , sample_selector_column_name = "sampleType"
  , sample_metadata_path        = "~/src/w4mclstrpeakpics/tests/testthat/input_sampleMetadata.tsv"
  , variable_metadata_path      = "~/src/w4mclstrpeakpics/tests/testthat/input_variableMetadata.tsv"
  , data_matrix_path            = "~/src/w4mclstrpeakpics/tests/testthat/input_dataMatrix.tsv"
  , output_pdf                  = "~/src/w4mclstrpeakpics/tests/testthat/output_assessment.pdf"
  , output_tsv                  = "~/src/w4mclstrpeakpics/tests/testthat/output_assessment.tsv"
  , output_rdata                = "~/src/w4mclstrpeakpics/tests/testthat/output_assessment.RData"
  )
  print("results stored in environment 'smoke_test_result'")
  return(smoke_test_result)
}

# print("smoke test with 'smoke_test_pool_peak_assessment()'")
