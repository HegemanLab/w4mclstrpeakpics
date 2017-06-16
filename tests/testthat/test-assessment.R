# test w4mclassfilter::w4m_filter_by_sample_class

#require(base)
#require(testthat)
#require(w4mclassfilter)

read_data_frame <- function(file_path, kind_string, failure_action = print) {
  # ---
  # read in the data frame
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
    failure_action(my.env$msg)
    return ( FALSE )
  }
  return (my.env)
}

#' @import testthat w4mclassfilter
#' @export
test_that("assessment test",{
  # set up variables
  variableMetadata_in  <- "input_variableMetadata.tsv"
  sampleMetadata_in <- "input_sampleMetadata.tsv"
  dataMatrix_in <- "input_dataMatrix.tsv"
  assessment_pdf <- "output_assessment.pdf"
  assessment_out <- "output_assessment.tsv"
  assessment_exp <- "expected_assessment.tsv"
  sample_selector_value       <- "pool"
  sample_selector_column_name <- "sampleType"
  # test input files
  data_matrix_input_env <- read_data_frame(dataMatrix_in, "data matrix input")
  expect_true(data_matrix_input_env$success, info = "read data matrix input")
  rm(data_matrix_input_env)
  sample_metadata_input_env <- read_data_frame(sampleMetadata_in, "sample metadata input")
  expect_true(sample_metadata_input_env$success, info = "read sample metadata input")
  rm(sample_metadata_input_env)
  variable_metadata_input_env <- read_data_frame(variableMetadata_in, "variable metadata input")
  expect_true(variable_metadata_input_env$success, info = "read variable metadata input")
  rm(variable_metadata_input_env)
  # filter, impute, and write output
  assessment_result <- pool_peak_assessment(
    sample_selector_value       = sample_selector_value
  , sample_selector_column_name = sample_selector_column_name
  , sample_metadata_path        = sampleMetadata_in
  , variable_metadata_path      = variableMetadata_in
  , data_matrix_path            = dataMatrix_in
  , output_pdf                  = assessment_pdf                 
  , output_tsv                  = assessment_out                 
  , output_rdata                = ""            
  )
  
  # read actual output file
  assessment_output_env <- read_data_frame(assessment_out, "assessment output")
  expect_true(assessment_output_env$success, info = "read assessment output")
  # read expected output file
  assessment_expected_env <- read_data_frame(assessment_exp, "assessment expected")
  expect_true(assessment_expected_env$success, info = "read assessment expected")
  # compare actuals with expecteds
  expect_equivalent(assessment_output_env$data, assessment_expected_env$data, info = "validate assessment")
})

