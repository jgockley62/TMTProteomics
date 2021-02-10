library(testthat)
if (!require(covr)) install.packages('covr')
library(covr)
#source('~/Desktop/Projects/TMT_Proteomics/code/01_Analysis_Functions.R')
#source('../R/01_Analysis_Functions.R')

library(testthat)
library(TMTProteomics)

# synapser::synLogin()

# Check TMT_Load()
## Input Tests
test_that("TMT_Express_Load checks input for type", {
  expect_identical(TMT_Express_Load(4,1), "Invalid Synapse ID: Input SynID is not a character string")
  expect_identical(TMT_Express_Load('NotAString',1), "Invalid Synapse ID: Input SynID is not a valid synapse ID")
  expect_identical(TMT_Express_Load('syn21266',1), "Invalid Synapse ID: Input SynID is not a valid synapse ID")
  #expect_identical(TMT_Express_Load('syn88888888',1), "Syn ID does not exist")
  expect_identical(TMT_Express_Load('syn2126645458',1), "Invalid Synapse ID: Input SynID is not a valid synapse ID")
  expect_identical(TMT_Express_Load('syn21641411',1), "User does not have READ access")
  expect_identical(TMT_Express_Load('syn21442609',1), "File to import must be a CSV File")
  expect_identical(TMT_Express_Load('syn21266454',"a"), "names is not a number")
})
## Output Tests
test_that("TMT_Express_Load checks output for type", {  
  expect_identical( is.data.frame( TMT_Express_Load('syn21266454', 1)), TRUE)
  expect_identical( is.data.frame( TMT_Express_Load('syn21266454', 0)), TRUE)
})

##covr <- file_coverage("code/01_Analysis_Functions.R", "tests/testthat.R")
##covr <- file_coverage("code/01_Analysis_Functions.R", "tests/test_01_Analysis_Functions.R")

