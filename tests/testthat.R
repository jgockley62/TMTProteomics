library(testthat)
source('code/01_Analysis_Functions.R')
synapser::synLogin()
#test_check("code/01_Analysis_Functions.R")
test_file("tests/testthat.R")

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
  expect_identical( is.data.frame( TMT_Express_Load('syn21266454')), TRUE)
})