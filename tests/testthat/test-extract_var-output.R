library(rtplots)
context("Tests for extract_var function output")

# Load in data from package
data( 'priming_data' )
# Rename for easy access
d = priming_data
colnames( d ) = c( 'S', 'Ac', 'Ch', 'RT', 'Cnd' )

test_that("extract_var correctly extracts variables", {

  # Create simple test input
  input = data.frame( RT = rep( seq( .1, .9, .2 ), 2 ),
                      Ch = rep( 0:1, each = 5 ),
                      S = 1,
                      Cnd = 1 )

  # Correct output (Default lookup)
  output = data.frame( t = input$RT, v = input$Ch )
  output$g = rep( "1", nrow( output ) )
  expect_equal( output, rtplots( input )$ad )

  # obj = rtplots( input, 'RT' )

  #
  # rtplots( input, c( 'RT', 'Ch' ) )

})
