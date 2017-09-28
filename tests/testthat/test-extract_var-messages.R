library(rtplots)
context("Tests of extract_var error/warning messages")

# Load in data from package
data( 'priming_data' )
# Rename for easy access
d = priming_data
colnames( d ) = c( 'S', 'Ac', 'Ch', 'RT', 'Cnd' )

test_that("extract_var produces correct errors/warnings", {

  # Not a data frame
  expect_error( rtplots( "Not a data frame" ), "Input variable must be a data frame" )

  # No variables
  input = data.frame( A = 1, B = 2 )
  expect_error( rtplots( input ), 'No feasible default variable names for RT found' )

  # Incorrect variable name(s)
  expect_error( rtplots( d, 'Y' ), 'Variable name not found in data frame' )
  expect_error( rtplots( d, c('RT','Y') ) )

  # Too many variable names
  expect_warning( rtplots( d, c('RT','Ch','S','Cnd') ),
                  paste(
                    'More than 3 variable names were provided.', '\n',
                    'Only the first 3 will be used.' ) )

  # 4 variables are acceptable for type == 'PVT'
  expect_silent( rtplots( d, c('RT','Ch','S','Cnd'), type = 'PVT' ) )

  # At least four variables are required for type == 'PVT'
  expect_error( rtplots( d, c('RT','Ch'), type = 'PVT' ),
                'Need a grouping factor and an additional covariate for PvT figures' )
  # Too many variable names
  expect_warning( rtplots( d, c('RT','Ch','S','Cnd','Ac'), type = 'PVT' ),
                  paste(
                    'More than 4 variable names were provided.', '\n',
                    'Only the first 4 will be used.' ) )

})

# Clean up workspace
rm( d )

