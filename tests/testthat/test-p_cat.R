context("Tests for p_cat function")

test_that("p_cat returns correct output", {

  # Define internal function
  p_cat = function( x, val ) {

    # Initialize output
    out = numeric( length( val ) )
    names( out ) = val

    # Tally frequencies for each choice
    freq = table( x )
    # Compute associated proportions
    prp = freq/sum( freq )

    # Extract values
    out[ names( prp ) ] = as.numeric( prp )

    return( out )
  }

  # Generate input
  set.seed( 387 )
  x = sample( 0:1, 10, replace = T )
  output = as.numeric( table(x)/length(x) )
  names( output ) = as.character( 0:1 )

  # Check with 2 values
  expect_equal( output, p_cat( x, 0:1 ) )

  # Check with 3 values
  output['2'] = 0
  expect_equal( output, p_cat( x, 0:2 ) )

  # Generate string input
  set.seed( 930 )
  val = c("old","new","unknown")
  x = sample( val[1:2], 10, replace = T )
  output = as.numeric( table(x)[val[1:2]]/length(x) )
  names( output ) = val[1:2]

  # Check with 2 values
  expect_equal( output, p_cat( x, val[1:2] ) )

  # Check with 3 values
  output['unknown'] = 0
  expect_equal( output, p_cat( x, val[1:3] ) )

})
