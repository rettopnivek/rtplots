#----------------------------------#
# Functions for determining number #
# of histogram bins                #
#----------------------------------#

# Index
# Lookup - 01:  bin_create
# Lookup - 02:  bin_cost_function
# Lookup - 03:  n_bins
# Lookup - 04:  create_grp_bins

# Lookup - 01
bin_create = function( N, x ) {
  # Purpose:
  # Creates a set of equally spaced bins and
  # computes the number of events that fall
  # within these bins.
  # Arguments:
  # N - The number of bins
  # x - A vector of data
  # Returns:
  # A list with the N + 1 intervals defining the
  # bins and the frequencies of the events falling
  # within those intervals.

  # The range of the observations
  R = range( x )

  # Generate the bins
  b = seq( R[1], R[2], length = N + 1 )

  # Compute the number of events that
  # fall within a bin
  ki = sapply( b, function(b) sum( x < b ) )[-1]
  ki = c( ki[1], diff( ki ) )

  return( list( b = b, ki = ki ) )
}

# Lookup - 02
bin_cost_function = function( N, x ) {
  # Purpose:
  # Computes a loss function that can be
  # minimized in order to determine the
  # optimal number of bins in a histogram.
  # Arguments:
  # N - The number of bins
  # x - A vector of values
  # Reference:
  # Shimazaki, H., & Shinomoto, S. (2007). A method for
  # selecting the bin size of a time histogram. Neural
  # computation, 19(6), 1503-1527.
  # Returns:
  # The value of the loss function.

  tmp = bin_create( N, x )
  b = tmp$b; ki = tmp$ki; rm(tmp);

  delta = b[2] - b[1]

  # Compute mean and variance
  # in fashion to avoid overflow
  K = 0 # Sample size
  M = 0 # Mean
  M2 = 0 # 2nd moment

  # Loop over observations
  for ( i in 1:length( ki ) ) {
    K1 = K; K = K + 1;
    term1 = ki[i] - M;
    term2 = term1 / K;
    term3 = term1 * term2 * K1;
    # Update mean
    M = M + term2;
    # Update second moment
    M2 = M2 + term3;
  }

  # Compute biased variance
  V = M2/K;

  C_delta = ( 2*M - V )/pow( delta, 2 )
  return( C_delta )
}

# Lookup - 03
n_bins = function(x) {
  # Purpose:
  # Computes the optimal number of bins for the
  # histogram summarizing a vector of data.
  # Arguments:
  # x - A vector of values
  # Returns:
  # The optimal number of bins for the histogram.

  mx = length( x )

  if ( mx > 1 ) {
    res = optimize( bin_cost_function, c( 1, mx ), x = x )
    out = round( res$minimum )
  } else out = 1

  return( out )
}

# Lookup - 04
create_grp_bins = function( input, vt,
                            T_B = mean, T_b = mean ) {
  # Purpose:
  # A function to create a set of intervals for the
  # bins used in histograms/PDFs/CDFs smoothed
  # over a grouping factor.
  # Arguments:
  # input - A list with the data frame of observations,
  #         as well as their characteristics (e.g.,
  #         the number, unique values, etc...)
  # vt    - The index for the particular choice value
  # T_B   - The type of univariate test statistic to compute
  #         when determining the number of bins to use over
  #         the grouping factor
  # T_b   - The type of univariate test statistic to compute
  #         when smoothing the bin intervals over the
  #         grouping factor
  # Returns:
  # A vector containing the intervals for the bins, smoothed
  # over the grouping factor.

  # Select response type
  sel = input$ad$v == input$val[vt]

  # Compute number of smoothing bins needed for
  # each subject
  nb = by( input$ad$t[sel],
           list( input$ad$g[sel] ),
           n_bins )

  # Compute number of bins for group
  B = round( T_B( as.vector( nb ) ) )

  # Compute the sequence of bins for each subject
  bins = by( input$ad$t[sel],
             list( input$ad$g[sel] ),
             function(x) bin_create( B, x )[[1]] )
  # Determine the bins for the group
  mat = matrix( unlist( bins ), length(bins), B + 1, byrow = T )
  grp_bin = apply( mat, 2, T_b )

  return( grp_bin )
}

