#---------------#
# CDF functions #
#---------------#

# Index
# Lookup - 01:  j_ecdf_one_level
# Lookup - 02:  j_cdf_by_bins
# Lookup - 03:  create_grp_j_ecdf
# Lookup - 04:  create_cdf_output
# Lookup - 05:  check_for_cdf_type

# Lookup - 01
j_ecdf_one_level = function( input, npfd ) {
  # Purpose:
  # Computes the joint empirical cumulative distribution
  # function for data from a single level of the grouping
  # factor.
  # Arguments:
  # input - The output from the 'extract_var' function
  #         (defined in 'utility.R')
  # npfd  - The output from the 'extract_npfd' function
  #         (defined in 'gen_dist_char.R')
  # Returns:
  # A list consisting of...
  # pd = A dataframe with the x and y plotting values;
  # xm = Set to null;
  # ym = Set to null.

  pd = data.frame( x = rep( NA, input$N ),
                   y = rep( NA, input$N ),
                   v = input$ad$v )

  for ( vl in input$val ) {
    sel = input$ad$v == vl
    if ( length( sel ) > 0 ) {
      pd$x[sel] = input$ad$t[ sel ];
      pd$v[sel] = input$ad$v[ sel ];
      pd$y[sel] = (1:sum(sel))/nrow(pd)
    }
  }
  out = list( pd = pd, xm = NULL, ym = NULL )

  return( out )
}

# Lookup - 02
j_ecdf_by_bins = function( x, bins ) {
  # Purpose:
  # Computes the joint cumulative probability
  # given a vector of observations over a
  # specified set of values.
  # Arguments:
  # x    - The vector of observations
  # bins - The set of bins to compute the ECDF over
  # Returns:
  # The empirical CDF at each bin value.

  out = sapply( bins, function(t) sum( x <= t )/length(x) )

  return( out )
}

# Lookup - 03
j_ecdf_group = function( input, npfd,
                         T_B = mean, T_b = mean,
                         T_x = mean ) {
  # Purpose:
  # Computes the joint cumulative probability over
  # multiple levels of a grouping factor.
  # Arguments:
  # input - The output from the 'extract_var' function
  #         (defined in 'utility.R')
  # npfd  - The output from the 'extract_npfd' function
  #         (defined in 'gen_dist_char.R')
  # T_B   - Test statistic for aggregating number of
  #         bins
  # T_b   - Test statistic for smoothing bin interval
  # T_x   - Test statistic for aggregating observations
  # Returns:
  # A list with...
  # pd = The x and y plotting values;
  # ym = A matrix with the full set of y-axis
  #      values.

  # Initializ output
  out = list( pd = NULL,
              xm = NULL, ym = NULL )

  # Initialize output
  ym = c();
  pd = matrix( NA, max( npfd$N ), 3 )
  colnames( pd ) = c( 'x', 'y', 'v' )
  pd = as.data.frame( pd )

  sta = 1; end = 0
  for ( i in 1:( input$n_v ) ) {

    # Create empty list
    ym = c( ym, list( NULL ) )

    # Compute number of bins to use
    sel = input$ad$v == input$val[i]
    bins = create_grp_bins( input, i, T_B, T_b )

    # Create index
    end = end + length( bins )
    ind = sta:end
    pd$x[ ind ] = bins;
    pd$v[ ind ] = input$val[i]

    ecdf = by( input$ad$t[sel],
               list( input$ad$g[sel] ),
               function(x) j_ecdf_by_bins( x, bins ) )

    mat = matrix( unlist( ecdf ), length( ecdf ), length(bins),
                  byrow = T )
    g_sel = names( ecdf )
    g_sel = as.character( npfd$G ) %in% g_sel

    # Weight by choice proportion
    mat = mat * npfd$P[g_sel,i]

    grp_ecdf = apply( mat, 2, T_x )
    pd$y[ ind ] = grp_ecdf
    ym[[i]] = mat

    sta = 1 + end
  }

  pd = pd[ 1:end, ]

  names( ym ) = as.character( input$val )
  out$pd = pd;
  out$ym = ym

  return( out )
}

# Lookup - 04
create_cdf_output = function( input, npfd, ... ) {
  # Purpose:
  # Computes the empirical CDF based on the number
  # grouping factor levels
  # Arguments:
  # input - The output from the 'extract_var' function
  #         (defined in 'utility.R')
  # npfd  - The output from the 'extract_npfd' function
  #         (defined in 'gen_dist_char.R')
  # ...   - Optional variables for the j_ecdf functions
  # Returns:
  # The list of plotting elements.

  plt = NULL
  if ( input$n_g == 1 ) plt = j_ecdf_one_level( input, npfd )
  if ( input$n_g > 1 ) plt = j_ecdf_group( input, npfd, ... )

  return( plt )
}

# Lookup - 05
check_for_cdf_type = function( type ) {
  # Purpose:
  # Checks a string character against a list of
  # possible labels for selecting the CDF option.
  # Arguments:
  # type = The input string
  # Returns:
  # A logical value, equal to 'TRUE' if one of the
  # labels matches the input string.

  cdf_types = c( 'CDF', 'cdf', 'CdF', 'cDF', 'CDf',
                 'cdF', 'cDf', 'Cdf', 'distribution',
                 'DF', 'df', 'ecdf', 'ECDF',
                 'eCDF', 'EcDF','ECdF','ECDf',
                 'Ecdf', 'eCdf', 'ecDf',
                 'ecdF', 'ECdf', 'eCDf', 'ecDF',
                 'EcdF', 'EcDf', 'eCdF' )

  return( type %in% cdf_types )
}
