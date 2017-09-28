#-------------------#
# Utility functions #
#-------------------#

# Index
# Lookup - 01:  p_cat [tested]
# Lookup - 02:  ln_interp
# Lookup - 03:  multiple_options
# Lookup - 04:  extract_var
# Lookup - 05:  pow
# Lookup - 06:  which_max_min
# Lookup - 07:  lower_upper

# Lookup - 01
p_cat = function( x, val ) {
  # Purpose:
  # Computes the proportion of times each value in
  # 'val' occurs.
  # Arguments:
  # x   - A vector of observations
  # val - A character vector with the possible outcomes
  # Returns:
  # A vector of the proportions.

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

# Lookup - 02
ln_interp = function( x, yPts, xPts ) {
  # Purpose:
  # Calculates the linear interpolation for a y-axis point that
  # lies on the line between given a pair of x-axis and y-axis
  # points.
  # Arguments:
  # x    - The x-axis point that corresponds to the desired y-axis
  #        point
  # yPts - The pair of y-axis points that the point lies between
  # xPts - The pair of x-axis points that the point lies between
  # Returns:
  # The predicted y-axis point.

  # Determine
  b = diff(yPts)/diff(xPts)
  y = yPts[1] + b*(x-xPts[1])

  return( y )
}

# Lookup - 03
multiple_options = function( var_options, feasible_var ) {
  # Purpose:
  # A function that asks users to indicate which variables
  # out of a feasible set they would like to use.
  # Arguments:
  # var_options  - The set of variable names
  # feasible_var - A logical vector indicating which
  #                variables are feasible
  # Returns:
  # The variable chosen by the user.

  cat( 'Multiple feasible variables found.', '\n' )
  cat( 'Select variable:', '\n' )

  for ( i in 1:sum( feasible_var ) ) {
    string = paste( i, ') ',
                    var_options[ feasible_var ][i],
                    sep = '' )
    cat( string, '\n' )
  }

  input = NULL

  while( !is.numeric( input ) ) {

    input = readline( 'Enter number: ' )
    input = as.numeric( input )
    if ( is.numeric( input ) ) {
      if ( input < 1 | input > i ) input = NULL
    }

  }

  return( var_options[ feasible_var ][input] )
}

# Lookup - 04
extract_var = function( df, label, keep, type ) {
  # Purpose:
  # Extracts response times, choice/accuracy,
  # and grouping variables from a data frame.
  # Arguments:
  # df    - A data frame
  # label - A character vector with up to
  #         3 labels
  # keep  - A logical vector of matching length
  #         to the number of observations in
  #         df
  # type  - The type of function (e.g., CDF,
  #         PDF, etc.)
  # Returns:
  # A list, consisting of...
  # t       = A vector of response times
  # v       = A vector of choice values
  # g       = A vector of indices for the grouping
  #           factor
  # val     = The set of possible choices
  # n_v     = The total number of possible choices
  # g_v     = The set of unique group indices
  # n_g     = The total number of elements for the
  #           grouping factor
  # cv_v    = The levels for an optional additional
  #           covariate
  # n_cv    = The number of levels for the additional
  #           covariate
  # N       = The total number of observations
  # ad      = A dataframe with the response times,
  #           choice values, and current grouping
  #           levels

  vnames = names( df ) # Extract variable names

  # Define defaults for covariate variable
  cv = NULL
  cv_v = NULL
  n_cv = NULL

  # If no labels for RT and choice variables are
  # provided, attempt to locate variables via
  # common variable names
  if ( is.null( label ) ) {

    ### Extract RT ###

    # Define set of common variable names for RT
    common_rt_names = c( 'RT', 'rt', 'Time', 'time',
                         't' )

    # Attempt to find default variable for RT
    rt_label_sel = common_rt_names %in% vnames
    if ( sum( rt_label_sel ) == 1 ) {
      var_name = common_rt_names[ rt_label_sel ]
      t = df[ , var_name ]
    }

    # If multiple labels
    if ( sum( rt_label_sel ) > 1 ) {
      # Ask user to identify variable to use
      var_name = multiple_options( common_rt_names, rt_label_sel )
      t = df[ , var_name ]
    }

    # If no labels
    if ( sum( rt_label_sel ) == 0 ) {
      stop( 'No feasible default variable names for RT found',
            call. = FALSE )
    }

    ### Extract choice/accuracy ###

    # Define set of common variable names for choice
    common_ch_names = c( 'Choice', 'Ch', 'ch', 'choice',
                         'Accuracy', 'accuracy', 'ac', 'Ac',
                         'Correct', 'correct', 'Co', 'Answer',
                         'answer' )

    # Attempt to find default variable for choice
    ch_label_sel = common_ch_names %in% vnames
    if ( sum( ch_label_sel ) == 1 ) {
      var_name = common_ch_names[ ch_label_sel ]
      v = df[ , var_name ]
    }

    # If multiple labels
    if ( sum( ch_label_sel ) > 1 ) {
      # Ask user to identify variable to use
      var_name = multiple_options( common_ch_names, ch_label_sel )
      v = df[ , var_name ]
    }

    # If no labels
    if ( sum( ch_label_sel ) == 0 ) {
      v = rep( 1, length( t ) )
    }

    g = rep( 1, length( t ) )
  }

  # If a single label for RTs is provided
  if ( length( label ) == 1 ) {

    # If label is present
    if ( label %in% names( df ) ) {
      t = df[ , label ]
      v = rep( 1, length( label ) )
      g = rep( 1, length( t ) )
    } else {
      stop( 'Variable name not found in data frame',
            call. = FALSE )
    }

  }

  # If labels for RT and choice/accuracy are provided
  if ( length( label ) == 2 ) {

    if ( all( label %in% names( df ) ) ) {
      t = df[ , label[1] ]
      v = df[ , label[2] ]
      g = rep( 1, length( t ) )
    } else {
      stop( 'Variable name(s) not found in data frame',
            call. = FALSE )
    }

  }

  if ( length( label ) >= 3 & !check_for_pvt_type(type) ) {

    if ( length( label ) > 3 ) {
      warn = paste(
        'More than 3 variable names were provided.', '\n',
        'Only the first 3 will be used.' )
      warning( warn, call. = FALSE )
      label = label[1:3]
    }

    if ( all( label %in% names( df ) ) ) {
      t = df[ , label[1] ]
      v = df[ , label[2] ]
      g = df[ , label[3] ]
    } else {
      stop( 'Variable name(s) not found in data frame',
            call. = FALSE )
    }

  }

  is_pvt = check_for_pvt_type(type)
  if ( is_pvt & length( label ) < 4 )
    stop( 'Need a grouping factor and an additional covariate for PvT figures',
          call. = FALSE )

  if ( length( label ) >= 4 & is_pvt ) {

    if ( length( label ) > 4 ) {
      warn = paste(
        'More than 4 variable names were provided.', '\n',
        'Only the first 4 will be used.' )
      warning( warn, call. = FALSE )
      label = label[1:4]
    }

    if ( all( label %in% names( df ) ) ) {
      t = df[ , label[1] ]
      v = df[ , label[2] ]
      g = df[ , label[3] ]
      cv = df[ , label[4] ]
    } else {
      stop( 'Variable name(s) not found in data frame',
            call. = FALSE )
    }

  }

  if ( is.null( keep ) ) {
    keep = rep( T, length( t ) )
  }

  val = unique( v );
  n_v = length( val )
  cur_g = g[ keep ]
  if ( all( is.na( cur_g ) ) )
    cur_g[ is.na( cur_g ) ] = 1
  ad = data.frame( t = t[ keep ],
                   v = v[ keep ],
                   g = cur_g )
  g_v = unique( cur_g )
  n_g = length( g_v )
  if ( !is.null( cv ) ) {
    ad$cv = cv[keep]
    cv_v = unique(cv[keep])
    n_cv = length( cv_v )
  }

  output = list(
    ad = ad,
    val = val,
    n_v = n_v,
    g_v = g_v,
    n_g = n_g,
    cv_v = cv_v,
    n_cv = n_cv,
    N = nrow( ad ) )

  return( output )
}

# Lookup - 05
pow = function( x, a ) {
  # Purpose:
  # A function to raise a variable to a power.
  # Arguments:
  # x - A continuous value
  # a - The power to raise the value
  # Returns:
  # The result of the value raised to the desired power.

  return( x^a )
}

# Lookup - 06
which_max_min = function( x, max = T ) {
  # Purpose:
  # A function that extracts the index for the
  # largest or smallest value in a vector.
  # Arguments:
  # x   - A vector of values
  # max - Logical; if true, the largest value is
  #       extracted
  # Returns:
  # The position of the largest/smallest value in
  # the vector.

  if ( max ) {
    out = which( x == max( x ) )
  } else {
    out = which( x == min( x ) )
  }

  return( out )
}

# Lookup - 07
lower_upper = function( int, dat ) {
  # Purpose:
  # Computes the lower and upper bounds for a set
  # of data within a desired increment.
  # Arguments:
  # int - The increment to round down/up to for
  #       the limits
  # dat - A vector of data
  # Returns:
  # A vector with a lower and upper limit.

  ll = int*floor(min(dat)/int)
  ll = c(ll,int*ceiling(max(dat)/int))

  return( ll )
}
