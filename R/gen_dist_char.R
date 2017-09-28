#------------------------------------#
# Functions for general distribution #
# characteristics                    #
#------------------------------------#

# Index
# Lookup - 01:  sort_input
# Lookup - 02:  extract_npfd

# Lookup - 01
sort_input = function( input, level ) {
  # Purpose:
  # A function that sorts response times based on
  # magnitude and the choice values.
  # Arguments:
  # input - The output from the 'extract_var' function
  #         (defined in 'utility.R')
  # level - The specific level of the grouping factor
  #         (else NULL)
  # Returns:
  # A list consisting of...
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
  # at      = A sorted dataframe witht the
  #           response times, choice values, and
  #           current grouping levels

  at = input$ad

  select_sort = is.null( input$n_cv )
  if ( !select_sort ) {
    if ( input$n_cv > 1 ) select_sort = F
  }

  # If no variable for covariates was included


  if ( select_sort ) {
    # Sort response times, conditions on choice value
    sta = 1; end = 0;
    for ( vl in input$val ) {

      sel = input$ad$v == vl
      if ( any( sel ) ) {
        end = sum(sel) + end
        ind = sta:end
        tmp = by( input$ad$t[sel],
                  list( input$ad$g[sel] ), function(x) sort(x) )
        at$t[ ind ] = as.vector( unlist( tmp ) )
        at$v[ ind ] = vl
        # Extract group labels
        n = sapply( tmp, length )
        dn = dimnames( tmp )
        gl = data.frame( dn = dn, n = n )
        at$g[ ind ] = as.vector( unlist( apply( gl, 1, function(x) rep(x[1],x[2]) ) ) )
        sta = 1+end
      }

    }
  } else {

    # Sort response times, conditions on choice value
    sta = 1; end = 0;
    for ( vl in input$val ) {

      sel = input$ad$v == vl
      if ( any( sel ) ) {
        end = sum(sel) + end
        ind = sta:end
        tmp = by( input$ad$t[sel],
                  list( input$ad$g[sel], input$ad$cv[sel] ), function(x) sort(x) )
        at$t[ ind ] = as.vector( unlist( tmp ) )
        at$v[ ind ] = vl
        # Extract group labels
        dn = dimnames( tmp )
        gl = expand.grid( dn[[1]], dn[[2]] )
        gl$n = sapply( tmp, length )
        at$g[ ind ] = as.vector( unlist( apply( gl, 1, function(x) rep(x[1],x[3] ) ) ) )
        at$cv[ ind ] = as.vector( unlist( apply( gl, 1, function(x) rep(x[2],x[3] ) ) ) )
        sta = 1+end
      }

    }

  }

  # If a specific level for the grouping factor is given
  if ( !is.null( level ) ) {
    # If there is more than just one level
    if ( input$n_g > 1 ) {
      # If the given level is present
      if ( level %in% input$g_v ) {
        at = at[ at$g == level, ]
        input$ad = at
        input$n_g = 1
        input$g_v = level
        input$N = nrow( at )
      } else {
        err = 'Given level does not exist in grouping factor'
        stop( err, call. = FALSE )
      }
    }
  }
  if ( is.null( level ) | input$n_g == 1 ) {
    input$ad = at
  }

  return( input )
}

# Lookup - 02
extract_npfd = function(input) {
  # Purpose:
  # Creates a data frame with the number of observations,
  # proportions, and frequencies for each grouping factor.
  # Arguments:
  # input - The output from the 'extract_var' function
  #         (defined in 'utility.R')
  # Returns:
  # A data frame with...
  # G = The grouping factor levels
  # N = The number of observations
  # P = The proportion of times each choice value was
  #     chosen
  # F = The frequency with which each choice value was
  #     chosen

  select_sort = is.null( input$n_cv )
  if ( !select_sort ) {
    if ( input$n_cv > 1 ) select_sort = F
  }

  if ( select_sort ) {

    if ( input$n_g > 1 ) {
      # Determine the number of observations per group factor
      npfd = data.frame( G = input$g_v, N = rep( NA, input$n_g ) )

      # Compute the number of observations per grouping factor
      I = rep( 1, input$N )
      tmp = by( I, list( input$ad$g ), sum )
      npfd$N = as.vector( unlist( tmp ) )

      # Determine the proportions/frequencies per value and
      # group factor
      tmp = by( input$ad$v, list( input$ad$g ),
                function(x) p_cat(x,as.character(input$val)) )
      npfd$G = names( tmp )
      mat = matrix( unlist( tmp ), input$n_g, input$n_v, byrow = T )
      colnames(mat) = paste( 'P',
                             as.character( input$val ), sep = '' )
      npfd$P = mat
      mat = round( mat * npfd$N )
      colnames( mat ) = paste( 'F',
                               as.character( input$val ), sep = '' )
      npfd$F = round( mat )
    } else {
      npfd = data.frame( G = input$g_v,
                         N = input$N )
      npfd$P = matrix( p_cat( input$ad$v, input$val ), 1,
                       input$n_v )
      colnames( npfd$P ) = paste( 'P',
                                  as.character( input$val ),
                                  sep = '' )
      npfd$F = round( npfd$P * npfd$N )
      colnames( npfd$F ) = paste( 'F',
                                  as.character( input$val ),
                                  sep = '' )

    }
  } else {

    # Compute the number of observations per grouping factor
    # and the additional covariate
    I = rep( 1, input$N )
    tmp = by( I, list( input$ad$g, input$ad$cv ), sum )
    npfd = data.frame( N = rep( NA, length( tmp ) ) )
    npfd$N = as.vector( unlist( tmp ) )

    # Extract labels
    dn = dimnames( tmp )
    gl = expand.grid( dn[[1]], dn[[2]] )
    gl$n = sapply( tmp, length )
    npfd$G = apply( gl, 1, function(x) rep(x[1],x[3] ) )
    npfd$CV = apply( gl, 1, function(x) rep(x[2],x[3] ) )

    # Determine the proportions/frequencies per value and
    # group factor
    tmp = by( input$ad$v, list( input$ad$g, input$ad$cv ),
              function(x) p_cat(x,as.character(input$val)) )
    mat = matrix( unlist( tmp ), nrow(npfd), input$n_v, byrow = T )
    colnames(mat) = paste( 'P',
                           as.character( input$val ), sep = '' )
    npfd$P = mat
    mat = round( mat * npfd$N )
    colnames( mat ) = paste( 'F',
                             as.character( input$val ), sep = '' )
    npfd$F = round( mat )

  }

  return( npfd )
}
