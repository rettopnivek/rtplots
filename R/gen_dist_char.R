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
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...

  at = input$ad

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
      at$g[ ind ] = input$ad$g[ sel ]
      sta = 1+end
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
  # input - ...
  # Returns:
  # A data frame with...
  # G = ...
  # N = ...
  # P = ...
  # F = ...

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
    mat = matrix( unlist( tmp ), input$n_g, 2, byrow = T )
    colnames(mat) = paste( 'P',
                           as.character( input$val ), sep = '' )
    npfd$P = mat
    mat = round( mat * npfd$N )
    colnames( mat ) = paste( 'F',
                             as.character( input$val ), sep = '' )
    npfd$F = mat
  } else {
    npfd = data.frame( G = input$g_v,
                       N = input$N )
    npfd$P = matrix( p_cat( input$ad$v, input$val ), 1,
                     input$n_v )
    colnames( npfd$P ) = paste( 'P',
                                as.character( input$val ),
                                sep = '' )
    npfd$F = npfd$P * npfd$N
    colnames( npfd$F ) = paste( 'F',
                                as.character( input$val ),
                                sep = '' )

  }

  return( npfd )
}
