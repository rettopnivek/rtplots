#---------------#
# PvT functions #
#---------------#

# * Incomplete documentation

# Index
# Lookup - 01:  pvt_one_level*
# Lookup - 02:  pvt_multi_levels*
# Lookup - 03:  create_pvt_output
# Lookup - 04:  check_for_pvt_type
# Lookup - 05:  pvt_default_options*

# Lookup - 01
pvt_one_level = function( input, npfd, f = NULL ) {
  # Purpose:
  # ...
  # Arguments:
  # input - ?
  # npfd  - ?
  # f     - ?
  # Returns:
  # ...

  # If no function for the test statistic is provided
  # use mean as default
  if ( is.null( f ) ) f = mean

  # Create output dataframe
  l = nrow( npfd ); n_v = input$n_v
  pd = data.frame( x = as.vector( npfd$P ),
                   v = rep( input$val, each = l ),
                   cv = rep( npfd$CV, n_v ) )

  inc = 1
  vrb = rep( ' ', nrow(pd) )
  for ( vl in input$val ) {
    sel = input$ad$v == vl
    if ( sum( sel ) > 0 ) {

      tmp = by( input$ad$t[sel], list( input$ad$cv[sel] ), f )
      if ( inc == 1 ) mat =  matrix( NA, nrow( pd ), length(tmp[[1]]) )

      sel = pd$v == vl & pd$cv %in% names( tmp )
      vrb[sel] = names( tmp )
      mat[sel,]  = matrix( unlist( tmp ), length( tmp ), length(tmp[[1]]),
                    byrow = T )
      inc = inc + 1
    }
  }
  if ( exists( 'mat' ) ) pd$y = mat

  out = list( pd = pd, xm = NULL, ym = NULL )

  return( out )
}


# Lookup - 02
pvt_multi_levels = function( input, npfd, f = NULL, T_x = NULL ) {
  # Purpose:
  # ...
  # Arguments:
  # input
  # npfd
  # T_x
  # Returns:
  # ...

  if ( is.null(f) ) f = mean
  if ( is.null(T_x) ) T_x = function(x) mean(x,na.rm = T)

  n_v = input$n_v; n_cv = input$n_cv; l = n_v * n_cv
  pd = data.frame( x = rep( NA, l ),
                   v = rep( input$val, each = n_cv ),
                   cv = rep( unique( npfd$CV ), n_v ) )

  n_t = length( f( input$ad$t ) )

  xm = matrix( NA, input$n_g, nrow( pd ) )
  ym = array( NA, dim = c( input$n_g, nrow( pd ), n_t ) )

  for ( g in 1:input$n_g ) {

    sel_1 = input$ad$g == input$g_v[g]
    sel_2 = npfd$G == input$g_v[g]

    inc = 1
    vrb = rep( ' ', nrow(pd) )
    for ( vl in input$val ) {
      sel = input$ad$v == vl
      if ( sum( sel ) > 0 ) {

        tmp = by( input$ad$t[sel_1 & sel], list( input$ad$cv[sel_1 & sel] ), f )
        if ( inc == 1 ) mat =  matrix( NA, nrow( pd ), length(tmp[[1]]) )

        sel = pd$v == vl & pd$cv %in% names( tmp )
        vrb[sel] = names( tmp )
        mat[sel,]  = matrix( unlist( tmp ), length( tmp ), length(tmp[[1]]),
                             byrow = T )
        inc = inc + 1
      }
    }
    if ( exists( 'mat' ) ) ym[g,,] = mat
    xm[g,] = as.vector( npfd$P[sel_2,] )

  }

  pd$x = apply( xm, 2, T_x )
  pd$y = apply( ym, 2:3, T_x )

  out = list( pd = pd, xm = xm, ym = ym )

  return( out )
}

# Lookup - 03
create_pvt_output = function( input, npfd, ... ) {
  # Purpose:
  # Computes the ? function based on the number
  # grouping factor levels
  # Arguments:
  # input - The output from the 'extract_var' function
  #         (defined in 'utility.R')
  # npfd  - The output from the 'extract_npfd' function
  #         (defined in 'gen_dist_char.R')
  # ...   - Optional variables for the pvt functions
  # Returns:
  # The list of plotting elements.

  plt = NULL
  if ( input$n_g == 1 ) plt = pvt_one_level( input, npfd, ... )
  if ( input$n_g > 1 ) plt = pvt_multi_levels( input, npfd, ... )

  return( plt )
}

# Lookup - 04
check_for_pvt_type = function( type ) {
  # Purpose:
  # Checks a string character against a list of
  # possible labels for selecting the PvT option.
  # Arguments:
  # type = The input string
  # Returns:
  # A logical value, equal to 'TRUE' if one of the
  # labels matches the input string.

  pvt_types = c( 'PVT', 'pvt', 'PvT', 'pVT', 'PVt',
                 'pvT', 'pVt', 'Pvt' )

  return( type %in% pvt_types )
}

# Lookup - 05
pvt_default_options = function( opt, pd, ver ) {
  # Purpose:
  #
  # Arguments:
  # opt -
  # pd  -
  # ver -
  # Returns:
  #

  R = nrow( pd )
  C = ncol( pd$y )

  if ( is.null(opt) ) {

    if ( ver == 'pch' ) {

      if ( length( unique( pd$cv ) ) <= 9 ) {
        inc = 1
        lev = rep( " ", nrow( pd ) )
        for ( u in unique( pd$cv ) ) {
          lev[ u == pd$cv ] = as.character( inc )
          inc = inc + 1
        }
        out = matrix( lev, R, C )
      } else {
        out = matrix( 19, R, C )
      }
    }

    if ( ver == 'cex' ) {
      out = matrix( 1, R, C )
    }

    if ( ver == 'col' ) {
      out = matrix( 'black', R, C )
    }

    if ( ver == 'bg' ) {
      out = matrix( 'black', R, C )
    }

    if ( ver == 'lty' ) {
      out = matrix( 1, R, C )
    }

    if ( ver == 'lwd' ) {
      out = matrix( 1, R, C )
    }

  } else out = opt

  return( out )
}
