#---------------#
# QPE functions #
#---------------#

# Index
# Lookup - 01:  ?
# Lookup - 02:  ?
# Lookup - 03:  ?
# Lookup - 04:  ?

# Lookup - 01
j_qpe_one_level = function( input, npfd, prb ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...

  l = length( prb ) * input$n_v
  pd = data.frame( x = rep( NA, l ),
                   y = rep( NA, l ),
                   v = rep( input$val, each = length( prb ) ) )

  for ( vl in input$val ) {
    sel1 = input$ad$v == vl
    sel2 = pd$v == vl
    if ( length( sel1 ) > 0 ) {
      pd$x[sel2] = quantile( input$ad$t[ sel1 ], prb );
      sel_v = grep( as.character( vl ), colnames( npfd$P ) )
      pd$y[sel2] = prb * npfd$P[ sel_v ]
    } else {
      pd$y[sel2] = 0;
    }
  }
  out = list( pd = pd, xm = NULL, ym = NULL )

  return( out )
}

# Lookup - 02
j_qpe_group = function( input, npfd, prb,
                        T_x = mean ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...

  # Initialize output
  l = length( prb ) * input$n_v
  pd = data.frame( x = rep( NA, l ),
                   y = rep( NA, l ),
                   v = rep( input$val, each = length( prb ) ) )
  ym = c() # Create empty lists
  xm = c()

  inc = 1
  for ( vl in input$val ) {
    sel1 = input$ad$v == vl
    sel2 = pd$v == vl
    if ( length( sel1 ) > 0 ) {

      xm = c( xm, list( NULL ) )
      ym = c( ym, list( NULL ) )

      aq = by( input$ad$t[ sel1 ], list( input$ad$g[ sel1 ] ),
               quantile, prob = prb )
      mat = matrix( unlist( aq ),
                    length( unique( input$ad$g[ sel1 ] ) ),
                    length( prb ), byrow = T )
      xm[[ inc ]] = mat

      grp_q = apply( mat, 2, T_x )

      mat = matrix( prb, length( unique( input$ad$g[ sel1 ] ) ),
                    length( prb ), byrow = T )
      g_sel = names( aq )
      g_sel = as.character( npfd$G ) %in% g_sel
      v_sel = grep( vl, colnames( npfd$P ) )
      mat = mat * npfd$P[ g_sel, v_sel ]
      ym[[ inc ]] = mat

      grp_p = apply( mat, 2, T_x )
      inc = inc + 1

      pd$x[sel2] = grp_q;
      pd$y[sel2] = grp_p;

    } else {
      pd$y[sel2] = 0;
    }
  }
  names( xm ) = as.character( input$val )
  names( ym ) = as.character( input$val )

  out = list( pd = pd, xm = xm, ym = ym )

  return( out )
}

# Lookup - 03
create_qpe_output = function( input, npfd, ... ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...

  arguments = list(...)
  if ( !( 'prob' %in% names( arguments ) ) )
    prb = c( 0, seq(.1,.9,.2), 1 ) else
      prb = arguments$prob

  plt = NULL
  if ( input$n_g == 1 ) plt = j_qpe_one_level( input, npfd, prb )
  if ( input$n_g > 1 ) plt = j_qpe_group( input, npfd, prb )

  return( plt )
}

# Lookup - 05
check_for_qpe_type = function( type ) {
  # Purpose:
  # ...
  # Arguments:
  # ...
  # Returns:
  # ...

  qpe_types = c( 'QPE', 'qpe', 'Qpe', 'qPe', 'qpE',
                 'QPe', 'QpE', 'qPE', 'quantile' )

  return( type %in% qpe_types )
}
