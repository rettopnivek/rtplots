#-------------------#
# Utility functions #
#-------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01: curve_defaults

# Lookup - 01
curve_defaults = function( plt ) {
  # Purpose:
  # A helper function that creates a set of default
  # plotting parameters for the cdf_curve and pdf_curve
  # furnctions.
  # Arguments:
  # plt - A named list of up to 3 internal lists, plt$ln,
  #       plt$pt1, and plt$pt2. plt$ln is a named list giving
  #       the 'lty', 'lwd', 'col', and 'type' graphical
  #       parameters for a line. plt$pt1 and plt$pt2 are named
  #       lists giving 'pch', 'cex', 'col', 'bg', and 'lwd'
  #       graphical parameters for a set of points
  # Returns:
  # A list of three elements, the default graphical
  # parameters for a line and two sets of points.

  # Initialize defaults
  Defaults_1 = list( lty = 1, lwd = 1, col = 'black', type = 'l' )
  Defaults_2 = list( pch = 21, cex = 1, col = 'black',
                     bg = 'white', lwd = 1, type = 'p' )
  Defaults_3 = list( pch = 24, cex = 1, col = 'black',
                     bg = 'white', lwd = 1, type = 'p' )

  if ( length( plt$ln ) != 0 ) {

    lnOptions = plt$ln

    newOpt = names( lnOptions )
    oldOpt = names( Defaults_1 )

    inc = 1
    for ( i in 1:length(lnOptions) ) {

      chng = which( oldOpt %in% newOpt[inc] )
      Defaults_1[[chng]] = lnOptions[[inc]]
      inc = inc + 1
    }

  }

  if ( length( plt$pt1 ) != 0 ) {

    pt1Options = plt$pt1

    newOpt = names( pt1Options )
    oldOpt = names( Defaults_2 )

    inc = 1
    for ( i in 1:length(pt1Options) ) {

      chng = which( oldOpt %in% newOpt[inc] )
      Defaults_2[[chng]] = pt1Options[[inc]]
      inc = inc + 1
    }

  }

  inc = 1
  if ( length( plt$pt2 ) != 0 ) {

    pt2Options = plt$pt2

    newOpt = names( pt2Options )
    oldOpt = names( Defaults_3 )

    for ( i in 1:length(pt2Options) ) {

      chng = which( oldOpt %in% newOpt[inc] )
      Defaults_3[[chng]] = pt2Options[[inc]]
      inc = inc + 1
    }

  }

  return( list( Defaults_1, Defaults_2, Defaults_3 ) )
}

# Lookup - 02
point_defaults = function( plt ) {
  # Purpose:
  #
  # Arguments:
  # plt
  # Returns:
  #

  # Initialize defaults
  Defaults_1 = list( lty = 1, lwd = 1, col = 'black', type = 'l' )
  Defaults_2 = list( pch = 21, cex = 1, col = 'black',
                     bg = 'white', lwd = 1,  type = 'p' )

  if ( length( plt$ln ) != 0 ) {

    lnOptions = plt$ln

    newOpt = names( lnOptions )
    oldOpt = names( Defaults_1 )

    inc = 1
    for ( i in 1:length(lnOptions) ) {

      chng = which( oldOpt %in% newOpt[inc] )
      Defaults_1[[chng]] = lnOptions[[inc]]
      inc = inc + 1
    }

  }

  if ( length( plt$pt ) != 0 ) {

    pt1Options = plt$pt

    newOpt = names( pt1Options )
    oldOpt = names( Defaults_2 )

    inc = 1
    for ( i in 1:length(pt1Options) ) {

      chng = which( oldOpt %in% newOpt[inc] )
      Defaults_2[[chng]] = pt1Options[[inc]]
      inc = inc + 1
    }

  }

  return( list( Defaults_1, Defaults_2 ) )
}

lnInterp = function(x,yPts,xPts) {
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
  # The predicted y-axis point

  # Determine
  b = diff(yPts)/diff(xPts)
  y = yPts[1] + b*(x-xPts[1])

  return( y )
}

