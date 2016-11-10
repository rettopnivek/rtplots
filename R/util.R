#-------------------#
# Utility functions #
#-------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01: lnInterp
# Lookup - 02: addEllipse
# Lookup - 03: pvt_plot_options

# Lookup - 01
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

# Lookup - 02
addEllipse = function(a, b, Xc, Yc, k = 0, deg = T, ... ) {
  # Purpose:
  # Draws an ellipse on an already existing plot.
  # Arguments:
  # a  - the length of the x-axis vertice of the ellipse
  # b  - the length of the y-axis vertice of the ellipse
  # k  - the angle between the x-axis in the major vertice
  #      for the ellipse
  # Xc - the center on the x-axis for the ellipse
  # Yc - the center on the y-axis for the ellipse

  if (deg) k = k * (pi/180) # Convert to radians

  t = seq(0,2*pi,length=100)

  x = Xc + a * cos(t) * cos(k) - b * sin(t) * sin(k)
  y = Yc + b * sin(t) * cos(k) + a * cos(t) * sin(k)
  polygon(x,y,...)

}

pvt_plot_options = function( plt, dmn, xv ) {
  # Purpose:
  #
  # Arguments:
  # plt -
  # dmn -
  # xv  -
  # Returns:
  # An updated version of the list plt.

  # If the y-axis values are stored in a matrix
  if ( length( dmn ) > 0 ) {

    # Extract the pch options
    if ( length( plt$pch ) == 0 ) {
      plt$pch = matrix( as.character(1:dmn[1]), dmn[1], dmn[2] )
    } else {
      if ( length( dim( plt$pch ) ) == 0 ) {
        plt$pch = matrix( plt$pch, dmn[1], dmn[2] )
      }
    }

    # Extract the col options
    if ( length( plt$col ) == 0 ) {
      plt$col = matrix( 'black', dmn[1], dmn[2] )
    } else {
      if ( length( dim( plt$col ) ) == 0 ) {
        plt$col = matrix( plt$col, dmn[1], dmn[2] )
      }
    }

    # Extract the bg options
    if ( length( plt$bg ) == 0 ) {
      plt$bg = matrix( 'white', dmn[1], dmn[2] )
    } else {
      if ( length( dim( plt$bg ) ) == 0 ) {
        plt$bg = matrix( plt$bg, dmn[1], dmn[2] )
      }
    }

  }

  # If the y-axis values are stored in a vector
  if ( length( dmn ) == 0 ) {

    # Extract the pch options
    if ( length( plt$pch ) == 0 ) {
      plt$pch = as.character(1:length(xv$P))
    } else {
      plt$pch = plt$pch
    }

    # Extract the col options
    if ( length( plt$col ) == 0 ) {
      plt$col = 'black'
    } else {
      plt$col = plt$col
    }

    # Extract the bg options
    if ( length( plt$bg ) == 0 ) {
      plt$bg = 'white'
    } else {
      plt$bg = plt$bg
    }

  }

  return( plt )
}
