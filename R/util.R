#-------------------#
# Utility functions #
#-------------------#

# Internal functions that should not be exported

# Index
# Lookup - 01: lnInterp

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
