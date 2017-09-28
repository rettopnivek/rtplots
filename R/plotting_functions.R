#----------------------------------#
# Functions for plotting data from #
# simple choice and response time  #
# tasks                            #
#----------------------------------#

# Useful functions for package creation
# library(devtools)
# library(roxygen2)

# * Incomplete documentation

# Index
# Lookup - 01:  rtplots*
# Lookup - 02:  is.rtplots
# Lookup - 03:  lines.rtplots
# Lookup - 04:  points.rtplots
# Lookup - 05:  plot.rtplots
# Lookup - 06:  uncertainty*
# Lookup - 07:  draw_ui
# Lookup - 08:  add_to_rtplot*
# Lookup - 09

# Lookup - 01
#' Creates a 'rtplots' Object
#'
#' Creates a 'rtplots' object from a data frame of
#' response time and choice/accuracy data, allowing
#' for subsequent plotting.
#'
#' @param df a data frame that includes (at a minimum) a column of
#'   response times.
#' @param label a character string giving the column name for
#'   the response time variable, and optionally the choice/accuracy
#'   variable, a grouping factor variable (e.g., subjects), and
#'   finally an additional covariate (e.g., conditions). If
#'   \code{NULL}, the algorithm will attempt to locate response
#'   time and choice/accuracy variables via a set of default names.
#' @param keep a logical vector, indicating which rows should be
#'   kept when extracting the variables.
#' @param type the type of function to compute. Options are
#'   \itemize{
#'     \item \code{CDF}: the cumulative distribution function.
#'     \item \code{QPE}: the quantile function.
#'     \item \code{PDF}: the probability density function.
#'     \item \code{PVT}: probability versus time.
#'   }
#' @param level an optional index indicating the specific level
#'   of the grouping factor to consider. This allows, for instance,
#'   selecting a single subject.
#' @param ... additional parameters for the density, distribution,
#'   and quantile functions.
#'
#' @details
#'
#' Forthcoming
#'
#' @return An object of class 'rtplots', a list consisting of...
#'   \itemize{
#'     \item \code{pd}; a data frame with the x and y values
#'       to plot and their associated choice/accuracy value.
#'     \item \code{xm}; If relevant, a matrix whose columns
#'       contain the individual x-axis values of the function
#'       which were collapsed over the grouping factor.
#'     \item \code{ym}; If relevant, a matrix whose columns
#'       contain the individual y-axis values of the function
#'       which were collapsed over the grouping factor.
#'     \item \code{npfd}; a data frame with the number of
#'       observations, proportions, and frequencies based on
#'       the choice/accuracy values.
#'     \item \code{type}; the type of function that was applied.
#'     \item \code{x}; the range of x-axis values.
#'     \item \code{y}; the range of y-axis values.
#'   }
#'
#' @references
#'
#' Forthcoming
#'
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' lbl = c( 'RT', 'Accuracy', 'Subjects' )
#' keep = priming_data$Condition == 0
#' obj = rtplots( priming_data, label = lbl, keep = keep )
#'
#' @export
rtplots = function( df, label = NULL, keep = NULL,
                    type = 'CDF', level = NULL, ... ) {

  # Check input
  if ( !is.data.frame( df ) )
    stop( 'Input variable must be a data frame' )

  output = extract_var( df, label, keep, type )
  # Sort RTs into ascending order
  output = sort_input( output, level )
  # Compute proportions and frequencies
  npfd = extract_npfd( output )

  if ( check_for_cdf_type(type) ) {
    ### CDF

    # Generate output
    tmp = create_cdf_output( output, npfd, ... )
    output$pd = tmp$pd
    output$xm = tmp$xm
    output$ym = tmp$ym
    output$npfd = npfd
    output$type = 'CDF'

  } else if ( check_for_qpe_type(type) ) {
    ### QPE

    # Generate output
    tmp = create_qpe_output( output, npfd, ... )
    output$pd = tmp$pd
    output$xm = tmp$xm
    output$ym = tmp$ym
    output$npfd = npfd
    output$type = 'QPE'

  } else if ( check_for_pdf_type(type) ) {
    ### PDF

    # Generate output
    tmp = create_pdf_output( output, npfd, ... )
    output$pd = tmp$pd
    output$xm = tmp$xm
    output$ym = tmp$ym
    output$npfd = npfd
    output$type = 'PDF'
  } else if ( check_for_pvt_type(type) ) {
    ### PVT

    # Generate output
    tmp = create_pvt_output( output, npfd, ... )
    output$pd = tmp$pd
    output$xm = tmp$xm
    output$ym = tmp$ym
    output$npfd = npfd
    output$type = 'PVT'
  } else {
    err = 'That function type is unknown.'
    stop( err, call. = FALSE )
  }

  # Plotting range
  output$x = range( na.omit( output$pd$x ) )
  output$y = range( output$pd$y )

  class( output ) = 'rtplots'

  return( output )
}

# Lookup - 02
#' @rdname rtplots
#' @export
is.rtplots = function(x) inherits(x, "rtplots")

# Lookup - 03
#' Lines Method for 'rtplots' Object
#'
#' Provides a method to 'rtplots' objects that allows
#' for drawing line segments on an existing plot.
#'
#' @param object a 'rtplots' object.
#' @param val the choice or accuracy value whose observations
#'   are to be plotted.
#' @param ind a index to select observations for a specific
#'   choice/accuracy value.
#' @param flip logical; if true, flips the y-axis values to
#'   be negative.
#' @param ... additional plotting parameters. See
#'   \code{\link[graphics]{lines}}.
#'
#' @export
lines.rtplots = function( object, val = NULL,
                          ind = 1, flip = F, ... ) {

  # If it isn't a probability versus time plot
  if ( object$type != 'PVT' ) {

    if ( !is.null( val ) ) {
      vl = val == object$pd$v
    } else {
      vl = object$val[ ind ] == object$pd$v
    }

    if ( !flip ) {
      lines( object$pd$x[vl], object$pd$y[vl], ... )
    } else {
      lines( object$pd$x[vl], -object$pd$y[vl], ... )
    }
  } else { # Special options for PvT figures

    # Extract choice value
    if ( is.null( val ) ) val = object$val[ ind ]

    # Check for multiple RTs per probability
    dmn = dim( object$pd$y )

    # If the points aren't flipped about the x-axis
    if ( !flip ) {

      # Subset of values matching desired choice value
      vl = which( object$pd$v %in% val )

      for ( i in 1:dmn[2] ) {
        x = object$pd$x[vl]
        y = object$pd$y[vl,i]
        ord = order(x)
        lines( x[ord], y[ord], ... )
      }

    } else { # Flip

      # Subset of values matching desired choice value
      vl = which( object$pd$v %in% val )

      for ( i in 1:dmn[2] ) {
        x = object$pd$x[vl]
        y = object$pd$y[vl,i]
        ord = order(x)
        lines( x[ord], -y[ord], ... )
      }

    }
  }
}

# Lookup - 04
#' Points Method for 'rtplots' Object
#'
#' Provides a method to 'rtplots' objects that allows
#' for drawing points on an existing plot.
#'
#' @param object a 'rtplots' object.
#' @param val the choice or accuracy value whose observations
#'   are to be plotted.
#' @param ind a index to select observations for a specific
#'   choice/accuracy value.
#' @param flip logical; if true, flips the y-axis values to
#'   be negative.
#' @param ... additional plotting parameters. See
#'   \code{\link[graphics]{points}}.
#'
#' @export
points.rtplots = function( object, val = NULL,
                           ind = 1, flip = F, ... ) {

  # If it isn't a probability versus time plot
  if ( object$type != 'PVT' ) {

    # If a choice value is given
    if ( !is.null( val ) ) {
      vl = val == object$pd$v
    } else {
      vl = object$val[ ind ] == object$pd$v
    }

    # If the points aren't flipped about the x-axis
    if ( !flip ) {
      points( object$pd$x[vl], object$pd$y[vl], ... )
    } else { # Flip
      points( object$pd$x[vl], -object$pd$y[vl], ... )
    }

  } else { # Special options for PvT figures

    # Extract choice value
    if ( is.null( val ) ) val = object$val[ ind ]

    # Check for multiple RTs per probability
    dmn = dim( object$pd$y )

    # If the points aren't flipped about the x-axis
    if ( !flip ) {

      # Subset of values matching desired choice value
      vl = which( object$pd$v %in% val )

      for ( i in 1:length(vl) ) {
        x = rep( object$pd$x[vl[i]], dmn[2] )
        y = object$pd$y[vl[i],]
        points( x, y, ... )
      }

    } else { # Flip

      # Subset of values matching desired choice value
      vl = which( object$pd$v %in% val )

      for ( i in 1:length(vl) ) {
        x = rep( object$pd$x[vl[i]], dmn[2] )
        y = object$pd$y[vl[i],]
        points( x, -y, ... )
      }

    }
  }
}

# Lookup - 05
#' Plot Method for 'rtplots' Object
#'
#' Creates a blank plot based on a 'rtplots' object.
#'
#' @param x a 'rtplots' object.
#' @param y the y-axis boundaries.
#' @param inc the interval for the ticks on the
#'   x and y-axes (depending on the plot).
#' @param type if equal to \code{'blank'}, a completely
#'   blank plot is generated. Otherwise, converted to \code{'n'}.
#' @param xlab a character string for the label on the x-axis.
#' @param ylab a character string for the label on the y-axis.
#' @param bty the type of box to draw around the plot. See
#'   \code{\link[graphics]{par}}.
#' @param ... additional plotting parameters. See
#'   \code{\link[graphics]{par}}.
#'
#' @export
plot.rtplots = function( x,
                         y = NULL,
                         inc = NULL,
                         type = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         bty = NULL, ... ) {

  # Extract 'rtplots' object
  object = x

  if ( object$type == 'CDF' ) {

    # Determine x and y-axis boundaries
    if ( is.null( inc ) ) inc = .2
    xl = lower_upper( inc[1], object$pd$x )
    if ( is.null( y ) ) yl = c(0,1) else yl = y

    # Default options for plot
    if ( is.null( xlab ) ) xlb = 'Time' else xlb = xlab
    if ( is.null( ylab ) ) ylb = 'Distribution function'
    else ylb = ylab
    if ( is.null( bty ) ) bx = 'l' else bx = bty

    plotYes = T

    # Plot type
    if ( !is.null( type ) ) {

      if ( type == 'blank' | type == 'Blank' ) {
        plot( xl, yl, type = 'n', xlab = ' ',
              ylab = ' ', xaxt = 'n', yaxt = 'n',
              bty = 'n' )
      }

      plotYes = F

    } else type = 'n'

    if ( plotYes )
      plot( xl, yl, type = type, xlab = xlb, ylab = ylb,
            bty = bx, ... )

  }

  if ( object$type == 'PDF' ) {

    # Determine x and y-axis boundaries
    if ( is.null( inc ) ) inc = .2
    xl = lower_upper( inc[1], object$pd$x )
    if ( is.null( y ) ) yl = c(0,2) else yl = y

    # Default options for plot
    if ( is.null( xlab ) ) xlb = 'Time' else xlb = xlab
    if ( is.null( ylab ) ) ylb = 'Density function'
    else ylb = ylab
    if ( is.null( bty ) ) bx = 'l' else bx = bty

    plotYes = T

    # Plot type
    if ( !is.null( type ) ) {

      if ( type == 'blank' | type == 'Blank' ) {
        plot( xl, yl, type = 'n', xlab = ' ',
              ylab = ' ', xaxt = 'n', yaxt = 'n',
              bty = 'n' )
      }

      plotYes = F

    } else type = 'n'

    if ( plotYes )
      plot( xl, yl, type = type, xlab = xlb, ylab = ylb,
            bty = bx, ... )

  }

  if ( object$type == 'QPE' ) {

    # Determine x and y-axis boundaries
    if ( is.null( inc ) ) inc = .2
    xl = lower_upper( inc, na.omit( object$pd$x ) )
    if ( is.null( y ) ) yl = c(0,1) else yl = y

    # Default options for plot
    if ( is.null( xlab ) ) xlb = 'Time' else xlb = xlab
    if ( is.null( ylab ) ) ylb = 'Quantile function'
    else ylb = ylab
    if ( is.null( bty ) ) bx = 'l' else bx = bty

    plotYes = T

    # Plot type
    if ( !is.null( type ) ) {

      if ( type == 'blank' | type == 'Blank' ) {
        plot( xl, yl, type = 'n', xlab = ' ',
              ylab = ' ', xaxt = 'n', yaxt = 'n',
              bty = 'n' )
      }

      plotYes = F

    } else type = 'n'

    if ( plotYes )
      plot( xl, yl, type = type, xlab = xlb, ylab = ylb,
            bty = bx, ... )

  }

  if ( object$type == 'PVT' ) {

    # Determine x and y-axis boundaries
    if ( is.null( inc ) ) inc = .2
    yl = lower_upper( inc, as.vector( na.omit( object$pd$y ) ) )
    if ( is.null( y ) ) xl = c(0,1) else xl = y

    # Default options for plot
    if ( is.null( xlab ) ) xlb = 'Probability' else xlb = xlab
    if ( is.null( ylab ) ) ylb = 'Time'
    else ylb = ylab
    if ( is.null( bty ) ) bx = 'l' else bx = bty

    plotYes = T

    # Plot type
    if ( !is.null( type ) ) {

      if ( type == 'blank' | type == 'Blank' ) {
        plot( xl, yl, type = 'n', xlab = ' ',
              ylab = ' ', xaxt = 'n', yaxt = 'n',
              bty = 'n' )
      }

      plotYes = F

    } else type = 'n'

    if ( plotYes )
      plot( xl, yl, type = type, xlab = xlb, ylab = ylb,
            bty = bx, ... )

  }

}

# Lookup - 06
#' Computes Uncertainty Intervals for 'rtplots' Objects
#'
#' A function that computes a specified uncertainty
#' interval (such as \eqn{\pm} 2 standard errors from
#' the mean) for a 'rtplots' object.
#'
#' @param object a 'rtplots' object.
#' @param f a function to compute the boundaries of
#'   the uncertainty interval given the desired
#'   range of coverage (takes at least two parameters,
#'   a vector \code{x} and a range \code{interval}).
#' @param alpha the desired width of coverate for the
#'   given uncertainty interval.
#' @param ... additional parameters for the
#'   function \code{f}.
#'
#' @return Forthcoming
#'
#' @export
uncertainty = function( object, f = NULL,
                        alpha = .95, ... ) {

  # If no function is provided, assume uncertainty
  # intervals should be calculated using the
  # standard error of the mean
  if ( is.null( f ) ) {
    f = function(x,interval) {
      s = sd( x )/sqrt( length(x) ) # Standard error
      m = mean(x)
      return( c( m - abs( qnorm(interval[1]) )*s,
                 m + abs( qnorm(interval[2] ) )*s ) )
    }
  }

  # Determine lower and upper boundaries
  interval = numeric(2);
  interval[1] = (1 - alpha)/2
  interval[2] = interval[1] + alpha

  # Initialize output
  ud = NULL
  pd = NULL

  # If there is more than one group
  if ( object$n_g > 1 ) {

    l = nrow( object$pd )

    ud = data.frame( x = rep( NA, l * 3 ),
                     y = rep( NA, l * 3 ),
                     loc = rep( NA, l * 3 ),
                     v = rep( NA, l * 3 ),
                     i = rep( NA, l * 3 ) )

    # Loop over choice/accuracy values
    sta = 1; end = 0
    for ( i in 1:(object$n_v) ) {

      # Compute uncertainty intervals
      if ( !is.null( object$xm ) ) {
        xui = apply( object$xm[[i]], 2, f, interval = interval,
                     ... )
      } else xui = NULL
      if ( !is.null( object$ym ) ) {
        yui = apply( object$ym[[i]], 2, f, interval = interval,
                     ... )
      } else yui = NULL

      sel = object$val[ i ] == object$pd$v

      # Center values
      end = end + sum( sel )
      ud$x[ sta:end ] = object$pd$x[ sel ]
      ud$y[ sta:end ] = object$pd$y[ sel ]
      ud$v[ sta:end ] = object$pd$v[ sel ]
      ud$loc[ sta:end ] = 'center'
      ud$i[ sta:end ] = 1:sum(sel)
      sta = end + 1

      # Lower boundary
      end = end + sum( sel )
      if ( is.null( xui ) )
        ud$x[ sta:end ] = object$pd$x[ sel ] else
          ud$x[ sta:end ] = xui[1,]
      if ( is.null( yui ) )
        ud$y[ sta:end ] = object$pd$y[ sel ] else
          ud$y[ sta:end ] = yui[1,]
      ud$v[ sta:end ] = object$pd$v[ sel ]
      ud$loc[ sta:end ] = 'lower'
      sta = end + 1

      # Upper boundary
      end = end + sum( sel )
      if ( is.null( xui ) )
        ud$x[ sta:end ] = object$pd$x[ sel ] else
          ud$x[ sta:end ] = xui[2,]
      if ( is.null( yui ) )
        ud$y[ sta:end ] = object$pd$y[ sel ] else
          ud$y[ sta:end ] = yui[2,]
      ud$v[ sta:end ] = object$pd$v[ sel ]
      ud$loc[ sta:end ] = 'upper'
      sta = end + 1

    }

    pd = data.frame( xl = rep( NA, l * 2 ),
                     xu = rep( NA, l * 2 ),
                     yl = rep( NA, l * 2 ),
                     yu = rep( NA, l * 2 ) )
    pd$xl = c( ud$x[ ud$loc == 'lower' ],
               ud$x[ ud$loc == 'center' ] )
    pd$xu = c( ud$x[ ud$loc == 'upper' ],
               ud$x[ ud$loc == 'center' ] )
    pd$yl = c( ud$y[ ud$loc == 'center' ],
               ud$y[ ud$loc == 'lower' ] )
    pd$yu = c( ud$y[ ud$loc == 'center' ],
               ud$y[ ud$loc == 'upper' ] )
    pd$v = rep( ud$v[ ud$loc == 'center' ], 2 )
    pd$i = rep( ud$i[ ud$loc == 'center' ], 2 )

  }

  out = list(
    pd = pd,
    ud = ud,
    xm = object$xm,
    ym = object$ym,
    val = object$val
  )

  # Define new class
  class( out ) = 'rtplots_ui';

  return( out )
}

# Lookup - 07
#' Draw Uncertainty Intervals
#'
#' A convenience function that draws uncertainty
#' intervals computed from a 'rtplots' object onto
#' an already existing plot.
#'
#' @param ui a 'rtplots_ui' object.
#' @param val the choice/accuracy value for which to
#'   draw the uncertainty intervals.
#' @param ind the index of the desired choice/accuracy
#'   value over which to draw uncertainty intervals.
#' @param flip logical; if \code{true} the y-axis
#'   values are flipped to be negative.
#' @param type the type of method used to draw the
#'   uncertainty interval. Options are
#'   \itemize{
#'     \item \code{segments}.
#'     \item \code{arrows}.
#'     \item \code{polygon}.
#'   }
#' @param ... additional plotting parameters. See
#'   \code{\link[graphics]{segments}},
#'   \code{\link[graphics]{arrows}}, and
#'   \code{\link[graphics]{polygon}}.
#'
#' @export
draw_ui = function( ui, val = NULL, ind = 1,
                    flip = F, type = 'segments', ... ) {

  if ( !is.null( val ) ) {
    vl = val == ui$pd$v
  } else {
    vl = ui$val[ ind ] == ui$pd$v
  }
  if (!flip) w = 1 else w = -1

  segment_names = c( 'segments', 's', 'seg' )
  if ( type %in% segment_names ) {
    segments( ui$pd$xl[vl], w*ui$pd$yl[vl],
              ui$pd$xu[vl], w*ui$pd$yu[vl], ... )
  }

  polygon_names = c( 'polygon', 'p', 'poly' )
  if ( type %in% polygon_names ) {

    polygon( c( ui$pd$xl[vl], rev( ui$pd$xu[vl] ) ),
             c( ui$pd$yu[vl], rev( ui$pd$yl[vl] ) ),
             ... )
  }

  arrow_names = c( 'arrows', 'a', 'arr' )
  if ( type %in% arrow_names ) {
    arrows( ui$pd$xl[vl], w*ui$pd$yl[vl],
            ui$pd$xu[vl], w*ui$pd$yu[vl], ... )
  }

}

# Lookup - 08
#' Add Points or Lines to a RT Plot
#'
#' A convenience function that draws points or
#' line segments onto an already existing plot. Useful
#' for adding summary statistics (e.g., the mean or
#' median).
#'
#' @param object a rtplots object or a list with the
#'   x and y-axis values to plot.
#' @param T_x the test statistic(s) to be added to the
#'   existing figure.
#' @param T_g the test statistic to use when collapsing
#'   over a grouping factor.
#' @param val Forthcoming.
#' @param ind Forthcoming.
#' @param type Forthcoming.
#' @param out Forthcoming.
#' @param ... Forthcoming.
#'
#' @export
add_to_rtplot = function( object, T_x = mean, T_g = mean, val = NULL,
                          ind = 1, type = 'p', out = T, ... ) {

  # Function to calculate values for y-axis
  # using linear interpolation
  f = function( x ) {

    ind = max( which( xa < x ) )
    if ( ind == -Inf ) {
      out = ya[1];
    } else {

      if ( ind == length(xa) ) {
        out = ya[ length(xa) ]
      } else {
        pts = c( ind, ind + 1 )
        out = ln_interp( x, ya[ pts ], xa[ pts ] )
      }

    }

    return( out )
  }

  if ( is.rtplots( object ) ) {

    # If it isn't a probability versus time plot
    if ( object$type != 'PVT' ) {

      # If there is no grouping factor

      if ( object$n_g <= 1 ) {

        if ( !is.null( val ) ) {
          vl = val == object$ad$v
          vly = val == object$pd$v
        } else {
          vl = object$val[ ind ] == object$ad$v
          vly = object$val[ ind ] == object$pd$v
        }

        ts = T_x( object$ad$t[vl] )

        # Extract plotting values
        xa = object$pd$x[vly] # x-axis
        ya = object$pd$y[vly] # y-axis

        # Determine y-axis values
        ny = sapply( ts, f )

        output = list( x = ts, y = ny )

      } else {

        if ( !is.null( val ) ) {
          vl = val == object$ad$v
          vly = val == object$pd$v
        } else {
          vl = object$val[ ind ] == object$ad$v
          vly = object$val[ ind ] == object$pd$v
        }

        dmn = length( T_x( object$ad$t[vl] ) )
        ts_all = by( object$ad$t[vl], list( object$ad$g[vl] ), T_x )
        ts_all = matrix( unlist( ts_all ), length( ts_all ), dmn, byrow = T )
        ts = apply( ts_all, 2, T_g )

        # Extract plotting values
        xa = object$pd$x[vly] # x-axis
        ya = object$pd$y[vly] # y-axis

        # Determine y-axis values
        ny = sapply( ts, f )

        output = list( x = ts, y = ny )

      }
    }
  } else {
    pass = F
    if ( is.list( object ) ) {
      if ( all( c('x','y') %in% names( object ) ) ) {
        pass = T
        output = object
      }
      if (!pass) stop( 'Please pass in an "rtplots" object',
                       .call = FALSE )
    }
  }

  if ( type == 's' | type == 'segments' ) {
    segments( output$x, 0, output$x, output$y, ... )
  }
  if ( type == 'p' | type == 'points' ) {
    points( output$x, output$y, ... )
  }

  if ( out ) return( output )
}

# Lookup - 09
#' Create a Probability versus Time Figure
#'
#' A convenience function that draws points, lines, or
#' line segments for a probability versus time figure
#' onto an already existing plot. Allows fine-tune
#' control of plotting parameters such as point type,
#' size, etc.
#'
#' @param object a rtplots object.
#'
#'
#' @export
pvt_plot = function( object, keep = NULL, ind = NULL, val = NULL, type = 'p',
                     pch = NULL, cex = NULL, col = NULL,
                     bg = NULL, lwd = NULL ) {

  if ( is.rtplots( object ) ) {

    if ( object$type == 'PVT' ) {

      # If a choice value is given
      if ( !is.null( val ) ) {
        vl = val == object$pd$v
      } else if ( !is.null( ind ) ) {
        vl = object$val[ ind ] == object$pd$v
      } else {
        vl = rep( T, nrow( object$pd ) )
      }

      # Extract relevant data
      if ( is.null( keep ) ) keep = rep( T, nrow( object$pd ) )

      # Extract relevant data
      pd = object$pd[ vl & keep, ]

      # Dimensions
      R = nrow( pd )
      C = dim( pd$y )[2]

      if ( type == 'p' | type == 'points' ) {

        # Extract plotting options
        pch = pvt_default_options( pch, pd, ver = 'pch' )
        cex = pvt_default_options( cex, pd, ver = 'cex' )
        col = pvt_default_options( col, pd, ver = 'col' )
        bg = pvt_default_options( bg, pd, ver = 'bg' )

        for ( r in 1:R ) {
          x = rep( pd$x[r], C )
          y = pd$y[r,]
          points( x, y, pch = pch[r,], col = col[r,], cex = cex[r,], bg = bg[r,] )
        }

      }

      if ( type == 'l' | type == 'lines' ) {

        # Extract plotting options
        lty = pvt_default_options( pch, pd, ver = 'lty' )
        lwd = pvt_default_options( cex, pd, ver = 'lwd' )
        col = pvt_default_options( col, pd, ver = 'col' )

        # Loop over test statistic values
        for ( c in 1:C ) {
          x = pd$x
          y = pd$y[,c]
          ord = order( x )
          lines( x[ord], y[ord], lty = lty[ ord, c ],
                 lwd = lwd[ ord, c ], col = col[ ord, c ] )
        }
      }

      if ( type == 's' | type == 'segments' ) {

        # Extract plotting options
        lty = pvt_default_options( pch, pd, ver = 'lty' )
        lwd = pvt_default_options( cex, pd, ver = 'lwd' )
        col = pvt_default_options( col, pd, ver = 'col' )

        for ( u in unique( pd$cv ) ) {

          sel = pd$cv == u
          for ( c in 1:C ) {
            lines( pd$x[sel], pd$y[sel,c],
                   lty = lty[ sel, c ], lwd = lwd[ sel, c ],
                   col = col[ sel, c ] )
          }

        }

      }

    } else stop( 'No values for a PvT plot', .call = FALSE )

  } else {
    stop( 'Not an rtplots object', .call = FALSE )
  }

}
