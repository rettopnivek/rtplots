#----------------------------------#
# Functions for plotting data from #
# simple choice and response time  #
# tasks                            #
#----------------------------------#

# Useful functions for package creation
# library(devtools)
# library(roxygen2)

# Index
# Lookup - 01:  rtplots
# Lookup - 02:  is.rtplots
# Lookup - 03:  lines.rtplots
# Lookup - 04:  points.rtplots
# Lookup - 05:  plot.rtplots

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
#'   the respons time variable, and optionally the choice/accuracy
#'   variable and a grouping factor variable (e.g., subjects).
#'   If \code{NULL}, the algorithm will attempt to locate response
#'   time and choice/accuracy variables via a set of default names.
#' @param keep a logical vector, indicating which rows should be
#'   kept when extracting the variables.
#' @param type the type of function to compute. Options are
#'   \itemize{
#'     \item \code{CDF}; the cumulative distribution function.
#'     \item \code{QPE}; the quantile function.
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

  output = extract_var( df, label, keep )
  output = sort_input( output, level )
  npfd = extract_npfd( output )

  ### Empirical cumulative distribution function ###

  if ( check_for_cdf_type(type) ) {
    tmp = create_cdf_output( output, npfd, ... )
    output$pd = tmp$pd
    output$xm = tmp$xm
    output$ym = tmp$ym
    output$npfd = npfd
    output$type = 'CDF'
  } else if ( check_for_qpe_type(type) ) {
    tmp = create_qpe_output( output, npfd, ... )
    output$pd = tmp$pd
    output$xm = tmp$xm
    output$ym = tmp$ym
    output$npfd = npfd
    output$type = 'QPE'
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

  if ( !is.null( val ) ) {
    vl = val == object$pd$v
  } else {
    vl = object$val[ ind ] == object$pd$v
  }

  if ( !flip ) {
    points( object$pd$x[vl], object$pd$y[vl], ... )
  } else {
    points( object$pd$x[vl], -object$pd$y[vl], ... )
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
#' @param type if eqaul to \code{'blank'}, a completely
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
