#----------------------------------#
# Functions for plotting data from #
# simple choice and response time  #
# tasks                            #
#----------------------------------#

# Index
# Lookup - 01:  cdf_curve
# Lookup - 02:  pdf_curve
# Lookup - 03:  hazard_curve
# Lookup - 04:  quantile_points
# Lookup - 05:  CAF_points
# Lookup - 06:  add_points
# Lookup - 07:  blankRTplot

### TO DO ###
# Add accuracy-latency plots
# Add helper functions for plotting variability (e.g. SEs)
# Check man pages for errors
# Add references (base R density function?)
# Write script testing the functions

# Useful functions for package creation
# library(devtools)
# library(roxygen2)

# Lookup - 01
#' CDF curves for response time and choice data.
#'
#' Draws a line for the empirical CDF of a set of response
#' times (conditioned on choice) on an already existing plot.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the CDF for responses times when choice = 0 is
#'   drawn. If 1, the CDF for responses times when choice = 1 is drawn.
#' @param grp an optional vector with a grouping factor (e.g. subjects).
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{jnt}}{If true, the joint distribution is used.}
#'     \item{\code{draw}}{If true, the curve is drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve is flipped about the
#'       x-axis.}
#'   }
#' @param ... additional plotting parameters.
#' @return A list consisting of ...
#' \describe{
#'   \item{\code{pv}}{a data frame with the plotting values
#'     used for the x-axis and the y-axis.}
#'   \item{\code{g}}{when a grouping factor is present, a list
#'     with the matrix of y-axis values per level, the vector
#'     of associated response times, and the total number of
#'     observations per each included level.}
#'   \item{\code{v}}{a list of additional variables, the total
#'     number of observations for all levels, the choice proportion
#'     for each level, and the choice selection.}
#'   \item{\code{i}}{a list of the input variables.}
#'   \item{\code{opt}}{a list of the options used.}
#'   }
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 4 & d$Subject == 1
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot()
#' cdf_curve( rt, ch )
#' cdf_curve( rt, ch, sel = 0, lty = 2 )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 0 # Not all subjects had responses for each choice
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( bty = 'l', ver = 'CDF', cex.axis = 1.5, cex.lab = 1.5 )
#' cdf_curve( rt, ch, grp = grp, lwd = 2 )
#' cdf_curve( rt, ch, sel = 0, grp = grp, lwd = 2, lty = 2 )
#' @export

cdf_curve = function( rt, ch, sel = 1, grp = NULL,
                      opt = list( ), ... ) {

  # Set options for joint distribution, drawing, output, and
  # whether curve should be flipped around x-axis
  if ( length( opt$jnt ) == 0 ) jnt = T else jnt = opt$jnt
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  # Save options
  optOut = list( jnt = jnt, draw = draw,
                 out = out, flip = flip )

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    # Total number of observations
    n = sum( ch == sel )
    if ( n < 1 ) stop('No observations')

    # Empirical estimate of CDF for response times
    p = (1:n)/n
    x = sort( rt[ ch == sel ] )

    # Adjust asymptote if estimating joint CDF
    adj = n/length(ch)
    if (jnt) p = p*adj;

    # Create list for output
    output = list(
      # Save plotting values
      pv = as.data.frame( cbind( x = x, y = p ) ),
      # Save grouping factor info
      g = NULL,
      # Save variables for calculations
      v = list( n = n, adj = adj, sel = sel ),
      # Save input
      i = as.data.frame( cbind( rt = rt, ch = ch ) ),
      # Save options
      opt = optOut
    )

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    # Determine number of observations per level of
    # grouping factor
    n = aggregate( ch == sel, list( grp ), sum )$x
    # Determine if there are sufficient observations
    cmp = round( mean(n) )
    if (cmp < 1) stop('Not enough observations')

    # Determine number of points to use to estimate CDF
    if ( cmp < 20 ) seqL = cmp
    if ( cmp >= 20 & cmp < 100 ) seqL = 20;
    if ( cmp >= 100 ) seqL = 100;

    # Condition on choice
    xAll = rt[ ch == sel ]
    g = grp[ ch == sel ]

    # Determine range to estimate CDF
    xRange = aggregate( xAll, list( g ),
                        function(x) c( min(x), max(x) ) )
    tmp = colMeans( xRange$x )

    # Generate sequence of times over which to calculate
    # curve
    xVal = seq( tmp[1], tmp[2], length = seqL  )

    # Determine empirical CDF for each subject
    pAll = matrix( NA, length( unique( g ) ), seqL  )
    for (i in 1:seqL ) {

      cnt = sapply( xAll, function(x) x < xVal[i] )
      pAll[,i] = aggregate( cnt, list( g ), sum )$x/n[ n > 0 ]

    }

    # Aggregate over subjects
    tot = aggregate( rep(1,length(rt)), list( grp ), sum )$x
    adj = n[ n > 0 ]/tot[ n > 0 ]

    # Adjust asymptote if estimating joint CDF
    if (jnt) pAll = pAll*adj;

    x = xVal
    p = colMeans( pAll )

    # Create list for output
    output = list(
      # Save plotting values
      pv = as.data.frame( cbind( x = x, y = p ) ),
      # Save grouping factor info
      g = list( w = pAll, x = xVal, n = n[ n > 0 ] ),
      # Save variables for calculations
      v = list( n = n, adj = adj, sel = sel ),
      # Save input
      i = as.data.frame( cbind( rt = rt, ch = ch, grp = grp ) ),
      # Save options
      opt = optOut
    )

  }

  if (draw) {
    if (flip) lines( x, -p, ... ) else lines( x, p, ... )
  }

  if (out) return( output )
}

# Lookup - 02
#' PDF curves for response time and choice data.
#'
#' Draws a line for the estimated PDF of a set response
#' times (conditioned on choice) on an already existing plot.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the PDF for responses times when choice = 0 is
#'   drawn. If 1, the PDF for responses times when choice = 1 is drawn.
#' @param grp an optional vector with a grouping factor (e.g. subjects).
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{jnt}}{If true, the joint distribution is used.}
#'     \item{\code{draw}}{If true, the curve is drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve is flipped about the
#'       x-axis.}
#'   }
#' @param ...  additional plotting parameters.
#' @return A list consisting of...
#' \describe{
#'   \item{\code{pv}}{a data frame with the plotting values
#'     used for the x-axis and the y-axis.}
#'   \item{\code{g}}{when a grouping factor is present, a list
#'     with the matrix of y-axis values per level, the vector
#'     of associated response times, and the total number of
#'     observations per each included level}
#'   \item{\code{v}}{a list of additional variables, the total
#'     number of observations for all levels, the choice proportion
#'     for each level, the choice selection, and when there
#'     is no grouping factor, the response times and their
#'     associated estimated density.}
#'   \item{\code{i}}{a list of the input variables.}
#'   \item{\code{opt}}{a list of the options used.}
#'   }
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 4 & d$Subject == 1
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot( yDim = c(0,4), ver='PDF' )
#' pdf_curve( rt, ch )
#' pdf_curve( rt, ch, sel = 0, lty = 2 )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 0 # Not all subjects had responses for each choice
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( bty = 'l', ver = 'PDF', yDim = c(0,4), cex.axis = 1.5, cex.lab = 1.5 )
#' pdf_curve( rt, ch, grp = grp, lwd = 2 )
#' pdf_curve( rt, ch, sel = 0, grp = grp, lwd = 2, lty = 2 )
#' @export

pdf_curve = function( rt, ch, sel = 1, grp = NULL,
                      opt = list( ), ... ) {

  # Set options for joint distribution, drawing, output, and
  # whether curve should be flipped around x-axis
  if ( length( opt$jnt ) == 0 ) jnt = T else jnt = opt$jnt
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  # Save options
  optOut = list( jnt = jnt, draw = draw,
                 out = out, flip = flip )

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    n = sum( ch == sel ) # Determine number of observations

    if ( n < 2 ) stop('Not enough observations')

    # Estimate density using base R
    xa = sort( rt[ ch == sel ] )
    dn = density( xa )
    x = dn$x; d = dn$y
    # Bound by minimum and maximum RT
    lb = min( rt[ ch == sel ] )
    ub = max(rt[ ch == sel ] )
    kp = x >= lb & x <= ub
    x = x[kp]; d = d[kp]

    # Weight density for joint pdf
    adj = mean( ch == sel )
    if (jnt) d = d*adj

    # Extract density estimates for the set of
    # observed response times
    df = approxfun(dn) # Approximates density function
    vl = df(xa) # Calculates density for observed data
    if (jnt) vl = vl*adj

    # Create list for output
    output = list(
      # Save plotting values
      pv = as.data.frame( cbind( x = x, y = d ) ),
      # Save grouping factor info
      g = NULL,
      # Save variables for calculations
      v = list( n = n, adj = adj, sel = sel, rt = xa, d = vl ),
      # Save input
      i = as.data.frame( cbind( rt = rt, ch = ch ) ),
      # Save options
      opt = optOut
    )

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    # Determine number of observations per level of
    # grouping factor
    n = aggregate( ch == sel, list( grp ), sum )
    # Determine if there are sufficient observations
    cmp = round( mean(n$x) )
    if (cmp < 2) stop('Not enough observations')

    # Condition on choice
    cur_t = rt[ ch == sel ]
    cur_grp = grp[ ch == sel ]

    # Skip subjects without a sufficient number of observations
    rmv = n[ n$x < 2, 1 ]
    if ( length( rmv ) > 0 ) {
      for (i in 1:length(rmv) ) cur_grp = cur_grp[ cur_grp != rmv[i] ]
    }

    # Determine range of response times over which to
    # estimate density
    x = seq( min(cur_t), quantile( cur_t, .99), length = 100 )
    # Levels of the grouping factor
    all_grp = sort( unique( cur_grp ) )

    # Calculate density for each level of grouping factor
    all_dn = matrix( 0, length( all_grp ), 100 )
    for (g in 1:length(all_grp) ) {

      dn = density( cur_t[ cur_grp == all_grp[g] ] )

      # Extract density estimates for each response time
      df = approxfun(dn) # Approximates density function
      all_dn[g,] = df(x);

    }
    all_dn[is.na( all_dn )] = 0

    # Weight density for joint pdf
    adj = aggregate( ch, list( grp ), function(x) mean(x == sel) )$x
    # Remove subjects without sufficient number of observations
    adj = adj[ n$x >= 2 ]
    if (jnt) all_dn = all_dn*adj;

    # Aggregate over subjects
    d = colMeans( all_dn )

    dn = all_dn

    # Create list for output
    output = list(
      # Save plotting values
      pv = as.data.frame( cbind( x = x, y = d ) ),
      # Save grouping factor info
      g = list( w = all_dn, x = x, n = n$x[ n$x >= 2 ] ),
      # Save variables for calculations
      v = list( n = n$x, adj = adj, sel = sel ),
      # Save input
      i = as.data.frame( cbind( rt = rt, ch = ch, grp = grp ) ),
      # Save options
      opt = optOut
    )

  }

  if (draw) {
    if (flip) lines( x, -d, ... ) else lines( x, d, ... )
  }

  if (out) return( output )
}

# Lookup - 03
#' Estimated the hazard function for response time and choice data.
#'
#' Draws a smoothed estimate of the hazard function for a set of
#' response times (conditioned on choice) using an algorithm
#' recommended by Luce (1986, see equations 4.1 and 4.2).
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the hazard function for responses times when choice = 0 is
#'   drawn. If 1, the hazard function for responses times when choice = 1 is
#'   drawn.
#' @param prb the sequence of cumulative probabilities used to define
#'   the intervals.
#' @param j a parameter controlling the degree of smoothing.
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{jnt}}{If true, the joint distribution is used.}
#'     \item{\code{draw}}{If true, the curve is drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve is flipped about the
#'       x-axis.}
#'   }
#' @param ...  additional plotting parameters.
#' @references Luce, R. D. (1986). Response Times: Their Role in Inferring
#'   Elementary Mental Organization. New York: Oxford University Press.
#' @return A list consisting of ...
#' \describe{
#'   \item{\code{pv}}{a data frame with the plotting values
#'     used for the x-axis and the y-axis.}
#'   \item{\code{g}}{when a grouping factor is present, a list
#'     with the matrix of y-axis values per level and the matrix
#'     of quantiles for each level.}
#'   \item{\code{v}}{a list of additional variables, the total
#'     number of observations for all levels, the choice proportion
#'     for each level, the choice selection, the smoothing
#'     parameter j, and the cumulative probabilities used to
#'     define the intervals.}
#'   \item{\code{i}}{a list of the input variables.}
#'   \item{\code{opt}}{a list of the options used.}
#'   }
#' @export
#' @examples
#' #' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 4 & d$Subject == 1
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot( yDim = c(0,5), ver='HF' )
#' hazard_curve( rt, ch )
#' hazard_curve( rt, ch, sel = 0, lty = 2 )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 4
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( bty = 'l', ver = 'HF', yDim = c(0,5), cex.axis = 1.5, cex.lab = 1.5 )
#' hazard_curve( rt, ch, grp = grp, lwd = 2 )
#' hazard_curve( rt, ch, sel = 0, grp = grp, lwd = 2, lty = 2 )

hazard_curve = function( rt, ch, sel = 1, prb = seq(.05,.95,.1),
                         j = 25, grp = NULL, opt = list( ), ... ) {

  # Set options for joint distribution, drawing, output, and
  # whether curve should be flipped around x-axis
  if ( length( opt$jnt ) == 0 ) jnt = F else jnt = opt$jnt
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  # Save options
  optOut = list( jnt = jnt, draw = draw,
                 out = out, flip = flip )

  # Equation 4.2 from Luce (1986)
  S = function( n, k, Z ) {

    Z_k = Z[k]
    if (k == 1) Z_km1 = 0 else Z_km1 = Z[k-1]

    out = ( n - k + 1 ) * ( Z_k - Z_km1 )

    return( out )
  }

  # Equation 4.1 from Luce (1986)
  l_hat = function(i,j,n,Z) {

    beg = i - j + 1; if (beg < 1) { beg = 1; j = i }
    denom = numeric( i - beg )

    inc = 1;
    for ( k in beg:i ) { denom[inc] = S( n, k, Z ); inc = inc + 1 }
    out = j/sum( denom )

    return( out )
  }

  # Function to determine total number of observations
  # less than or equal to a particular value
  f = function(q) return( sum( Z <= q ) )

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    n = sum( ch == sel ) # Determine number of observations
    if ( n < length(prb) ) stop('Not enough observations')

    # Condition on choice
    x = rt[ ch == sel ]

    # Determine number of intervals
    Z = sort( x ) # Sort observations
    q = quantile( x, prob = prb )
    int = sapply( q, f )

    # Estimate a smoothed version of the hazard function
    h = sapply( int, l_hat, j = j, n = n, Z = Z )
    x = Z[ int ];

    # Adjustment for joint distributions
    adj = mean( ch == sel )
    if (jnt) h = h*adj

    # Create list for output
    output = list(
      # Save plotting values
      pv = as.data.frame( cbind( x = x, y = h ) ),
      # Save grouping factor info
      g = NULL,
      # Save variables for calculations
      v = list( n = n, adj = adj, sel = sel, j = j, prb = prb ),
      # Save input
      i = as.data.frame( cbind( rt = rt, ch = ch ) ),
      # Save options
      opt = optOut
    )

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    # Determine number of observations per level of
    # grouping factor
    n = aggregate( ch == sel, list( grp ), sum )
    # Adjust estimate of curve based on
    # average number of observations
    cmp = round( mean(n$x) )
    if (cmp < 2) stop('Not enough observations')

    # Condition on choice
    xAll = rt[ ch == sel ]
    gAll = grp[ ch == sel ]
    nAll = n$x[ n$x > 0 ]

    # Calculate group quantiles
    qAll = aggregate( xAll, list( gAll ), quantile, prob = prb )
    colnames( qAll ) = c('G','Q')
    q = colMeans( qAll$Q )

    hAll = matrix( 0, nrow( qAll ), ncol( qAll$Q ) )

    for ( i in 1:nrow( qAll ) ) {

      # Determine number of intervals
      Z = sort( xAll[ gAll == qAll$G[i] ] ) # Sort observations
      int = sapply( q, f )

      ind = min( which( int > 0 ) ):length(q)

      # Estimate a smoothed version of the hazard function
      hAll[i,ind] = sapply( int[ind], l_hat, j = j, n = nAll[i], Z = Z )
      print( hAll[i,] )

    }

    h = colMeans( hAll )
    x = q

    # Adjustment for joint distributions
    adj = mean( ch == sel )
    if (jnt) h = h*adj

    # Create list for output
    output = list(
      # Save plotting values
      pv = as.data.frame( cbind( x = x, y = h ) ),
      # Save grouping factor info
      g = list( w = hAll, x = qAll ),
      # Save variables for calculations
      v = list( n = n, adj = adj, sel = sel, j = j, prb = prb ),
      # Save input
      i = as.data.frame( cbind( rt = rt, ch = ch, grp = grp ) ),
      # Save options
      opt = optOut
    )

  }

  if (draw) {
    if (flip) lines( x, -h, ... ) else lines( x, h, ... )
  }

  if (out) return( output )
}

# Lookup - 04
#' Quantile-probability estimates
#'
#' Draws the estimated quantiles based on the corresponding
#' cumulative probabilities for a set of response times
#' (conditioned on choice).
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the quantile-probability estimates for responses
#'   times when choice = 0 are drawn. If 1, the quantile-probability
#'   estimates for responses times when choice = 1 are drawn.
#' @param prb the sequence of cumulative probabilities for which the
#'   quantiles should be determined.
#' @param grp an optional vector with a grouping factor (e.g. subjects).
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{jnt}}{If true, the joint distribution is used.}
#'     \item{\code{draw}}{If true, the curve is drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve is flipped about the
#'       x-axis.}
#'     \item{\code{pts}}{If true, draw points instead of lines.}
#'   }
#' @param ...  additional plotting parameters.
#' @return A list consisting of...
#' \describe{
#'   \item{\code{pv}}{a data frame with the plotting values
#'     used for the x-axis and the y-axis.}
#'   \item{\code{g}}{when a grouping factor is present, a list
#'     with the matrix of y-axis values per level and the
#'     associated matrix of quantiles.}
#'   \item{\code{v}}{a list of additional variables, the choice
#'     proportion for each level and the choice selection.}
#'   \item{\code{i}}{a list of the input variables.}
#'   \item{\code{opt}}{a list of the options used.}
#'   }
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 6 & d$Subject == 20
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot( xDim = c(0,1) )
#' quantile_points( rt, ch, opt = list( pts = F ) )
#' quantile_points( rt, ch, pch = 21, bg = 'white' )
#' quantile_points( rt, ch, sel = 0, opt = list( pts = F ), lty = 2 )
#' quantile_points( rt, ch, sel = 0, pch = 24, bg = 'white' )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 6
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( xDim = c(0,1), bty = 'l', cex.axis = 1.5, cex.lab = 1.5 )
#' quantile_points( rt, ch, grp = grp, opt = list( pts = F ), lwd = 2 )
#' quantile_points( rt, ch, grp = grp, pch = 21, bg = 'white', cex = 1.2 )
#' quantile_points( rt, ch, sel = 0, grp = grp, opt = list( pts = F ), lwd = 2, lty = 2 )
#' quantile_points( rt, ch, sel = 0, grp = grp, pch = 24, bg = 'white', cex = 1.2 )
#' @export

quantile_points = function( rt, ch, sel = 1,
                            prb = seq( .1, .9, .2 ), grp = NULL,
                            opt = list( ), ... ) {

  # Set options for joint distribution, drawing, output, and
  # whether curve should be flipped around x-axis
  if ( length( opt$jnt ) == 0 ) jnt = T else jnt = opt$jnt
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  if ( length( opt$pts ) == 0 ) pts = T  else pts = opt$pts
  # Save options
  optOut = list( jnt = jnt, draw = draw,
                 out = out, flip = flip, pts = pts )

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    # Total number of observations
    n = sum( ch == sel )
    if ( n < 1 ) stop('No observations')

    # Condition on choice
    x = rt[ ch == sel ]

    # Calculate adjustment for joint distribution
    adj = mean( ch == sel )

    # Calculate quantiles
    q = quantile( x, prob = prb )
    if (jnt) p = prb*adj else p = prb

    output = list( pv = cbind( x = q, y = p ),
                   g = NULL,
                   v = list(adj = adj, sel = sel),
                   i = list( rt = rt, ch = ch ),
                   opt = optOut )

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    # Determine number of observations per level of
    # grouping factor
    n = aggregate( ch == sel, list( grp ), sum )$x
    # Determine if there are sufficient observations
    cmp = round( mean(n) )
    if (cmp < 1) stop('Not enough observations')

    # Condition on choice
    xAll = rt[ ch == sel ]
    grpAll = grp[ ch == sel ]

    # Calculate quantiles over levels of grouping factor
    allQ = aggregate( xAll, list( grpAll ), quantile, prob = prb )
    colnames( allQ ) = c('G','Q')

    # Calculate adjustment for joint distribution over
    # levels of grouping factor
    adj = aggregate( ch, list( grp ), function(x) mean(x == sel) )
    y = matrix( prb, nrow(allQ), length(prb), byrow = T )
    if (jnt) y = y*adj$x;

    q = colMeans( allQ$Q )
    p = colMeans( y )

    output = list( pv = cbind( x = q, y = p ),
                   g = list( w = y, x = allQ$Q ),
                   v = list(adj = adj$x, sel = sel),
                   i = list( rt = rt, ch = ch, grp = grp ),
                   opt = optOut )
  }

  if (draw) {
    if (pts) {
      if (flip) points( q, -p, ... ) else points( q, p, ... )
    } else {
      if (flip) lines( q, -p, ... ) else lines( q, p, ... )
    }
  }

  if (out) return( output )
}

# Lookup - 05
#' Conditional accuracy function
#'
#' Draws the conditional accuracy at different response time quantiles.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param prb the sequence of cumulative probabilities for which the
#'   quantiles should be determined.
#' @param grp an optional vector with a grouping factor (e.g. subjects).
#' @param opt a list of named options:
#'   \describe{
#'     \item{\code{jnt}}{If true, the joint distribution is used.}
#'     \item{\code{draw}}{If true, the curve is drawn.}
#'     \item{\code{out}}{If true, output is returned.}
#'     \item{\code{flip}}{If true, the curve is flipped about the
#'       x-axis.}
#'     \item{\code{pts}}{If true, draw points instead of lines.}
#'   }
#' @param ...  additional plotting parameters.
#' @return A list consisting of...
#' \describe{
#'   \item{\code{pv}}{a data frame with the plotting values
#'     used for the x-axis and the y-axis.}
#'   \item{\code{g}}{when a grouping factor is present, a list
#'     with the matrix of y-axis values per level and the
#'     associated matrix of quantiles.}
#'   \item{\code{v}}{a list of additional variables, the choice
#'     proportion for each level and the choice selection.}
#'   \item{\code{i}}{a list of the input variables.}
#'   \item{\code{opt}}{a list of the options used.}
#'   }
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 12 & d$Subject == 7
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot( xDim = c(0,1), ver = 'CAF' )
#' caf_points( rt, ch, opt = list( pts = F ) )
#' caf_points( rt, ch, pch = 21, bg = 'white' )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 12
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( xDim = c(0,1.2), ver = 'CAF', cex.axis = 1.5, cex.lab = 1.5 )
#' caf_points( rt, ch, grp = grp, opt = list( pts = F ), lwd = 2 )
#' caf_points( rt, ch, grp = grp, pch = 21, bg = 'white', cex = 1.2 )
#' @export

caf_points = function( rt, ch, prb = seq( .1, .9, .2 ),
                       grp = NULL, opt = list( ), ... ) {

  # Set options for joint distribution, drawing, output, and
  # whether curve should be flipped around x-axis
  if ( length( opt$jnt ) == 0 ) jnt = T else jnt = opt$jnt
  if ( length( opt$draw ) == 0 ) draw = T else draw = opt$draw
  if ( length( opt$out ) == 0 ) out = F else out = opt$out
  if ( length( opt$flip ) == 0 ) flip = F else flip = opt$flip
  if ( length( opt$pts ) == 0 ) pts = T  else pts = opt$pts
  # Save options
  optOut = list( jnt = jnt, draw = draw,
                 out = out, flip = flip, pts = pts )

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    q = quantile( rt, prob = prb )
    ca = sapply( q, function(q) sum( ch[rt < q] )/sum( rt < q ) )

    output = list( pv = cbind( x = q, y = ca ),
                   g = NULL,
                   v = list( sel = sel ),
                   i = list( rt = rt, ch = ch ),
                   opt = optOut )

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    allQ = aggregate( rt, list( grp ), quantile, prob = prb )
    colnames( allQ ) = c( 'G', 'Q' )

    all_grp = sort( unique( grp ) )
    allCAF = matrix( NA, nrow( allQ ), ncol( allQ$Q ) )
    for (i in 1:nrow(allQ)) {
      chg = ch[ grp == all_grp[i] ]
      ct = rt[ grp == all_grp[i] ]
      allCAF[i,] = sapply( allQ$Q[i,],
                           function(q) {
                             sum( chg[ct < q] )/sum( ct < q )
                           } )
    }

    q = colMeans( allQ$Q )
    ca = colMeans( allCAF )

    output = list( pv = cbind( x = q, y = ca ),
                   g = list( w = allCAF, x = allQ$Q ),
                   v = list(sel = sel),
                   i = list( rt = rt, ch = ch, grp = grp ),
                   opt = optOut )

  }

  if (draw) {
    if (pts) {
      if (flip) points( q, -ca, ... ) else points( q, ca, ... )
    } else {
      if (flip) lines( q, -ca, ... ) else lines( q, ca, ... )
    }
  }

  if (out) return( output )
}

# Lookup - 06
#' Add additional points.
#'
#' Adds additional points to an existing plot, drawing a given
#' test statistic for response times conditioned on choice using
#' output from a curve function.
#'
#' @param output the list output from \code{cdf_curve},
#'   \code{pdf_curve}, or \code{hazard_curve}.
#' @param T_x a function to calculate a text statistic over
#'   a vector (e.g. mean, median, etc.).
#' @param out A logical value, indicating if output should be returned.
#' @param ...  additional plotting parameters.
#' @return A list consisting of the x and y-axis plotting points.
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 15 & d$Subject == 1
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot()
#' out = cdf_curve( rt, ch, opt = list( out = T ) )
#' add_points( out, pch = 19, col = 'blue' )
#' add_points( out, T_x = median, pch = 19, col = 'red' )
#' legend( 'topleft', c('Mean','Median'), fill = c('blue','red'), bty = 'n' )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 15
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot()
#' out = cdf_curve( rt, ch, opt = list( out = T ) )
#' T_x = function(x) quantile(x,prob=seq(.1,.9,.1))
#' add_points( out, T_x = T_x, pch = 21, bg = 'white' )
#' @export

add_points = function( output, T_x = mean, out = F, ... ) {

  # Extract choice to condition on
  sel = output$v$sel
  # Extract input
  rt = output$i$rt; ch = output$i$ch
  # Extract options on whether to flip
  flip = output$opt$flip

  # Determine if there is a grouping variable
  grp = output$i$grp

  # No grouping factor
  if ( length( grp ) == 0 ) {

    # Condition on choice
    x = rt[ ch == sel ]
    # Calculate statistic
    ts = T_x( x )

  }
  # Grouping factor
  if ( length( grp ) == length( rt ) ) {

    # Condition on choice
    xAll = rt[ ch == sel ]
    g = grp[ ch == sel ]

    # Calculate statistic
    tsAll = aggregate( xAll, list( g ), T_x )

    # Collapse over group levels
    if ( is.matrix( tsAll$x ) ) {
      ts = colMeans( tsAll$x )
    } else {
      ts = mean( tsAll$x )
    }

  }

  # Extract plotting values
  xa = output$pv$x # x-axis
  ya = output$pv$y # y-axis

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
        out = lnInterp( x, ya[ pts ], xa[ pts ] )
      }

    }

    return( out )
  }

  # Determine y-axis values
  ny = sapply( ts, f )

  # Save output
  new_output = list( x = ts, y = ny )

  # Add points
  if (flip) points( ts, -ny, ... ) else points( ts, ny, ... )

  if (out) return( new_output )
}


# Lookup - 07
#' Blank response time plot
#'
#' Creates a blank response time plot with standard labels.
#'
#' @param xDim the minimum and maximum x-axis values.
#' @param yDim the minimum and maximum y-axis values.
#' @param ver the type of plot to draw (i.e. 'CDF', 'PDF',
#'   'HF', 'CAF', or 'blank' ).
#' @param unit the unit for the response times (e.g. 'ms' or 's').
#' @param ... additional plotting parameters.
#' @return A blank response time plot.
#' @examples
#' layout( rbind( 1:3, 4:6 ) )
#' blankRTplot(); blankRTplot(ver='PDF');
#' blankRTplot(ver='CAF'); blankRTplot(ver='HF');
#' blankRTplot(ver='blank')
#' @export

blankRTplot = function( xDim = c(0,2), yDim = c(0,1),
                        ver = 'CDF', unit = 's', ... ) {

  # Define label for x-axis
  xl = paste( 'RT (', unit, ')', sep = '' )

  if ( ver == 'CDF' | ver == 'QPE' ) yl = 'Cumulative probability'
  if ( ver == 'PDF' ) yl = 'Density'
  if ( ver == 'CAF' ) yl = 'Conditional accuracy'
  if ( ver == 'HF' ) yl = 'Hazard function'
  if ( ver == 'blank' ) { xl = ' '; yl = ' ' }

  plot( xDim, yDim, type = 'n', ylab = yl, xlab = xl, ... )

}
