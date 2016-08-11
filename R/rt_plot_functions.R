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
# Lookup - 06:  blankRTplot

### TO DO ###
# Add references
# Add examples (hazard function)
# Add aggregation over group factor (hazard function)
# Add accuracy-latency plots
# Add helper functions for plotting variability (e.g. SEs)

# Useful functions for package creation
# library(devtools)
# library(roxygen2)

# Lookup - 01
#' Joint CDF curves for response time and choice data.
#'
#' Draws a line for the estimated joint CDF of a set response
#' times and choices on an already existing plot.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the CDF for responses times when choice = 0 is
#'   drawn. If 1, the CDF for responses times when choice = 1 is drawn.
#' @param grp an optional vector with a grouping factor (e.g. subjects).
#' @param plt a named list that allows specification of graphical
#'   parameters. The named list plt$ln can be used to set the graphical
#'   elements (e.g.. the options \code{lwd} and \code{lty}) for the
#'   drawn line. The named lists plt$pt1 and plt$pt2 can be used to
#'   set the graphical elements (e.g. the option \code{pch}) for the
#'   points denoting the median and mean.
#' @param opt logical vector; indicates if 1) the joint distribution
#'   should be used, 2) the line should be drawn, and 3) if output
#'   should be returned.
#' @return A list consisting of...
#' \describe{
#'   \item{\code{CDF}}{a matrix with the estimated cumulative
#'     probabilities and the corresonding resopnse times.}
#'   \item{\code{Median}}{the response time and cumulative
#'     probability for the median}
#'   \item{\code{Mean}}{the response time and cumulative
#'     probability for the mean}
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
#' plt = list( ln = list( lty = 2 ) )
#' cdf_curve( rt, ch, sel = 0, plt = plt )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 4
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( bty = 'l', ver = 'CDF', cex.axis = 1.5, cex.lab = 1.5 )
#' plt = list( ln = list( lwd = 2 ), pt1 = list( bg = 'grey' ) )
#' cdf_curve( rt, ch, grp = grp, plt = plt )
#' plt = list( ln = list( lty = 2, lwd = 2 ), pt1 = list( bg = 'grey' ) )
#' cdf_curve( rt, ch, sel = 0, grp = grp, plt = plt )
#' @export

cdf_curve = function( rt, ch, sel = 1, grp = NULL, plt = NULL,
                      opt = c(T,T,F) ) {

  # Set options for joint distribution, drawing, and output
  jnt = opt[1];
  draw = opt[2];
  out = opt[3];

  # Determine line and point characteristics
  tmp = curve_defaults( plt )
  lnDefaults = tmp[[1]]
  medDefaults = tmp[[2]]
  avgDefaults = tmp[[3]]

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    # Empirical estimate of CDF for response times
    n = sum( ch == sel )
    p = (1:n)/n
    x = sort( rt[ ch == sel ] )

    # Determine median
    medVal = max( which( p <= .5 ) )
    xMed = x[ medVal ];

    # Determine mean
    xAvg = mean( x )
    intrvl = c( max( which( x <= xAvg ) ), min( which( x > xAvg ) ) )

    # Adjust asymptote if estimating joint CDF
    if (jnt) p = p*(n/length(ch));

    # Determine y-axis values
    pMed = p[ medVal ];
    pAvg = lnInterp( xAvg, p[intrvl], x[intrvl] )
  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    xAll = rt[ ch == sel ]
    g = grp[ ch == sel ]

    # Adjust estimate of curve based on sample size
    n = aggregate( xAll, list( g ), length )$x
    if ( min(n) < 20 ) prb = (1:round( mean( n ) ))/round( mean( n ) )
    if ( min(n) >= 20 & min(n) < 100 ) prb = seq( 0, 1, .05 )
    if ( min(n) >= 100 & min(n) < 200 ) prb = seq( 0, 1, .02 )
    if ( min(n) >= 200 ) prb = seq( 0, 1, .01 )

    allQ = aggregate( xAll, list( g ), quantile,
                   prob = prb )
    colnames( allQ ) = c('G','Q')

    adj = aggregate( ch, list( grp ), function(x) mean(x == sel) )
    p = matrix( prb, nrow(allQ), length(prb), byrow = T )
    if (jnt) p = p*adj$x;

    if ( length( dim(allQ$Q) ) > 0 ) x = colMeans( allQ$Q ) else
      x = mean( allQ$Q )
    p = colMeans( p )

    xAvg = mean( aggregate( xAll, list( g ), mean )$x )
    intrvl = c( max( which( x <= xAvg ) ), min( which( x > xAvg ) ) )
    pAvg = lnInterp( xAvg, p[intrvl], x[intrvl] )

    xMed = mean( aggregate( xAll, list( g ), median )$x )
    intrvl = c( max( which( x <= xMed ) ), min( which( x > xMed ) ) )
    pMed = lnInterp( xMed, p[intrvl], x[intrvl] )

  }

  if (draw) {
    lines( x, p, lty = lnDefaults$lty, lwd = lnDefaults$lwd,
           col = lnDefaults$col, type = lnDefaults$type )
    points( xMed, pMed, pch = medDefaults$pch,
            bg = medDefaults$bg, lwd = medDefaults$lwd,
            col = medDefaults$col, cex = medDefaults$cex,
            type = medDefaults$type )
    points( xAvg, pAvg, pch = avgDefaults$pch,
            bg = avgDefaults$bg, lwd = avgDefaults$lwd,
            col = avgDefaults$col, cex = avgDefaults$cex,
            type = avgDefaults$type )
  }

  if (out) return( list( CDF = cbind(x,p), Median = c( xMed, pMed ),
                         Mean = c( xAvg, pAvg ) ) )
}

# Lookup - 02
#' Joint PDF curves for response time and choice data.
#'
#' Draws a line for the estimated joint PDF of a set response
#' times and choices on an already existing plot.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the PDF for responses times when choice = 0 is
#'   drawn. If 1, the PDF for responses times when choice = 1 is drawn.
#' @param grp an optional vector with a grouping factor (e.g. subjects).
#' @param plt a named list that allows specification of graphical
#'   parameters. The named list plt$ln can be used to set the graphical
#'   elements (e.g.. the options \code{lwd} and \code{lty}) for the
#'   drawn line. The named lists plt$pt1 and plt$pt2 can be used to
#'   set the graphical elements (e.g. the option \code{pch}) for the
#'   points denoting the median and mean.
#' @param opt logical vector; indicates if 1) the joint distribution
#'   should be used, 2) the line should be drawn, and 3) if output
#'   should be returned.
#' @return A list consisting of...
#' \describe{
#'   \item{\code{PDF}}{a matrix with the estimated density}
#'   \item{\code{Mode}}{the response time and likelihood for the mode}
#'   \item{\code{Mean}}{the response time and likelihood for the mean}
#'   }
#' @examples
#' # Load in example data
#' data(priming_data)
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 5 & d$Subject == 16
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot(yDim=c(0,3),ver='PDF')
#' pdf_curve( rt, ch )
#' plt = list( ln = list( col = 'blue' ) )
#' pdf_curve( rt, ch, sel = 0, plt = plt )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 5
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( yDim = c(0,2), bty = 'l', ver = 'CDF', cex.axis = 1.5,
#'   cex.lab = 1.5 )
#' plt = list( ln = list( lwd = 2 ), pt1 = list( bg = 'grey', cex = 1.5 ),
#'   pt2 = list( type = 'n' ) )
#' pdf_curve( rt, ch, grp = grp, plt = plt )
#' plt = list( ln = list( lty = 2, lwd = 2 ),
#'   pt1 = list( pch = 15, cex = 1.5 ),
#'   pt2 = list( type = 'n' ) )
#' pdf_curve( rt, ch, sel = 0, grp = grp, plt = plt )
#' @export

pdf_curve = function( rt, ch, sel = 1, grp = NULL, plt = NULL,
                      opt = c(T,T,F) ) {

  # Set options for joint distribution, drawing, and output
  jnt = opt[1];
  draw = opt[2];
  out = opt[3];

  # Determine line and point characteristics
  tmp = curve_defaults( plt )
  lnDefaults = tmp[[1]]
  modeDefaults = tmp[[2]]
  avgDefaults = tmp[[3]]

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    # Estimate density using base R
    xa = rt[ ch == sel ]
    dn = density( xa )
    x = dn$x; d = dn$y

    if (jnt) d = d*mean( ch == sel )

    # Determine mode
    modeVal = which( d == max( d ) )
    xMode = x[ modeVal ];
    dMode = d[ modeVal ];

    # Determine mean
    xAvg = mean( xa )
    avgVal = max( which( x <= xAvg ) )
    xAvg = x[ avgVal ];
    dAvg = d[ avgVal ];

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    cur_t = rt[ ch == sel ]
    cur_grp = grp[ ch == sel ]

    x = seq( min(cur_t), quantile( cur_t, .99), length = 100 )
    all_grp = sort( unique( cur_grp ) )
    all_dn = matrix( 0, length( all_grp ), 100 )

    for (g in 1:length(all_grp) ) {
      dn = density( cur_t[ cur_grp == all_grp[g] ] )
      # Extract density estimates for
      df = approxfun(dn) # Approximates density function

      all_dn[g,] = df(x);

    }
    all_dn[is.na( all_dn )] = 0

    adj = aggregate( ch, list( grp ), function(x) mean(x == sel) )
    if (jnt) all_dn = all_dn*adj$x;

    # Aggregate over subjects
    d = colMeans( all_dn )

    # Determine mode
    modeVal = which( d == max( d ) )
    xMode = x[ modeVal ];
    dMode = d[ modeVal ];

    # Determine mean
    xAvg = mean( cur_t )
    avgVal = max( which( x <= xAvg ) )
    xAvg = x[ avgVal ];
    dAvg = d[ avgVal ];

    dn = all_dn

  }

  if (draw) {
    lines( x, d, lty = lnDefaults$lty, lwd = lnDefaults$lwd,
           col = lnDefaults$col, type = lnDefaults$type )
    points( xMode, dMode, pch = modeDefaults$pch,
            bg = modeDefaults$bg, lwd = modeDefaults$lwd,
            col = modeDefaults$col, cex = modeDefaults$cex,
            type = modeDefaults$type )
    points( xAvg, dAvg, pch = avgDefaults$pch,
            bg = avgDefaults$bg, lwd = avgDefaults$lwd,
            col = avgDefaults$col, cex = avgDefaults$cex,
            type = avgDefaults$type )
  }

  if (out) return( list(
    PDF = dn, Mode = c( xMode, dMode ), Mean = c( xAvg, dAvg ) ) )
}

# Lookup - 03
#' Estimated hazard function for response time and choice data.
#'
#' Draws a line for the estimated joint PDF of a set response
#' times and choices on an already existing plot.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the PDF for responses times when choice = 0 is
#'   drawn. If 1, the PDF for responses times when choice = 1 is drawn.
#' @param opt logical vector; indicates if 1) the joint distribution
#'   should be used, 2) the line should be drawn, and 3) if output
#'   should be returned.
#' @param ...  additional plotting parameters.
#' @return A matrix giving the selected response times, density
#'   estimates, cumulative probabilities, and the estimated hazard
#'   function.
#' @export

hazard_curve = function(rt, ch, sel = 1, opt = c(T,T,F), ...) {

  # Set options for joint distribution, drawing, and output
  jnt = opt[1];
  draw = opt[2];
  out = opt[3];

  # Extract relevant response times
  x = rt[ ch == sel ]

  adj = mean( ch == sel )

  # Estimate PDF
  dn = density(x)
  # Extract density estimates for observed response times
  df = approxfun(dn) # Approximates density function
  g = df( sort(x) )

  # Estimate distribution function values
  G = (1:length(x))/length(x)

  # Adjust values for joint distribution
  if (jnt) { g = g*adj; G = G*adj }

  h = g/(1-G) # Calculate empirical hazard function

  if (draw) {

    lines( sort(x), h, ... )

  }

  if (out) return( as.data.frame( cbind( rt = sort(x),
                                         g = g, G = G, h = h ) ) )
}

# Lookup - 04
#' Quantile-probability estimates
#'
#' Draws the estimated quantiles based on the corresponding joint
#' cumulative probabilities for a set of response times and choices.
#'
#' @param rt vector of response times.
#' @param ch a vector of binary choices (i.e. 0 or 1).
#' @param sel If 0, the quantile-probability estimates for responses
#'   times when choice = 0 are drawn. If 1, the quantile-probability
#'   estimates for responses times when choice = 1 are drawn.
#' @param prb the sequence of cumulative probabilities for which the
#'   quantiles should be determined.
#' @param plt a named list that allows specification of graphical
#'   parameters. The named list plt$ln can be used to set the graphical
#'   elements (e.g.. the options \code{lwd} and \code{lty}) for the
#'   drawn line. The named list plt$pt can be used to
#'   set the graphical elements (e.g. the option \code{pch}) for the
#'   drawn points.
#' @param opt logical vector; indicates if 1) the joint distribution
#'   should be used, 2) the line should be drawn, and 3) if output
#'   should be returned.
#' @return A matrix giving the quantile estimates and the corresponding
#'   cumulative probabilities.
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 6 & d$Subject == 20
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot(xDim = c(0,1), ver='QPE')
#' quantile_points( rt, ch )
#' plt = list( ln = list( type = 'n' ), pt = list( pch = 18,
#'   col = 'red' ) )
#' quantile_points( rt, ch, sel = 0, plt = plt )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 6
#' rt = d$RT[sel]; ch = d$Choice[sel]; grp = d$Subject[sel]
#' blankRTplot( xDim = c(0,1), bty = 'l', ver = 'QPE', cex.axis = 1.5,
#'   cex.lab = 1.5 )
#' plt = list( ln = list( type = 'n' ),
#'   pt = list( bg = 'grey', pch = 21, cex = 1.5 ) )
#' quantile_points( rt, ch, sel = 1, plt = plt )
#' plt = list( ln = list( type = 'n' ),
#'   pt = list( bg = 'white', pch = 22, cex = 1.5 ) )
#' quantile_points( rt, ch, sel = 0, plt = plt )
#' @export

quantile_points = function( rt, ch, sel = 1, prb = seq( .1, .9, .2 ),
                            grp = NULL, plt = NULL, opt = c(T,T,F) ) {

  # Set options for joint distribution, drawing, and output
  jnt = opt[1];
  draw = opt[2];
  out = opt[3];

  # Determine line and point characteristics
  tmp = point_defaults( plt )
  lnDefaults = tmp[[1]]
  ptDefaults = tmp[[2]]

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    x = rt[ ch == sel ]
    adj = mean( ch == sel )

    q = quantile( x, prob = prb )
    if (jnt) y = prb*adj else y = prb

  }
  # If there is a grouping variable
  if ( length( grp ) == length( rt ) ) {

    xAll = rt[ ch == sel ]
    grpAll = grp[ ch == sel ]

    allQ = aggregate( xAll, list( grpAll ), quantile, prob = prb )
    colnames( allQ ) = c('G','Q')

    adj = aggregate( ch, list( grp ), function(x) mean(x == sel) )
    y = matrix( prb, nrow(allQ), length(prb), byrow = T )
    if (jnt) y = y*adj$x;

    q = colMeans( allQ )
    y = colMeans( y )

  }

  if (draw) {
    lines( q, y, lty = lnDefaults$lty, lwd = lnDefaults$lwd,
           col = lnDefaults$col, type = lnDefaults$type )
    points( q, y, pch = ptDefaults$pch,
            bg = ptDefaults$bg, lwd = ptDefaults$lwd,
            col = ptDefaults$col, cex = ptDefaults$cex,
            type = ptDefaults$type )
  }

  if (out) return( cbind( q = q, p = y ) )
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
#' @param plt a named list that allows specification of graphical
#'   parameters. The named list plt$ln can be used to set the graphical
#'   elements (e.g.. the options \code{lwd} and \code{lty}) for the
#'   drawn line. The named list plt$pt can be used to
#'   set the graphical elements (e.g. the option \code{pch}) for the
#'   drawn points.
#' @param opt logical vector; indicates if 1) the line should be drawn,
#'   and 2) if output should be returned.
#' @return A matrix giving the quantile estimates and the corresponding
#'   cumulative probabilities.
#' @examples
#' # Load in example dataset
#' data("priming_data")
#' d = priming_data
#' layout( cbind(1,2) )
#' # Single subject
#' sel = d$Condition == 2 & d$Subject == 6
#' rt = d$RT[sel]; ch = d$Accuracy[sel]
#' blankRTplot( xDim = c(.2,1), yDim=c(.5,1), ver='CAF')
#' plt = list( pt = list( pch = 19 ) )
#' CAF_points( rt, ch, plt = plt )
#' # Aggregating over multiple subjects
#' sel = d$Condition == 2
#' rt = d$RT[sel]; ch = d$Accuracy[sel]; grp = d$Subject[sel]
#' blankRTplot( xDim = c(.2,1), yDim = c(.5,1),
#'   bty = 'l', ver = 'QPE', cex.axis = 1.5,
#'   cex.lab = 1.5 )
#' plt = list( pt = list( bg = 'grey', pch = 23, cex = 1.2 ) )
#' CAF_points( rt, ch, grp = grp, plt = plt )
#' @export

CAF_points = function( rt, ch, prb = seq( .1, .9, .2 ),
                            grp = NULL, plt = NULL, opt = c(T,F) ) {

  # Set options for joint distribution, drawing, and output
  draw = opt[1];
  out = opt[2];

  # Determine line and point characteristics
  tmp = point_defaults( plt )
  lnDefaults = tmp[[1]]
  ptDefaults = tmp[[2]]

  # If there is no grouping variable
  if ( length( grp ) == 0 ) {

    q = quantile( rt, prob = prb )
    CAF = sapply( q, function(q) sum( ch[rt < q] )/sum( rt < q ) )

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
    CAF = colMeans( allCAF )

  }

  if (draw) {
    lines( q, CAF, lty = lnDefaults$lty, lwd = lnDefaults$lwd,
           col = lnDefaults$col, type = lnDefaults$type )
    points( q, CAF, pch = ptDefaults$pch,
            bg = ptDefaults$bg, lwd = ptDefaults$lwd,
            col = ptDefaults$col, cex = ptDefaults$cex,
            type = ptDefaults$type )
  }

  if (out) return( cbind( prb = prb, q = q, CAF = CAF ) )
}

# Lookup - 06
#' Blank response time plot
#'
#' Creates a blank response time plot with standard labels.
#'
#' @param xDim the minimum and maximum x-axis values.
#' @param yDim the minimum and maximum y-axis values.
#' @param ver the type of plot to draw (e.g. 'CDF', 'PDF',
#'   'QPE', or 'CAF').
#' @param unit the unit for the response times (e.g. 'ms' or 's').
#' @return A blank response time plot.
#' @export

blankRTplot = function( xDim = c(0,2), yDim = c(0,1),
                        ver = 'CDF', unit = 's', ... ) {

  if (ver == 'CDF') {
    plot( xDim, yDim, type = 'n',
          ylab = 'Cumulative probability',
          xlab = paste('RT (',unit,')',sep=''), ... )
  }
  if (ver == 'PDF') {
    plot( xDim, yDim, type = 'n',
          ylab = 'Density',
          xlab = paste('RT (',unit,')',sep=''), ... )
  }
  if (ver == 'QPE') {
    plot( xDim, yDim, type = 'n',
          ylab = 'Cumulative probability',
          xlab = paste('RT (',unit,')',sep=''), ... )
  }
  if (ver == 'CAF') {
    plot( xDim, yDim, type = 'n',
          ylab = 'Conditional accuracy',
          xlab = paste('RT (',unit,')',sep=''), ... )
  }

}
